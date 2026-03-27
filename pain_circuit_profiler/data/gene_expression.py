"""
data/gene_expression.py — Allen Brain Atlas ISH gene expression interface
=========================================================================
Queries the Allen Mouse Brain Atlas (ISH) for gene expression per brain region.

DATA SOURCE:
    Allen Mouse Brain Atlas — in situ hybridization (ISH) dataset
    Lein et al., Nature 2007 | https://mouse.brain-map.org/

KEY CONCEPT:
    For each gene, the Allen atlas contains ISH experiments (one or more
    section data sets). Each section data set is analyzed to produce
    expression_energy per brain region:

        expression_energy = (expressing pixel intensity sum) / (total voxel volume)

    This is comparable across regions and experiments, making it suitable
    for ranking regions by relative expression level.

API USED:
    Allen Brain Atlas RMA (RESTful Model Access) API
    https://api.brain-map.org/api/v2/data/query.json
    We use requests directly rather than allensdk wrappers for this module
    because allensdk does not expose a high-level ISH cache interface.

MAIN ENTRY POINT:
    loader = GeneExpressionLoader(cache_dir='allen_connectivity_cache')
    expr_df = loader.load(genes=['Oprm1', 'Pdyn'], structure_ids=[12345, 67890])
"""

import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)

# Allen Brain Atlas RMA endpoint
_RMA_URL = "https://api.brain-map.org/api/v2/data/query.json"

# Product abbreviation for adult mouse ISH (Allen Mouse Brain Atlas)
_ISH_PRODUCT = "Mouse"

# Maximum rows per single API response page
_PAGE_SIZE = 2000


def _rma_query(criteria: str, only: list = None, include: str = None,
               num_rows: int = _PAGE_SIZE, start_row: int = 0) -> list:
    """
    Execute a single Allen RMA query and return the result list.

    Parameters
    ----------
    criteria : str
        RMA model + criteria string, e.g.
        "model::Gene,rma::criteria,[acronym$eq'Oprm1']"
    only : list of str, optional
        Fields to return (reduces payload size).
    include : str, optional
        Related models to include (e.g. 'genes').
    num_rows : int
        Page size.
    start_row : int
        Pagination offset.

    Returns
    -------
    list of dict
        API result objects.

    Raises
    ------
    requests.HTTPError
        If the API returns a non-200 status code.
    ValueError
        If the API returns success=false.
    """
    params = {
        "criteria":  criteria,
        "num_rows":  num_rows,
        "start_row": start_row,
    }
    if only:
        params["only"] = ",".join(only)
    if include:
        params["include"] = include

    response = requests.get(_RMA_URL, params=params, timeout=30)
    response.raise_for_status()

    payload = response.json()
    if not payload.get("success", True):
        raise ValueError(f"Allen API error: {payload.get('msg', 'unknown error')}")

    return payload.get("msg", [])


def _rma_query_all(criteria: str, only: list = None, include: str = None) -> list:
    """
    Paginated RMA query — fetches ALL results across multiple pages.
    Automatically handles large result sets that exceed _PAGE_SIZE.
    """
    results   = []
    start_row = 0

    while True:
        page = _rma_query(criteria, only=only, include=include,
                          num_rows=_PAGE_SIZE, start_row=start_row)
        results.extend(page)

        if len(page) < _PAGE_SIZE:
            break   # Last page reached

        start_row += _PAGE_SIZE
        time.sleep(0.1)   # Polite rate limiting

    return results


class GeneExpressionLoader:
    """
    Retrieves Allen ISH gene expression data per brain structure.

    Caches API results as JSON files in cache_dir to avoid redundant
    network calls across runs.

    Parameters
    ----------
    cache_dir : str
        Directory for caching downloaded expression data.
    """

    def __init__(self, cache_dir: str = "allen_connectivity_cache"):
        self.cache_dir = Path(cache_dir) / "ish_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _cache_path(self, key: str) -> Path:
        """Return path to a cache file for a given query key."""
        safe_key = key.replace("/", "_").replace(" ", "_")
        return self.cache_dir / f"{safe_key}.json"

    def _read_cache(self, key: str):
        """Return cached result or None if not cached."""
        path = self._cache_path(key)
        if path.exists():
            with open(path) as f:
                return json.load(f)
        return None

    def _write_cache(self, key: str, data) -> None:
        """Persist query result to cache."""
        with open(self._cache_path(key), "w") as f:
            json.dump(data, f)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_section_datasets(self, gene_acronym: str) -> list:
        """
        Find all ISH section data set IDs for one gene in adult mouse brain.

        The Allen ISH atlas contains multiple section data sets per gene
        (different planes of section, different animals). We use all of them
        and later average expression across data sets.

        Parameters
        ----------
        gene_acronym : str
            Allen gene acronym (case-sensitive), e.g. 'Oprm1', 'Pdyn'.

        Returns
        -------
        list of int
            Section data set IDs. Empty list if gene not found.
        """
        cache_key = f"datasets_{gene_acronym}"
        cached    = self._read_cache(cache_key)
        if cached is not None:
            return cached

        # Query: find non-failed section data sets for this gene in Mouse product
        criteria = (
            f"model::SectionDataSet,"
            f"rma::criteria,"
            f"[failed$eq'false'],"
            f"products[abbreviation$eq'{_ISH_PRODUCT}'],"
            f"genes[acronym$eq'{gene_acronym}']"
        )
        try:
            results = _rma_query_all(criteria, only=["id", "plane_of_section_id"])
        except Exception as exc:
            logger.error("Error fetching datasets for '%s': %s", gene_acronym, exc)
            return []

        dataset_ids = [r["id"] for r in results]
        logger.info(
            "  Gene '%s': %d ISH section data sets found.", gene_acronym, len(dataset_ids)
        )

        self._write_cache(cache_key, dataset_ids)
        return dataset_ids

    def get_expression_by_structure(
        self, dataset_ids: list, structure_ids: list
    ) -> pd.DataFrame:
        """
        Retrieve expression_energy per structure for a set of ISH data sets.

        Parameters
        ----------
        dataset_ids : list of int
            ISH section data set IDs (from get_section_datasets).
        structure_ids : list of int
            Allen structure IDs to query.

        Returns
        -------
        pd.DataFrame
            Columns: dataset_id, structure_id, expression_energy, expression_density.
            One row per (dataset, structure) pair.
        """
        if not dataset_ids or not structure_ids:
            return pd.DataFrame()

        # Build cache key from dataset IDs AND structure IDs — both must match for
        # the cache to be valid. Omitting structure_ids caused stale hits when the
        # same gene was queried for a different region set.
        structs_hash = abs(hash(tuple(sorted(structure_ids)))) % (10**8)
        cache_key = (
            f"expr_{'_'.join(str(i) for i in sorted(dataset_ids[:5]))}"
            f"_n{len(dataset_ids)}_s{structs_hash}"
        )
        cached    = self._read_cache(cache_key)
        if cached is not None:
            return pd.DataFrame(cached)

        ids_str  = ",".join(str(i) for i in dataset_ids)
        structs_str = ",".join(str(i) for i in structure_ids)

        criteria = (
            f"model::StructureUnionize,"
            f"rma::criteria,"
            f"section_data_set[id$in{ids_str}],"
            f"structure[id$in{structs_str}]"
        )
        only = [
            "section_data_set_id",
            "structure_id",
            "expression_energy",
            "expression_density",
            "sum_expressing_pixel_intensity",
        ]

        try:
            results = _rma_query_all(criteria, only=only)
        except Exception as exc:
            logger.error("Error fetching expression data: %s", exc)
            return pd.DataFrame()

        df = pd.DataFrame(results)
        if df.empty:
            return df

        # Standardize column names (API returns camelCase or hyphenated)
        df = df.rename(columns={
            "section_data_set_id":          "dataset_id",
            "sectionDataSetId":             "dataset_id",
            "expressionEnergy":             "expression_energy",
            "expressionDensity":            "expression_density",
            "structureId":                  "structure_id",
        })

        # Ensure numeric dtypes
        for col in ["expression_energy", "expression_density"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        self._write_cache(cache_key, df.to_dict(orient="records"))
        return df

    def load(self, genes: list, structure_ids: list) -> pd.DataFrame:
        """
        Main entry point: build a (region × gene) expression matrix.

        Fetches ISH expression data for each gene, averages across section
        data sets, and returns a clean DataFrame ready for analysis.

        Parameters
        ----------
        genes : list of str
            Allen gene acronyms, e.g. ['Oprm1', 'Oprd1', 'Pdyn'].
        structure_ids : list of int
            Allen structure IDs for target brain regions.

        Returns
        -------
        pd.DataFrame
            Index = Allen structure IDs, Columns = gene acronyms.
            Values = mean expression_energy across ISH data sets.
            Missing values filled with 0.0.
        """
        # Initialize empty result matrix
        expr_matrix = pd.DataFrame(
            0.0,
            index   = structure_ids,
            columns = genes,
        )

        for gene in genes:
            logger.info("Loading expression for gene: %s", gene)
            dataset_ids = self.get_section_datasets(gene)

            if not dataset_ids:
                logger.warning("No ISH data sets found for '%s' — column will be zero.", gene)
                continue

            raw = self.get_expression_by_structure(dataset_ids, structure_ids)

            if raw.empty:
                logger.warning("No expression data returned for '%s'.", gene)
                continue

            # Average expression_energy across data sets for each structure
            avg = (
                raw
                .groupby("structure_id")["expression_energy"]
                .mean()
                .reindex(structure_ids, fill_value=0.0)
            )
            expr_matrix[gene] = avg.values

        return expr_matrix
