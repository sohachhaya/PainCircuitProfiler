"""
data/connectivity.py — Allen Mouse Connectivity Atlas interface
===============================================================
Queries the Allen Connectivity Atlas for anterograde tracing data.
Builds a region × region projection strength matrix.

DATA SOURCE:
    Allen Mouse Brain Connectivity Atlas (Oh et al., Nature 2014)
    https://connectivity.brain-map.org/

KEY CONCEPT:
    Each "experiment" = one mouse with AAV-GFP injected into a specific brain region.
    After 3 weeks, axons fluoresce wherever the injected neurons project.
    projection_energy = normalized fluorescence density in a target region.

MAIN ENTRY POINT:
    loader = ConnectivityLoader(cache_dir='allen_connectivity_cache')
    matrix, metadata = loader.load(regions=PAIN_REGIONS)
"""

import logging
import numpy as np
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

logger = logging.getLogger(__name__)


class ConnectivityLoader:
    """
    Retrieves and caches Allen Connectivity Atlas projection data.

    The class handles:
        1. Structure ID lookup (acronym → Allen numeric ID)
        2. Experiment retrieval per injection region
        3. Projection energy queries per experiment × target region pair
        4. Averaging across biological replicates (multiple mice per region)
        5. Returning a clean pandas DataFrame matrix

    Parameters
    ----------
    cache_dir : str
        Directory for storing downloaded Allen data files.
        Created automatically if it does not exist.
    resolution : int
        Atlas resolution in microns. 100 is sufficient for region-level analysis.
    """

    def __init__(self, cache_dir: str = "allen_connectivity_cache", resolution: int = 100):
        self.cache_dir  = cache_dir
        self.resolution = resolution
        self._mcache    = None   # Lazy initialization — only connect when needed

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @property
    def mcache(self) -> MouseConnectivityCache:
        """Lazy-load the MouseConnectivityCache (connects to Allen on first use)."""
        if self._mcache is None:
            self._mcache = MouseConnectivityCache(
                manifest_file=f"{self.cache_dir}/mouse_connectivity_manifest.json",
                resolution=self.resolution,
            )
            logger.info("Connected to Allen Mouse Connectivity Atlas.")
        return self._mcache

    @property
    def structure_tree(self):
        """Cached structure tree (all brain region IDs and hierarchy)."""
        if not hasattr(self, "_structure_tree"):
            self._structure_tree = self.mcache.get_structure_tree()
        return self._structure_tree

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_structure_ids(self, acronyms: list) -> dict:
        """
        Look up Allen numeric IDs for a list of region acronyms.

        Queries each acronym individually so that one unknown acronym does
        not abort the entire lookup (allensdk raises KeyError for missing
        entries when querying a batch).

        Returns
        -------
        dict
            {acronym: structure_id} — missing acronyms are omitted with a warning.
        """
        found = {}
        for acr in acronyms:
            try:
                results = self.structure_tree.get_structures_by_acronym([acr])
                if results:
                    found[acr] = results[0]["id"]
                else:
                    logger.warning("Region '%s' not found in Allen structure tree — skipping.", acr)
            except (KeyError, IndexError):
                logger.warning("Region '%s' not found in Allen structure tree — skipping.", acr)
        return found

    def get_experiments(self, structure_id: int) -> list:
        """
        Retrieve all injection experiments for one source region.

        Returns
        -------
        list of dict
            Each dict: experiment metadata (id, coordinates, product, etc.).
        """
        return self.mcache.get_experiments(injection_structure_ids=[structure_id])

    def build_matrix(
        self,
        acronym_to_id: dict,
        label_map: dict,
        hemisphere_ids: list = None,
    ) -> pd.DataFrame:
        """
        Build a connectivity matrix: rows = source, columns = target.

        For each source region, retrieves all experiments and queries projection
        energy into each target region. Values are averaged across experiments
        (biological replicates from multiple mice).

        Parameters
        ----------
        acronym_to_id : dict
            {acronym: allen_id} for all regions to include.
        label_map : dict
            {acronym: display_label} for row/column naming.
        hemisphere_ids : list, optional
            Filter by hemisphere: [1]=left, [2]=right, [3]=both.
            Default: both hemispheres (no filter).

        Returns
        -------
        pd.DataFrame
            Symmetric-shaped matrix (n_regions × n_regions).
            Rows = source (injection site), Columns = target (projection site).
            Values = mean projection_energy.
        """
        all_ids    = list(acronym_to_id.values())
        all_labels = [label_map.get(acr, acr) for acr in acronym_to_id]
        records    = []

        for acronym, source_id in acronym_to_id.items():
            source_label = label_map.get(acronym, acronym)
            experiments  = self.get_experiments(source_id)

            if not experiments:
                logger.warning("No experiments found for '%s' (%s).", acronym, source_label)
                continue

            exp_ids = [e["id"] for e in experiments]
            logger.info(
                "  %s (%s): %d experiments", source_label, acronym, len(exp_ids)
            )

            try:
                unionizes = self.mcache.get_structure_unionizes(
                    experiment_ids      = exp_ids,
                    is_injection        = False,      # Exclude injection site fluorescence
                    structure_ids       = all_ids,    # Target = all pain regions
                    include_descendants = False,
                )
            except Exception as exc:
                logger.error("Error querying unionizes for '%s': %s", acronym, exc)
                continue

            if unionizes.empty:
                continue

            # Apply hemisphere filter if requested
            if hemisphere_ids and "hemisphere_id" in unionizes.columns:
                unionizes = unionizes[unionizes["hemisphere_id"].isin(hemisphere_ids)]

            # Tag with source for later grouping
            unionizes["source_acronym"] = acronym
            unionizes["source_label"]   = source_label
            records.append(unionizes)

        if not records:
            logger.error("No connectivity data retrieved — returning empty matrix.")
            return pd.DataFrame(0.0, index=all_labels, columns=all_labels)

        master = pd.concat(records, ignore_index=True)

        # Translate target structure IDs → labels
        id_to_label                = {v: label_map.get(k, k) for k, v in acronym_to_id.items()}
        master["target_label"]     = master["structure_id"].map(id_to_label)
        master                     = master.dropna(subset=["target_label"])

        # Average projection_energy across experiments for the same source→target pair
        avg = (
            master
            .groupby(["source_label", "target_label"])["projection_energy"]
            .mean()
            .reset_index()
        )

        # Pivot to matrix format
        matrix = avg.pivot_table(
            index      = "source_label",
            columns    = "target_label",
            values     = "projection_energy",
            fill_value = 0.0,
        )
        matrix.columns.name = None
        matrix.index.name   = None

        # Enforce consistent row/column ordering
        matrix = matrix.reindex(index=all_labels, columns=all_labels, fill_value=0.0)

        return matrix

    def load(self, regions: dict, hemisphere_ids: list = None) -> tuple:
        """
        Main entry point: load connectivity matrix for a set of brain regions.

        Parameters
        ----------
        regions : dict
            Region definitions from config.py (PAIN_REGIONS), or a custom dict
            with structure {acronym: {label: str, ...}}.
        hemisphere_ids : list, optional
            Restrict to specific hemisphere(s).

        Returns
        -------
        matrix : pd.DataFrame
            Connectivity matrix (display labels as row/column names).
        metadata : dict
            {acronym: {id, label, n_experiments}} for each found region.
        """
        acronyms      = list(regions.keys())
        label_map     = {acr: info["label"] for acr, info in regions.items()}
        acronym_to_id = self.get_structure_ids(acronyms)

        # ---- Report any dropped regions clearly ----
        dropped = [acr for acr in acronyms if acr not in acronym_to_id]
        if dropped:
            print("\n⚠  REGIONS NOT FOUND IN ALLEN CCF — excluded from analysis:")
            print("   (Check spelling against the Allen Brain Atlas structure tree)")
            print("   https://atlas.brain-map.org/\n")
            for acr in dropped:
                print(f"   ✗  '{acr}'  ({label_map.get(acr, '?')})")
            print()

        logger.info(
            "Loading connectivity matrix for %d / %d regions (%d dropped)...",
            len(acronym_to_id), len(acronyms), len(dropped),
        )

        matrix = self.build_matrix(acronym_to_id, label_map, hemisphere_ids)

        # Build metadata dict for reporting
        # Note: experiments were already fetched inside build_matrix; we query
        # them once more here only for the count. The SDK caches locally so
        # this does not make extra network requests.
        metadata = {}
        for acr, sid in acronym_to_id.items():
            exps = self.get_experiments(sid)
            metadata[acr] = {
                "id":            sid,
                "label":         label_map[acr],
                "n_experiments": len(exps),
            }

        # Also record dropped regions in metadata so callers can inspect them
        for acr in dropped:
            metadata[acr] = {
                "id":            None,
                "label":         label_map.get(acr, acr),
                "n_experiments": 0,
                "dropped":       True,
                "reason":        "Acronym not found in Allen CCF structure tree",
            }

        return matrix, metadata
