"""
analysis/hub_metrics.py — Hub metric computation for pain circuit analysis
==========================================================================
Integrates connectivity and gene expression data to identify circuit "hubs":
regions that both project strongly to high-expression targets AND are
themselves highly interconnected within the pain network.

SCIENTIFIC RATIONALE:
    A region is considered a pain circuit hub if it satisfies two criteria:
    (1) Projection relevance: it sends strong axonal projections to regions
        where a target molecule (e.g., Oprm1) is highly expressed.
        Interpretation: its downstream targets are opioid-rich → modulating
        these projections could shift opioid signaling.
    (2) Network centrality: among regions meeting criterion (1), it is highly
        interconnected with other pain regions (high in-degree + out-degree).
        Interpretation: it can broadly influence/be influenced by the network.

    Hub Score = projection_score × interconnectivity_score
    (geometric mean, so a region must score on BOTH dimensions)

    This integrates both molecular accessibility (expression) and
    circuit architecture (connectivity) in a single ranked output.

MAIN ENTRY POINT:
    results = run_analysis(connectivity_matrix, expression_df, id_to_label={allen_id: label})
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# Core metric functions
# =============================================================================

def find_expression_hubs(
    expression_series: pd.Series,
    percentile: float = 75,
    min_energy: float = 0.0,
) -> list:
    """
    Identify regions with high expression of a target gene.

    "High expression" is defined as: expression_energy > percentile threshold
    AND expression_energy > min_energy (absolute floor).

    Parameters
    ----------
    expression_series : pd.Series
        Index = region labels, values = expression_energy for one gene.
    percentile : float
        Percentile threshold (0–100). Regions above this are "high expression".
    min_energy : float
        Absolute minimum expression_energy required.

    Returns
    -------
    list of str
        Region labels classified as expression hubs.
    """
    if expression_series.empty or expression_series.max() == 0:
        logger.warning("Expression series is empty or all zeros.")
        return []

    threshold = np.percentile(expression_series.values, percentile)
    threshold = max(threshold, min_energy)

    hubs = expression_series[expression_series >= threshold].index.tolist()
    logger.debug(
        "Expression hubs (top %g%%): %d regions above %.4f", 100 - percentile, len(hubs), threshold
    )
    return hubs


def compute_projection_to_hubs(
    connectivity: pd.DataFrame,
    hubs: list,
) -> pd.Series:
    """
    For each region, compute its total projection strength to expression hubs.

    Sums projection_energy across all hub target regions.
    Regions not in the connectivity matrix are assigned 0.

    Parameters
    ----------
    connectivity : pd.DataFrame
        Square connectivity matrix (rows=source, columns=target).
    hubs : list of str
        Target regions identified as expression hubs.

    Returns
    -------
    pd.Series
        Index = source region labels, values = summed projection to hubs.
    """
    # Only include hub columns that exist in the connectivity matrix
    valid_hubs = [h for h in hubs if h in connectivity.columns]

    if not valid_hubs:
        logger.warning("None of the expression hubs are in the connectivity matrix.")
        return pd.Series(0.0, index=connectivity.index)

    # Sum projection energy across all hub target columns
    projection_scores = connectivity[valid_hubs].sum(axis=1)
    return projection_scores


def compute_interconnectivity(
    connectivity: pd.DataFrame,
    regions: list,
) -> pd.Series:
    """
    For each region in `regions`, compute its interconnectivity within the group.

    Interconnectivity score for region R = mean projection energy received FROM
    other regions in the group + mean projection energy sent TO other regions.
    (i.e., normalized bidirectional connectivity within the candidate set)

    Parameters
    ----------
    connectivity : pd.DataFrame
        Full connectivity matrix.
    regions : list of str
        Subset of regions to evaluate interconnectivity within.

    Returns
    -------
    pd.Series
        Index = region labels (all regions in connectivity matrix),
        values = interconnectivity score (0 for regions not in `regions`).
    """
    valid_regions = [r for r in regions if r in connectivity.index and r in connectivity.columns]
    n = len(valid_regions)

    scores = pd.Series(0.0, index=connectivity.index)

    if n < 2:
        logger.warning("Fewer than 2 valid hub candidates — interconnectivity not meaningful.")
        return scores

    for region in valid_regions:
        # Outgoing: how strongly does this region project to OTHER candidates?
        outgoing = connectivity.loc[region, [r for r in valid_regions if r != region]].mean()
        # Incoming: how strongly do OTHER candidates project to this region?
        incoming = connectivity.loc[[r for r in valid_regions if r != region], region].mean()
        scores[region] = (outgoing + incoming) / 2.0

    return scores


def compute_hub_scores(
    connectivity: pd.DataFrame,
    expression_series: pd.Series,
    expression_percentile: float = 75,
    projection_percentile: float = 50,
) -> pd.DataFrame:
    """
    Compute hub scores integrating connectivity and gene expression.

    Pipeline:
        1. expression_hubs   ← regions above expression_percentile
        2. projection_scores ← sum of projections to expression_hubs (per region)
        3. hub_candidates    ← regions above projection_percentile in step 2
        4. interconnect      ← bidirectional connectivity within hub_candidates
        5. hub_score         ← geometric mean(projection_score, interconnect)

    Parameters
    ----------
    connectivity : pd.DataFrame
        Connectivity matrix (display labels).
    expression_series : pd.Series
        Expression energy for one gene. Index = region display labels.
    expression_percentile : float
        Top X% expression defines expression hubs (default: top 25%).
    projection_percentile : float
        Top X% projection defines hub candidates (default: top 50%).

    Returns
    -------
    pd.DataFrame
        Columns: expression_energy, projection_to_hubs, interconnectivity, hub_score.
        Sorted by hub_score descending.
    """
    # Step 1: Find expression hubs
    expression_hubs = find_expression_hubs(expression_series, percentile=expression_percentile)

    # Step 2: Projection scores to expression hubs
    projection_scores = compute_projection_to_hubs(connectivity, expression_hubs)

    # Step 3: Identify hub candidates (regions with strong projection to hubs)
    if projection_scores.max() > 0:
        proj_threshold    = np.percentile(
            projection_scores[projection_scores > 0].values, projection_percentile
        )
        hub_candidates    = projection_scores[projection_scores >= proj_threshold].index.tolist()
    else:
        hub_candidates = []
        logger.warning("All projection scores are zero.")

    # Step 4: Interconnectivity within hub candidates
    interconnect = compute_interconnectivity(connectivity, hub_candidates)

    # Step 5: Hub score = geometric mean of normalized projection × normalized interconnect
    # Normalize to [0, 1] range for each component before combining
    proj_norm        = _safe_normalize(projection_scores)
    interconn_norm   = _safe_normalize(interconnect)

    # Geometric mean: both dimensions must be non-zero to produce a high score
    hub_score        = np.sqrt(proj_norm * interconn_norm)

    # Assemble output DataFrame — reindex to common set of regions
    all_regions = connectivity.index.tolist()
    results = pd.DataFrame({
        "expression_energy":   expression_series.reindex(all_regions, fill_value=0.0),
        "is_expression_hub":   [r in expression_hubs for r in all_regions],
        "projection_to_hubs":  projection_scores.reindex(all_regions, fill_value=0.0),
        "is_hub_candidate":    [r in hub_candidates for r in all_regions],
        "interconnectivity":   interconnect.reindex(all_regions, fill_value=0.0),
        "hub_score":           hub_score.reindex(all_regions, fill_value=0.0),
    }, index=all_regions)

    return results.sort_values("hub_score", ascending=False)


def _safe_normalize(series: pd.Series) -> pd.Series:
    """Normalize a Series to [0, 1]; returns zeros if max == 0."""
    max_val = series.max()
    if max_val == 0:
        return pd.Series(0.0, index=series.index)
    return series / max_val


# =============================================================================
# Full pipeline
# =============================================================================

def run_analysis(
    connectivity: pd.DataFrame,
    expression_df: pd.DataFrame,
    id_to_label: dict = None,
    expression_percentile: float = 75,
    projection_percentile: float = 50,
) -> dict:
    """
    Run the full hub analysis for all genes in expression_df.

    For each gene, computes hub scores and returns a summary DataFrame
    plus per-gene detailed results.

    Parameters
    ----------
    connectivity : pd.DataFrame
        Connectivity matrix with region display labels.
    expression_df : pd.DataFrame
        Rows = region IDs (or labels if id_to_label is None).
        Columns = gene acronyms.
    id_to_label : dict, optional
        {allen_id: display_label}. If provided, translates expression_df index.
    expression_percentile : float
        Top X% for expression hub definition.
    projection_percentile : float
        Top X% for hub candidate definition.

    Returns
    -------
    dict with keys:
        "per_gene"  : {gene: pd.DataFrame of hub scores}
        "summary"   : pd.DataFrame — one row per region, columns = hub_score per gene
        "top_hubs"  : pd.DataFrame — regions ranked by mean hub score across all genes
    """
    # Translate expression_df index from Allen IDs to display labels if needed
    if id_to_label is not None:
        expression_df = expression_df.rename(index=id_to_label)

    # Restrict to regions present in both connectivity and expression
    shared_regions = [r for r in connectivity.index if r in expression_df.index]
    if not shared_regions:
        logger.error(
            "No shared regions between connectivity (%d) and expression (%d) DataFrames.",
            len(connectivity), len(expression_df),
        )
        return {"per_gene": {}, "summary": pd.DataFrame(), "top_hubs": pd.DataFrame()}

    conn_sub = connectivity.loc[shared_regions, shared_regions]
    expr_sub = expression_df.loc[shared_regions]

    per_gene = {}
    genes    = expression_df.columns.tolist()

    for gene in genes:
        if gene not in expr_sub.columns:
            continue
        logger.info("Running hub analysis for: %s", gene)
        scores = compute_hub_scores(
            conn_sub,
            expr_sub[gene],
            expression_percentile=expression_percentile,
            projection_percentile=projection_percentile,
        )
        per_gene[gene] = scores

    if not per_gene:
        return {"per_gene": {}, "summary": pd.DataFrame(), "top_hubs": pd.DataFrame()}

    # Summary: hub_score for each gene as separate columns
    summary = pd.DataFrame(
        {gene: per_gene[gene]["hub_score"] for gene in per_gene},
        index=shared_regions,
    )
    summary.columns = [f"hub_{g}" for g in summary.columns]

    # Top hubs: ranked by mean hub score across all genes
    summary["mean_hub_score"] = summary.mean(axis=1)
    top_hubs = summary.sort_values("mean_hub_score", ascending=False)

    return {
        "per_gene":  per_gene,
        "summary":   summary,
        "top_hubs":  top_hubs,
    }


def print_top_hubs(results: dict, top_n: int = 10) -> None:
    """
    Print a formatted summary of the top hub regions.

    Parameters
    ----------
    results : dict
        Output from run_analysis().
    top_n : int
        Number of top regions to display.
    """
    top = results["top_hubs"].head(top_n)

    print(f"\n{'='*65}")
    print(f"TOP {top_n} PAIN CIRCUIT HUBS")
    print(f"{'='*65}")
    print(f"Ranked by mean hub score across all analyzed genes.\n")
    print(f"  {'Region':<25}  {'Mean Hub Score':>15}  {'Genes Driving Score'}")
    print("  " + "-" * 60)

    gene_cols = [c for c in top.columns if c.startswith("hub_") and c != "mean_hub_score"]

    for region, row in top.iterrows():
        # Which genes contribute the highest hub score for this region?
        gene_scores = {col.replace("hub_", ""): row[col] for col in gene_cols}
        top_genes   = sorted(gene_scores, key=gene_scores.get, reverse=True)[:3]
        gene_str    = ", ".join(f"{g}({gene_scores[g]:.2f})" for g in top_genes)
        print(f"  {region:<25}  {row['mean_hub_score']:>15.4f}  {gene_str}")

    print()
