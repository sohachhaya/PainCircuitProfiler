"""
PainCircuitProfiler
===================
A tool for integrating Allen Brain Atlas connectivity and gene expression
data to analyse pain processing circuits in the mouse brain.

Typical usage
-------------
from pain_circuit_profiler import profile_pain_circuit
results = profile_pain_circuit()

Or step-by-step:
    from pain_circuit_profiler.data.connectivity   import ConnectivityLoader
    from pain_circuit_profiler.data.gene_expression import GeneExpressionLoader
    from pain_circuit_profiler.analysis.hub_metrics import run_analysis
    from pain_circuit_profiler.visualization.plots  import save_all_figures
"""

import logging
from pathlib import Path

from .config import PAIN_REGIONS, GENES, ANALYSIS_DEFAULTS
from .data.connectivity    import ConnectivityLoader
from .data.gene_expression import GeneExpressionLoader
from .analysis.hub_metrics import run_analysis, print_top_hubs
from .visualization.plots  import (
    plot_expression_heatmap,
    plot_connectivity_heatmap,
    plot_network,
    plot_hub_scores,
    plot_integrated_dashboard,
    save_all_figures,
)

__version__ = "0.1.0"
__author__  = "PainCircuitProfiler"

logging.basicConfig(
    level  = logging.INFO,
    format = "%(asctime)s  %(levelname)-8s  %(name)s — %(message)s",
    datefmt= "%H:%M:%S",
)


def profile_pain_circuit(
    regions: dict = None,
    genes:   list = None,
    cache_dir:   str = None,
    output_dir:  str = None,
    expression_percentile: float = None,
    projection_percentile: float = None,
    save_figures: bool = True,
) -> dict:
    """
    Run the full PainCircuitProfiler pipeline in one call.

    Steps
    -----
    1. Load connectivity matrix from Allen Mouse Connectivity Atlas.
    2. Look up Allen structure IDs for the region set.
    3. Load ISH gene expression data from Allen Mouse Brain Atlas.
    4. Run hub metric analysis (expression × connectivity integration).
    5. Optionally save all figures to output_dir.

    Parameters
    ----------
    regions : dict, optional
        Region definitions. Defaults to PAIN_REGIONS from config.py.
        Format: {acronym: {label: str, category: str, description: str}}.
    genes : list of str, optional
        Gene acronyms to query. Defaults to list(GENES.keys()).
    cache_dir : str, optional
        Local cache directory. Defaults to ANALYSIS_DEFAULTS['cache_dir'].
    output_dir : str, optional
        Figure output directory. Defaults to ANALYSIS_DEFAULTS['output_dir'].
    expression_percentile : float, optional
        Percentile above which a region is an "expression hub". Default: 75.
    projection_percentile : float, optional
        Percentile above which a region is a "hub candidate". Default: 50.
    save_figures : bool
        Whether to save figures to output_dir. Default: True.

    Returns
    -------
    dict with keys:
        connectivity  : pd.DataFrame — region × region projection matrix
        expression    : pd.DataFrame — region × gene expression matrix
        conn_metadata : dict         — experiment counts per region
        hub_results   : dict         — output of run_analysis()
        figure_paths  : dict         — {figure_name: file_path} (if saved)
    """
    regions    = regions    or PAIN_REGIONS
    genes      = genes      or list(GENES.keys())
    cache_dir  = cache_dir  or ANALYSIS_DEFAULTS["cache_dir"]
    output_dir = output_dir or ANALYSIS_DEFAULTS["output_dir"]
    expr_pct   = expression_percentile or ANALYSIS_DEFAULTS["expression_percentile"]
    proj_pct   = projection_percentile or ANALYSIS_DEFAULTS["projection_percentile"]

    logger = logging.getLogger(__name__)
    logger.info("=" * 55)
    logger.info("PainCircuitProfiler v%s", __version__)
    logger.info("Regions: %d  |  Genes: %d", len(regions), len(genes))
    logger.info("=" * 55)

    # ---- Step 1: Connectivity ----
    logger.info("Step 1/4  Loading connectivity data...")
    conn_loader = ConnectivityLoader(cache_dir=cache_dir)
    connectivity, conn_metadata = conn_loader.load(regions)

    # Build id_to_label for expression index translation
    acronym_to_id = conn_loader.get_structure_ids(list(regions.keys()))
    id_to_label   = {v: regions[k]["label"] for k, v in acronym_to_id.items()}

    # ---- Step 2: Gene expression ----
    logger.info("Step 2/4  Loading gene expression data...")
    expr_loader = GeneExpressionLoader(cache_dir=cache_dir)
    expression  = expr_loader.load(genes, structure_ids=list(acronym_to_id.values()))

    # ---- Step 3: Hub analysis ----
    logger.info("Step 3/4  Running hub analysis...")
    hub_results = run_analysis(
        connectivity,
        expression,
        id_to_label          = id_to_label,
        expression_percentile = expr_pct,
        projection_percentile = proj_pct,
    )
    print_top_hubs(hub_results)

    # Rename expression index from Allen IDs → display labels (for plotting)
    expression.index = [id_to_label.get(i, str(i)) for i in expression.index]

    # ---- Step 4: Figures ----
    figure_paths = {}
    if save_figures:
        logger.info("Step 4/4  Saving figures to %s/...", output_dir)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        figure_paths = save_all_figures(connectivity, expression, hub_results, output_dir)
    else:
        logger.info("Step 4/4  Skipping figure saving (save_figures=False).")

    return {
        "connectivity":   connectivity,
        "expression":     expression,
        "conn_metadata":  conn_metadata,
        "hub_results":    hub_results,
        "figure_paths":   figure_paths,
    }
