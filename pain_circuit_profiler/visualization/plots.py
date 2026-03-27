"""
visualization/plots.py — All figure generation for PainCircuitProfiler
=======================================================================
Provides five plotting functions covering:
    1. Gene expression heatmap across pain regions
    2. Connectivity heatmap (log-scaled projection energies)
    3. Network diagram (connectivity + expression overlaid)
    4. Hub score rankings per gene (horizontal bar chart)
    5. Integrated dashboard (4-panel figure combining all views)

DESIGN PRINCIPLES:
    - Every function returns a matplotlib Figure object so callers can
      further customize or embed figures in notebooks.
    - All functions accept an optional `output_path` to save the figure.
    - A consistent visual style (color palette, font sizes) is set once
      at module load and applied throughout.
    - Category colors from config.py are used consistently across all plots
      so the same region always appears in the same color.

DEPENDENCIES:
    matplotlib, seaborn, networkx, numpy, pandas
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns

from ..config import CATEGORY_COLORS, PAIN_REGIONS, GENES, ANALYSIS_DEFAULTS

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Module-level style — applied to all figures in this module
# ---------------------------------------------------------------------------
sns.set_style("ticks")
sns.set_context("notebook", font_scale=1.15)

_FIGSIZE_SQUARE    = (10, 9)
_FIGSIZE_WIDE      = (14, 7)
_FIGSIZE_DASHBOARD = (18, 16)
_DPI               = ANALYSIS_DEFAULTS["figure_dpi"]
_CONN_CMAP         = ANALYSIS_DEFAULTS["connectivity_cmap"]
_EXPR_CMAP         = ANALYSIS_DEFAULTS["expression_cmap"]


def _node_color_for_region(region_label: str) -> str:
    """Return category color for a region label, grey if not found."""
    for acronym, info in PAIN_REGIONS.items():
        if info["label"] == region_label or acronym == region_label:
            return CATEGORY_COLORS.get(info["category"], "#95a5a6")
    return "#95a5a6"


def _save(fig: plt.Figure, output_path) -> None:
    """Save figure if output_path is given."""
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=_DPI, bbox_inches="tight")
        logger.info("Figure saved: %s", output_path)


def _add_colorbar_legend(ax, cmap, vmin, vmax, label, shrink=0.7):
    """Add a standalone colorbar to an axes."""
    norm   = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm     = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar   = plt.colorbar(sm, ax=ax, shrink=shrink, pad=0.02)
    cbar.set_label(label, fontsize=10)
    return cbar


# =============================================================================
# 1. Gene Expression Heatmap
# =============================================================================

def plot_expression_heatmap(
    expression_df: pd.DataFrame,
    title: str = "Opioid Gene Expression Across Pain Regions",
    output_path: str = None,
) -> plt.Figure:
    """
    Heatmap of gene expression energy: regions (rows) × genes (columns).

    Color encodes expression_energy. Regions are ordered by anatomical
    category (matching the category groupings in config.py).
    Tick labels on the y-axis are color-coded by brain region category.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Index = region display labels, columns = gene acronyms.
        Values = mean expression_energy (from GeneExpressionLoader).
    title : str
        Figure title.
    output_path : str, optional
        If given, save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if expression_df.empty:
        logger.warning("plot_expression_heatmap: empty DataFrame — skipping.")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        return fig

    fig, ax = plt.subplots(figsize=(max(8, len(expression_df.columns) * 1.5),
                                    max(6, len(expression_df) * 0.55)))

    sns.heatmap(
        expression_df,
        ax         = ax,
        cmap       = _EXPR_CMAP,
        annot      = True,
        fmt        = ".3f",
        linewidths = 0.4,
        linecolor  = "whitesmoke",
        cbar_kws   = {"label": "Expression Energy (ISH)", "shrink": 0.75},
        annot_kws  = {"size": 8},
    )

    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Brain Region", fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right", fontsize=10)

    # Color-code y-axis tick labels by anatomical category
    for tick_label in ax.get_yticklabels():
        region = tick_label.get_text()
        tick_label.set_color(_node_color_for_region(region))
        tick_label.set_fontsize(10)

    # Category legend
    legend_patches = [
        mpatches.Patch(color=color, label=cat)
        for cat, color in CATEGORY_COLORS.items()
    ]
    ax.legend(
        handles    = legend_patches,
        loc        = "upper left",
        bbox_to_anchor = (1.18, 1.0),
        fontsize   = 9,
        title      = "Category",
        title_fontsize = 10,
        framealpha = 0.9,
    )

    # Gene annotation below x-axis: receptor vs peptide
    gene_roles = {g: info.get("pathway", "") for g, info in GENES.items()}
    for tick_label in ax.get_xticklabels():
        gene = tick_label.get_text()
        if "receptor" in gene_roles.get(gene, "").lower():
            tick_label.set_fontstyle("italic")

    plt.tight_layout()
    _save(fig, output_path)
    return fig


# =============================================================================
# 2. Connectivity Heatmap
# =============================================================================

def plot_connectivity_heatmap(
    connectivity: pd.DataFrame,
    title: str = "Pain Circuit Connectivity\n(Allen Mouse Connectivity Atlas)",
    log_transform: bool = True,
    output_path: str = None,
) -> plt.Figure:
    """
    Heatmap of axonal projection strength: source regions (rows) × target regions (cols).

    Values are log₁₀-transformed by default because projection energies
    span several orders of magnitude — linear scale renders weak connections
    invisible.

    Parameters
    ----------
    connectivity : pd.DataFrame
        Square connectivity matrix (display labels).
    title : str
        Figure title.
    log_transform : bool
        Apply log₁₀(x + 1e-6) transform before plotting.
    output_path : str, optional
        Save path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if connectivity.empty:
        logger.warning("plot_connectivity_heatmap: empty matrix.")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        return fig

    n   = len(connectivity)
    fig, ax = plt.subplots(figsize=(max(9, n * 0.7), max(8, n * 0.65)))

    data = np.log10(connectivity.values + 1e-6) if log_transform else connectivity.values
    data_df = pd.DataFrame(data, index=connectivity.index, columns=connectivity.columns)

    sns.heatmap(
        data_df,
        ax         = ax,
        cmap       = _CONN_CMAP,
        annot      = True,
        fmt        = ".2f",
        linewidths = 0.4,
        linecolor  = "lightgray",
        cbar_kws   = {
            "label":  "Log₁₀(Projection Energy + 1e⁻⁶)" if log_transform
                      else "Projection Energy",
            "shrink": 0.75,
        },
        annot_kws  = {"size": 7.5},
    )

    ax.set_title(title, fontsize=13, fontweight="bold", pad=18)
    ax.set_xlabel("Target Region  (receives projection)", fontsize=11, labelpad=8)
    ax.set_ylabel("Source Region  (injection site)", fontsize=11, labelpad=8)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)

    # Color-code both axes by category
    for tick_label in ax.get_xticklabels():
        tick_label.set_color(_node_color_for_region(tick_label.get_text()))
    for tick_label in ax.get_yticklabels():
        tick_label.set_color(_node_color_for_region(tick_label.get_text()))
        tick_label.set_fontsize(9)

    # Footnote
    fig.text(
        0.5, 0.01,
        "Rows = injection site  |  Columns = projection target  |  "
        "Color = connection strength  |  Log scale: each unit = 10× difference",
        ha="center", fontsize=8, style="italic", color="gray",
    )

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    _save(fig, output_path)
    return fig


# =============================================================================
# 3. Network Diagram (connectivity + expression integrated)
# =============================================================================

def plot_network(
    connectivity: pd.DataFrame,
    expression_series: pd.Series = None,
    gene_name: str = "",
    edge_threshold: float = None,
    layout: str = "spring",
    output_path: str = None,
) -> plt.Figure:
    """
    Directed network diagram with expression data overlaid on nodes.

    Visual encoding:
        Node color  : expression_energy for target gene (low=pale → high=deep red)
                      If no expression data given, nodes are colored by category.
        Node size   : out-degree (number of significant projections sent)
        Edge width  : projection_energy (stronger = thicker)
        Edge color  : projection_energy (lighter = weaker, darker = stronger)
        Arrow heads : direction of projection

    Parameters
    ----------
    connectivity : pd.DataFrame
        Connectivity matrix (display labels).
    expression_series : pd.Series, optional
        Index = region labels, values = expression_energy for one gene.
        If None, nodes are colored by anatomical category.
    gene_name : str
        Gene name shown in title and colorbar.
    edge_threshold : float, optional
        Minimum projection_energy to draw an edge.
        Defaults to ANALYSIS_DEFAULTS['network_edge_threshold'].
    layout : str
        'spring', 'circular', or 'kamada_kawai'.
    output_path : str, optional
        Save path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    threshold = edge_threshold or ANALYSIS_DEFAULTS["network_edge_threshold"]

    # Build directed graph
    G = nx.DiGraph()
    nodes = connectivity.index.tolist()
    G.add_nodes_from(nodes)

    edge_weights = {}
    for src in nodes:
        for tgt in nodes:
            if src == tgt:
                continue
            w = float(connectivity.loc[src, tgt])
            if w > threshold:
                G.add_edge(src, tgt, weight=w)
                edge_weights[(src, tgt)] = w

    # Node positions
    if layout == "spring":
        pos = nx.spring_layout(G, weight="weight", seed=42, k=2.5)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    else:
        pos = nx.kamada_kawai_layout(G)

    # ---- Node visual properties ----
    if expression_series is not None and not expression_series.empty:
        # Color by expression energy
        expr_vals = expression_series.reindex(nodes, fill_value=0.0)
        vmin, vmax = expr_vals.min(), expr_vals.max()
        if vmax == vmin:
            vmax = vmin + 1e-6
        expr_cmap  = cm.get_cmap("YlOrRd")
        node_colors = [expr_cmap((v - vmin) / (vmax - vmin)) for v in expr_vals]
        colorbar_needed = True
    else:
        node_colors     = [_node_color_for_region(n) for n in nodes]
        colorbar_needed = False
        vmin, vmax      = 0, 1

    # Node size proportional to out-degree (minimum size = 900)
    out_degrees = dict(G.out_degree())
    max_deg     = max(out_degrees.values()) if out_degrees else 1
    node_sizes  = [900 + (out_degrees.get(n, 0) / max(max_deg, 1)) * 2200 for n in nodes]

    # ---- Edge visual properties ----
    edges = list(G.edges())
    raw_w = [edge_weights.get(e, 0) for e in edges]
    if raw_w:
        max_w       = max(raw_w)
        edge_widths = [1.0 + (w / max_w) * 7.0 for w in raw_w]
        edge_colors = [cm.get_cmap("Greys")(0.3 + 0.65 * (w / max_w)) for w in raw_w]
    else:
        edge_widths = []
        edge_colors = []

    # ---- Draw ----
    fig, ax = plt.subplots(figsize=(13, 10))

    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_color = node_colors,
        node_size  = node_sizes,
        alpha      = 0.92,
    )
    nx.draw_networkx_labels(
        G, pos, ax=ax,
        font_size   = 7.5,
        font_weight = "bold",
        font_color  = "white",
    )
    if edges:
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            edge_color        = edge_colors,
            width             = edge_widths,
            arrows            = True,
            arrowsize         = 18,
            arrowstyle        = "-|>",
            connectionstyle   = "arc3,rad=0.15",
            min_source_margin = 30,
            min_target_margin = 30,
        )
    else:
        ax.text(
            0.5, 0.1,
            f"No edges above threshold ({threshold})\nTry lowering edge_threshold.",
            ha="center", transform=ax.transAxes, fontsize=11, color="gray",
        )

    # ---- Colorbar for expression ----
    if colorbar_needed:
        _add_colorbar_legend(
            ax, "YlOrRd", vmin, vmax,
            label=f"{gene_name} Expression Energy (ISH)",
        )
    else:
        # Category legend
        legend_patches = [
            mpatches.Patch(color=color, label=cat)
            for cat, color in CATEGORY_COLORS.items()
        ]
        ax.legend(handles=legend_patches, loc="upper left",
                  fontsize=8, title="Category", title_fontsize=9, framealpha=0.9)

    # ---- Node size legend ----
    size_legend = [
        plt.scatter([], [], s=900,  c="gray", alpha=0.7, label="1 projection"),
        plt.scatter([], [], s=2000, c="gray", alpha=0.7, label="~4 projections"),
        plt.scatter([], [], s=3100, c="gray", alpha=0.7, label="~7 projections"),
    ]
    ax.legend(
        handles        = size_legend,
        loc            = "lower left",
        fontsize       = 8,
        title          = "Node size = out-degree",
        title_fontsize = 9,
        framealpha     = 0.9,
        scatterpoints  = 1,
    )

    gene_str = f" — {gene_name}" if gene_name else ""
    ax.set_title(
        f"Pain Circuit Network{gene_str}\n"
        f"(edges: projection energy > {threshold}  |  "
        f"width ∝ strength  |  node color = expression)",
        fontsize=12, fontweight="bold", pad=12,
    )
    ax.axis("off")
    plt.tight_layout()
    _save(fig, output_path)
    return fig


# =============================================================================
# 4. Hub Score Rankings
# =============================================================================

def plot_hub_scores(
    per_gene_results: dict,
    genes: list = None,
    top_n: int = 12,
    output_path: str = None,
) -> plt.Figure:
    """
    Horizontal bar chart of hub scores, one panel per gene.

    Bars show hub_score per region (top N regions per gene).
    Bars are colored by anatomical category.
    Regions classified as expression hubs are starred (*).

    Parameters
    ----------
    per_gene_results : dict
        {gene: pd.DataFrame} — output of run_analysis()['per_gene'].
    genes : list, optional
        Subset of genes to plot. Defaults to all genes.
    top_n : int
        Number of top-scoring regions to show per gene.
    output_path : str, optional
        Save path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if genes is None:
        genes = list(per_gene_results.keys())
    genes = [g for g in genes if g in per_gene_results]

    if not genes:
        logger.warning("plot_hub_scores: no gene results to plot.")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No hub results available", ha="center", va="center")
        return fig

    n_genes = len(genes)
    fig, axes = plt.subplots(
        1, n_genes,
        figsize=(max(7, n_genes * 5.5), max(6, top_n * 0.45)),
        squeeze=False,
    )
    axes = axes[0]

    for ax, gene in zip(axes, genes):
        df = per_gene_results[gene].head(top_n).copy()
        if df.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(gene, fontweight="bold")
            continue

        regions   = df.index.tolist()
        scores    = df["hub_score"].values
        is_hub    = df["is_expression_hub"].values
        bar_colors = [_node_color_for_region(r) for r in regions]

        bars = ax.barh(range(len(regions)), scores, color=bar_colors,
                       edgecolor="white", linewidth=0.5, alpha=0.88)

        # Mark expression hubs with a star
        for i, (region, hub) in enumerate(zip(regions, is_hub)):
            label = f"★ {region}" if hub else f"   {region}"
            ax.text(
                -0.002, i, label,
                va="center", ha="right", fontsize=8.5,
                color=_node_color_for_region(region),
                fontweight="bold" if hub else "normal",
            )

        ax.set_yticks([])
        ax.set_xlabel("Hub Score", fontsize=10)
        ax.set_xlim(-0.02, max(scores) * 1.15 if max(scores) > 0 else 1)

        # Gene info from config
        gene_info = GENES.get(gene, {})
        gene_full = gene_info.get("name", gene)
        ax.set_title(
            f"{gene}\n{gene_full}",
            fontsize=10, fontweight="bold", pad=8,
        )

        # Shade every other bar region for readability
        for i in range(len(regions)):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color="gray", alpha=0.04, zorder=0)

        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)

    fig.suptitle(
        "Pain Circuit Hub Scores by Gene\n"
        "★ = High-expression region  |  Color = Anatomical category  |  "
        "Score = projection × interconnectivity",
        fontsize=11, fontweight="bold", y=1.02,
    )

    # Category legend on the right
    legend_patches = [
        mpatches.Patch(color=color, label=cat)
        for cat, color in CATEGORY_COLORS.items()
    ]
    fig.legend(
        handles    = legend_patches,
        loc        = "lower center",
        ncol       = len(CATEGORY_COLORS),
        fontsize   = 8,
        title      = "Anatomical Category",
        title_fontsize = 9,
        framealpha = 0.9,
        bbox_to_anchor = (0.5, -0.04),
    )

    plt.tight_layout()
    _save(fig, output_path)
    return fig


# =============================================================================
# 5. Integrated Dashboard (4-panel)
# =============================================================================

def plot_integrated_dashboard(
    connectivity: pd.DataFrame,
    expression_df: pd.DataFrame,
    hub_results: dict,
    gene: str = None,
    output_path: str = None,
) -> plt.Figure:
    """
    Four-panel integrated figure suitable for a paper supplement or README.

    Layout:
        ┌──────────────┬──────────────────┐
        │  (A)          │  (B)              │
        │  Connectivity │  Expression       │
        │  Heatmap      │  Heatmap          │
        ├──────────────┼──────────────────┤
        │  (C)          │  (D)              │
        │  Network      │  Hub Score        │
        │  Diagram      │  Rankings         │
        └──────────────┴──────────────────┘

    Parameters
    ----------
    connectivity : pd.DataFrame
        Connectivity matrix.
    expression_df : pd.DataFrame
        Expression matrix (regions × genes).
    hub_results : dict
        Output from run_analysis().
    gene : str, optional
        Gene to highlight in network (panel C). Defaults to first gene.
    output_path : str, optional
        Save path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    # Select gene for network panel
    genes = list(hub_results.get("per_gene", {}).keys())
    if gene is None:
        gene = genes[0] if genes else None

    fig = plt.figure(figsize=_FIGSIZE_DASHBOARD, constrained_layout=True)
    fig.suptitle(
        "PainCircuitProfiler — Integrated Analysis Dashboard",
        fontsize=16, fontweight="bold", y=1.01,
    )

    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35)
    ax_conn   = fig.add_subplot(gs[0, 0])
    ax_expr   = fig.add_subplot(gs[0, 1])
    ax_net    = fig.add_subplot(gs[1, 0])
    ax_hub    = fig.add_subplot(gs[1, 1])

    # ---- Panel A: Connectivity heatmap ----
    log_conn = np.log10(connectivity.values + 1e-6)
    log_df   = pd.DataFrame(log_conn, index=connectivity.index, columns=connectivity.columns)
    sns.heatmap(
        log_df, ax=ax_conn, cmap=_CONN_CMAP, annot=False,
        linewidths=0.3, linecolor="lightgray",
        cbar_kws={"label": "Log₁₀(Proj. Energy)", "shrink": 0.8},
    )
    ax_conn.set_title("(A)  Connectivity Matrix", fontweight="bold", fontsize=11)
    ax_conn.set_xticklabels(ax_conn.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    ax_conn.set_yticklabels(ax_conn.get_yticklabels(), rotation=0, fontsize=7)
    for tl in ax_conn.get_xticklabels():
        tl.set_color(_node_color_for_region(tl.get_text()))
    for tl in ax_conn.get_yticklabels():
        tl.set_color(_node_color_for_region(tl.get_text()))

    # ---- Panel B: Expression heatmap ----
    if not expression_df.empty:
        sns.heatmap(
            expression_df, ax=ax_expr, cmap=_EXPR_CMAP, annot=False,
            linewidths=0.3, linecolor="whitesmoke",
            cbar_kws={"label": "Expression Energy", "shrink": 0.8},
        )
        ax_expr.set_title("(B)  Gene Expression (ISH)", fontweight="bold", fontsize=11)
        ax_expr.set_xticklabels(ax_expr.get_xticklabels(), rotation=35, ha="right", fontsize=8)
        ax_expr.set_yticklabels(ax_expr.get_yticklabels(), rotation=0, fontsize=7)
        for tl in ax_expr.get_yticklabels():
            tl.set_color(_node_color_for_region(tl.get_text()))
    else:
        ax_expr.text(0.5, 0.5, "Expression data\nnot loaded",
                     ha="center", va="center", transform=ax_expr.transAxes,
                     fontsize=12, color="gray")
        ax_expr.set_title("(B)  Gene Expression (ISH)", fontweight="bold", fontsize=11)

    # ---- Panel C: Network ----
    threshold = ANALYSIS_DEFAULTS["network_edge_threshold"]
    G = nx.DiGraph()
    nodes = connectivity.index.tolist()
    G.add_nodes_from(nodes)
    edge_weights = {}
    for src in nodes:
        for tgt in nodes:
            if src != tgt:
                w = float(connectivity.loc[src, tgt])
                if w > threshold:
                    G.add_edge(src, tgt, weight=w)
                    edge_weights[(src, tgt)] = w

    pos = nx.spring_layout(G, weight="weight", seed=42, k=2.0)

    if gene and not expression_df.empty and gene in expression_df.columns:
        expr_vals  = expression_df[gene].reindex(nodes, fill_value=0.0)
        vmin, vmax = expr_vals.min(), expr_vals.max()
        vmax       = vmax if vmax > vmin else vmin + 1e-6
        expr_cmap  = cm.get_cmap("YlOrRd")
        node_colors = [expr_cmap((v - vmin) / (vmax - vmin)) for v in expr_vals]
    else:
        node_colors = [_node_color_for_region(n) for n in nodes]
        vmin, vmax  = 0, 1

    edges   = list(G.edges())
    raw_w   = [edge_weights.get(e, 0) for e in edges]
    max_w   = max(raw_w) if raw_w else 1
    e_widths = [0.8 + (w / max_w) * 4.5 for w in raw_w]
    e_colors = [cm.get_cmap("Greys")(0.3 + 0.6 * (w / max_w)) for w in raw_w]

    nx.draw_networkx_nodes(G, pos, ax=ax_net, node_color=node_colors,
                           node_size=700, alpha=0.93)
    nx.draw_networkx_labels(G, pos, ax=ax_net, font_size=5.5,
                            font_weight="bold", font_color="white")
    if edges:
        nx.draw_networkx_edges(G, pos, ax=ax_net, edge_color=e_colors,
                               width=e_widths, arrows=True, arrowsize=10,
                               connectionstyle="arc3,rad=0.15",
                               min_source_margin=15, min_target_margin=15)

    gene_str = f" — {gene} expression" if gene else ""
    ax_net.set_title(f"(C)  Network{gene_str}", fontweight="bold", fontsize=11)
    ax_net.axis("off")

    # Small colorbar for expression in network panel
    if gene and not expression_df.empty and gene in expression_df.columns:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        sm   = cm.ScalarMappable(cmap="YlOrRd", norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax_net, shrink=0.55, pad=0.02)
        cbar.set_label(f"{gene} expression", fontsize=8)

    # ---- Panel D: Hub score bar chart ----
    per_gene = hub_results.get("per_gene", {})
    if gene and gene in per_gene:
        hub_df   = per_gene[gene].head(10)
        regions  = hub_df.index.tolist()
        scores   = hub_df["hub_score"].values
        is_hubs  = hub_df["is_expression_hub"].values
        colors   = [_node_color_for_region(r) for r in regions]

        ax_hub.barh(range(len(regions)), scores, color=colors,
                    edgecolor="white", linewidth=0.5, alpha=0.88)
        ax_hub.set_yticks(range(len(regions)))
        ax_hub.set_yticklabels(
            [f"★ {r}" if h else r for r, h in zip(regions, is_hubs)],
            fontsize=8,
        )
        for tick, region in zip(ax_hub.get_yticklabels(), regions):
            tick.set_color(_node_color_for_region(region))

        ax_hub.set_xlabel("Hub Score", fontsize=9)
        ax_hub.set_title(f"(D)  Hub Rankings — {gene}", fontweight="bold", fontsize=11)
        ax_hub.spines[["top", "right"]].set_visible(False)
        ax_hub.invert_yaxis()

    elif hub_results.get("top_hubs") is not None and not hub_results["top_hubs"].empty:
        # Fallback: mean hub score across all genes
        top = hub_results["top_hubs"].head(10)
        regions = top.index.tolist()
        scores  = top["mean_hub_score"].values
        colors  = [_node_color_for_region(r) for r in regions]
        ax_hub.barh(range(len(regions)), scores, color=colors, edgecolor="white", alpha=0.88)
        ax_hub.set_yticks(range(len(regions)))
        ax_hub.set_yticklabels(regions, fontsize=8)
        ax_hub.set_xlabel("Mean Hub Score", fontsize=9)
        ax_hub.set_title("(D)  Hub Rankings (mean)", fontweight="bold", fontsize=11)
        ax_hub.spines[["top", "right"]].set_visible(False)
        ax_hub.invert_yaxis()
    else:
        ax_hub.text(0.5, 0.5, "Hub results\nnot available",
                    ha="center", va="center", transform=ax_hub.transAxes, color="gray")
        ax_hub.set_title("(D)  Hub Rankings", fontweight="bold", fontsize=11)

    _save(fig, output_path)
    return fig


# =============================================================================
# Convenience: save all figures at once
# =============================================================================

def save_all_figures(
    connectivity: pd.DataFrame,
    expression_df: pd.DataFrame,
    hub_results: dict,
    output_dir: str = "outputs",
) -> dict:
    """
    Generate and save all standard figures to output_dir.

    Parameters
    ----------
    connectivity : pd.DataFrame
    expression_df : pd.DataFrame
    hub_results : dict
        Output from run_analysis().
    output_dir : str
        Directory for saved figures. Created if absent.

    Returns
    -------
    dict
        {figure_name: file_path} mapping.
    """
    out  = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {}

    logger.info("Saving expression heatmap...")
    fig = plot_expression_heatmap(
        expression_df,
        output_path=str(out / "01_expression_heatmap.png"),
    )
    paths["expression_heatmap"] = str(out / "01_expression_heatmap.png")
    plt.close(fig)

    logger.info("Saving connectivity heatmap...")
    fig = plot_connectivity_heatmap(
        connectivity,
        output_path=str(out / "02_connectivity_heatmap.png"),
    )
    paths["connectivity_heatmap"] = str(out / "02_connectivity_heatmap.png")
    plt.close(fig)

    genes = list(hub_results.get("per_gene", {}).keys())
    for gene in genes[:3]:   # Save network for up to 3 genes
        logger.info("Saving network diagram for %s...", gene)
        expr_series = (
            expression_df[gene]
            if not expression_df.empty and gene in expression_df.columns
            else None
        )
        fig = plot_network(
            connectivity, expr_series, gene_name=gene,
            output_path=str(out / f"03_network_{gene}.png"),
        )
        paths[f"network_{gene}"] = str(out / f"03_network_{gene}.png")
        plt.close(fig)

    logger.info("Saving hub score rankings...")
    fig = plot_hub_scores(
        hub_results.get("per_gene", {}),
        output_path=str(out / "04_hub_scores.png"),
    )
    paths["hub_scores"] = str(out / "04_hub_scores.png")
    plt.close(fig)

    logger.info("Saving integrated dashboard...")
    fig = plot_integrated_dashboard(
        connectivity, expression_df, hub_results,
        gene=genes[0] if genes else None,
        output_path=str(out / "05_integrated_dashboard.png"),
    )
    paths["integrated_dashboard"] = str(out / "05_integrated_dashboard.png")
    plt.close(fig)

    logger.info("All figures saved to: %s/", output_dir)
    return paths
