# PainCircuitProfiler

**Integrating Allen Brain Atlas connectivity and gene expression data to map opioid-relevant hubs in the mouse pain circuit.**

![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Data](https://img.shields.io/badge/data-Allen%20Brain%20Atlas-orange)

---

## Overview

PainCircuitProfiler queries two public Allen Brain Atlas datasets and integrates them to answer a specific circuit-level question:

> *Which brain regions in the ascending pain pathway are best positioned to influence opioid signalling — i.e., they project strongly to opioid-receptor-rich areas AND are highly interconnected within the pain network?*

The tool produces a **hub score** for each region that combines:
1. **Molecular accessibility** — how strongly does the region project to areas with high opioid receptor/peptide expression (from Allen ISH data)?
2. **Network centrality** — how embedded is the region within the pain circuit (from Allen Connectivity Atlas data)?

Regions with high hub scores are candidates for broad opioid modulation of pain circuits — making them targets for pharmacological, optogenetic, or chemogenetic interrogation.

---

## Scientific Background

The ascending pain pathway is not a simple relay chain. It is a distributed network where nociceptive signals are simultaneously processed for sensory discrimination (thalamus → S1), affective valence (PBN → amygdala), and cognitive/attentional components (thalamus → ACC). Descending modulation from the PAG can suppress or amplify all of these signals.

Opioid drugs act throughout this network — the spatial overlap between opioid receptor expression and circuit connectivity determines where and how opioids modulate pain.

**Key regions analysed:**

| Category | Regions | Pain role |
|----------|---------|-----------|
| Prefrontal Cortex | ACAd, ACAv, vlORBm | Top-down modulation, pain affect |
| Forebrain | ACB, CLA, EP | Pain–reward interface, bilateral coordination |
| Amygdala | LA, BLA, BMA, PA | Affective/emotional pain processing |
| Thalamus | DORsm, DORpm | Spinothalamic relay, higher-order association |
| Midbrain | PAG, VTA, VTN, RAmb | Pain modulation hub, descending control |
| Hindbrain | PB, NTS | Ascending relay, visceral nociception |

**Genes analysed:**

| Gene | Protein | Pathway |
|------|---------|---------|
| Oprm1 | Mu-opioid receptor | Opioid signalling |
| Oprd1 | Delta-opioid receptor | Opioid signalling |
| Oprk1 | Kappa-opioid receptor | Opioid signalling |
| Penk | Proenkephalin | Endogenous opioid |
| Pdyn | Prodynorphin | Endogenous opioid |
| Pomc | POMC / β-endorphin | Endogenous opioid |

---

## Features

- **Allen Connectivity Atlas integration** — queries anterograde tracing experiments per region, averages across biological replicates, returns a projection energy matrix
- **Allen ISH gene expression** — queries structure-level expression energy for any gene, with local JSON caching
- **Hub metric analysis** — integrates expression and connectivity into a single ranked output per gene
- **Five publication-quality figures** — expression heatmap, connectivity heatmap, network diagram, hub score rankings, integrated 4-panel dashboard
- **Modular and extensible** — swap in custom region lists, add genes, or replace Allen expression data with your own experimental data

---

## Installation

```bash
git clone https://github.com/yourusername/PainCircuitProfiler.git
cd PainCircuitProfiler
pip install -r requirements.txt
```

> **Note:** `allensdk` requires Python ≥ 3.8. A virtual environment is recommended.

---

## Quick Start

**Option 1 — Single function call:**
```python
from pain_circuit_profiler import profile_pain_circuit

results = profile_pain_circuit()
# Downloads data, runs analysis, saves 5 figures to outputs/
```

**Option 2 — Step by step (full control):**
```python
from pain_circuit_profiler.config import PAIN_REGIONS, GENES
from pain_circuit_profiler.data.connectivity    import ConnectivityLoader
from pain_circuit_profiler.data.gene_expression import GeneExpressionLoader
from pain_circuit_profiler.analysis.hub_metrics import run_analysis, print_top_hubs
from pain_circuit_profiler.visualization.plots  import plot_integrated_dashboard

# 1. Connectivity
conn_loader  = ConnectivityLoader(cache_dir='allen_connectivity_cache')
connectivity, metadata = conn_loader.load(regions=PAIN_REGIONS)

# 2. Gene expression
acronym_to_id = conn_loader.get_structure_ids(list(PAIN_REGIONS.keys()))
id_to_label   = {v: PAIN_REGIONS[k]['label'] for k, v in acronym_to_id.items()}
expr_loader   = GeneExpressionLoader(cache_dir='allen_connectivity_cache')
expression    = expr_loader.load(list(GENES.keys()), list(acronym_to_id.values()))
expression.index = [id_to_label.get(i, str(i)) for i in expression.index]

# 3. Hub analysis
hub_results = run_analysis(connectivity, expression)
print_top_hubs(hub_results)

# 4. Visualize
fig = plot_integrated_dashboard(connectivity, expression, hub_results, gene='Oprm1')
fig.savefig('outputs/dashboard.png', dpi=150, bbox_inches='tight')
```

See `notebooks/example_analysis.ipynb` for a full walkthrough.

---

## Repository Structure

```
PainCircuitProfiler/
├── pain_circuit_profiler/
│   ├── __init__.py              # profile_pain_circuit() convenience function
│   ├── config.py                # Region definitions, gene list, analysis defaults
│   │                            # ← Edit this file to customise the analysis
│   ├── data/
│   │   ├── connectivity.py      # ConnectivityLoader — Allen Connectivity Atlas
│   │   └── gene_expression.py  # GeneExpressionLoader — Allen ISH via RMA API
│   ├── analysis/
│   │   └── hub_metrics.py       # Hub scoring pipeline; run_analysis()
│   └── visualization/
│       └── plots.py             # Five plotting functions + save_all_figures()
├── notebooks/
│   └── example_analysis.ipynb  # Full walkthrough with markdown explanations
├── outputs/                     # Generated figures (created on first run)
├── allen_connectivity_cache/    # Downloaded Allen data (created on first run)
├── requirements.txt
└── README.md
```

---

## Output Figures

| Filename | Description |
|----------|-------------|
| `01_expression_heatmap.png` | ISH expression energy: regions × genes |
| `02_connectivity_heatmap.png` | Projection energy matrix (log₁₀ scale) |
| `03_network_<gene>.png` | Directed network; node colour = expression |
| `04_hub_scores.png` | Ranked bar charts per gene |
| `05_integrated_dashboard.png` | 4-panel summary figure |

---

## Customising the Analysis

### Use a different region set
Edit `PAIN_REGIONS` in `config.py`, or pass a custom dict at runtime:
```python
my_regions = {
    'NTS': {'label': 'NTS', 'category': 'Hindbrain', 'description': '...'},
    'RVM': {'label': 'RVM', 'category': 'Hindbrain', 'description': '...'},
}
results = profile_pain_circuit(regions=my_regions)
```

### Query different genes
```python
results = profile_pain_circuit(genes=['Grm5', 'Cnr1', 'Trpv1'])
```

### Replace Allen expression data with your own
`run_analysis()` accepts any DataFrame with region labels as index and gene/condition names as columns. You can substitute patch-seq data, TRAP2 profiles, or bulk RNA-seq:
```python
import pandas as pd
my_data = pd.read_csv('my_trap2_data.csv', index_col=0)  # regions × conditions
hub_results = run_analysis(connectivity, my_data)
```

### Restrict to one hemisphere
```python
connectivity_ipsi, _ = conn_loader.load(PAIN_REGIONS, hemisphere_ids=[2])
```

---

## Data Sources

- **Connectivity:** Oh, S.W. et al. (2014). A mesoscale connectome of the mouse brain. *Nature*, 508, 207–214. https://connectivity.brain-map.org
- **Expression:** Lein, E.S. et al. (2007). Genome-wide atlas of gene expression in the adult mouse brain. *Nature*, 445, 168–176. https://mouse.brain-map.org
- **Allen SDK:** Copyright 2015–2024 Allen Institute for Brain Science. https://github.com/AllenInstitute/AllenSDK

---

## License

MIT License. Data from the Allen Brain Atlas is subject to the [Allen Institute Terms of Use](https://alleninstitute.org/terms-of-use/).
