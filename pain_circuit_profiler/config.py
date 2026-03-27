"""
config.py — Central configuration for PainCircuitProfiler
==========================================================
All region definitions, gene lists, and analysis parameters live here.
To customize the tool for a different circuit or gene set, edit this file.
"""

# =============================================================================
# PAIN CIRCUIT REGIONS
# =============================================================================
# Each entry: Allen Atlas acronym → display metadata
# Categories follow broad anatomical groupings used in pain literature.
# 'description' briefly states the region's role in pain processing.
#
# NOTE: Allen acronyms must match the Allen CCF ontology exactly.
#   If a region is not found, a warning is issued and it is skipped.
#   Reference: https://atlas.brain-map.org/

PAIN_REGIONS = {
    # Prefrontal / Cingulate (top-down modulation)
    "ACAd": {
        "label":       "Ant. Cingulate (dorsal)",
        "category":    "Prefrontal Cortex",
        "description": "Dorsal ACC; pain unpleasantness, attentional modulation"
    },
    "ACAv": {
        "label":       "Ant. Cingulate (ventral)",
        "category":    "Prefrontal Cortex",
        "description": "Ventral ACC; affective valence, fear extinction"
    },
    "ORBvl": {
        "label":       "Orbital Cortex (vl)",
        "category":    "Prefrontal Cortex",
        "description": "Ventrolateral orbital cortex; nociceptive integration"
    },
    # Subcortical forebrain
    "CLA": {
        "label":       "Claustrum",
        "category":    "Forebrain",
        "description": "Bilateral coordination of cortical pain representations"
    },
    "ACB": {
        "label":       "Nucleus Accumbens",
        "category":    "Forebrain",
        "description": "Pain-reward interface; motivational aspects of chronic pain"
    },
    "EP": {
        "label":       "Endopiriform Nucleus",
        "category":    "Forebrain",
        "description": "Olfactory cortex adjacent; nociceptive-affective integration"
    },
    # Amygdala complex — emotional/affective pain
    "LA": {
        "label":       "Lateral Amygdala",
        "category":    "Amygdala",
        "description": "Pain-fear conditioning; nociceptive input from thalamus"
    },
    "BLA": {
        "label":       "Basolateral Amygdala",
        "category":    "Amygdala",
        "description": "Affective pain, opioid-sensitive pain modulation"
    },
    "BMA": {
        "label":       "Basomedial Amygdala",
        "category":    "Amygdala",
        "description": "Connects BLA to hypothalamus; pain-stress interface"
    },
    "PA": {
        "label":       "Post. Amygdala",
        "category":    "Amygdala",
        "description": "Pain-related aversive memory consolidation"
    },
    # Thalamus — sensory relay
    "DORsm": {
        "label":       "Thalamus (sensorimotor)",
        "category":    "Thalamus",
        "description": "Sensorimotor thalamic group; spinothalamic relay"
    },
    "DORpm": {
        "label":       "Thalamus (polymodal)",
        "category":    "Thalamus",
        "description": "Polymodal associative thalamus; affective pain relay"
    },
    # Midbrain
    "PAG": {
        "label":       "Periaqueductal Gray",
        "category":    "Midbrain",
        "description": "Master pain modulation hub; opioid-rich descending control"
    },
    "VTA": {
        "label":       "Ventral Tegmental Area",
        "category":    "Midbrain",
        "description": "Dopaminergic; pain-reward prediction error, opioid reward"
    },
    "VTN": {
        "label":       "Ventral Tegmental Nucleus",
        "category":    "Midbrain",
        "description": "Adjacent to VTA; nociceptive processing"
    },
    "RAmb": {
        "label":       "Raphe Magnus (b)",
        "category":    "Midbrain",
        "description": "Serotonergic descending pain modulation"
    },
    # Hindbrain — ascending relays
    "PB": {
        "label":       "Parabrachial Nucleus",
        "category":    "Hindbrain",
        "description": "First major relay above spinal cord; spino-parabrachio-amygdaloid path"
    },
    "NTS": {
        "label":       "Nucleus Tractus Solitarius",
        "category":    "Hindbrain",
        "description": "Visceral nociception relay; vagal pain signaling"
    },
}

# Region acronyms grouped by category (for consistent plot ordering)
REGION_CATEGORIES = [
    "Prefrontal Cortex",
    "Forebrain",
    "Amygdala",
    "Thalamus",
    "Midbrain",
    "Hindbrain",
]

# Sorted region list: categories in order, alphabetical within each category
def get_ordered_regions():
    """Return list of region acronyms ordered by category."""
    ordered = []
    for cat in REGION_CATEGORIES:
        cat_regions = [r for r, v in PAIN_REGIONS.items() if v["category"] == cat]
        ordered.extend(sorted(cat_regions))
    # Append any remaining regions not in known categories
    remaining = [r for r in PAIN_REGIONS if r not in ordered]
    ordered.extend(remaining)
    return ordered


# =============================================================================
# TARGET GENES
# =============================================================================
# Format: acronym → metadata
# Add or remove genes here to extend the analysis.
GENES = {
    # Opioid receptors
    "Oprm1": {
        "name":     "Mu-opioid receptor",
        "function": "Primary target for morphine and endogenous endorphins",
        "pathway":  "Opioid signaling",
        "color":    "#e74c3c",    # Red
    },
    "Oprd1": {
        "name":     "Delta-opioid receptor",
        "function": "Modulates mu-receptor signaling; chronic pain plasticity",
        "pathway":  "Opioid signaling",
        "color":    "#e67e22",    # Orange
    },
    "Oprk1": {
        "name":     "Kappa-opioid receptor",
        "function": "Spinal/supraspinal analgesia; dysphoria with systemic activation",
        "pathway":  "Opioid signaling",
        "color":    "#f39c12",    # Amber
    },
    # Endogenous opioid peptides
    "Penk": {
        "name":     "Proenkephalin (Met/Leu-enkephalin precursor)",
        "function": "Endogenous ligand for delta and mu receptors",
        "pathway":  "Endogenous opioid",
        "color":    "#3498db",    # Blue
    },
    "Pdyn": {
        "name":     "Prodynorphin (dynorphin precursor)",
        "function": "Endogenous kappa-receptor ligand; kappa system in stress/pain",
        "pathway":  "Endogenous opioid",
        "color":    "#9b59b6",    # Purple
    },
    "Pomc": {
        "name":     "POMC (beta-endorphin precursor)",
        "function": "Beta-endorphin precursor; stress-induced analgesia",
        "pathway":  "Endogenous opioid",
        "color":    "#27ae60",    # Green
    },
}


# =============================================================================
# ANALYSIS PARAMETERS (defaults, all overridable at call time)
# =============================================================================
ANALYSIS_DEFAULTS = {
    # Expression threshold: regions above this percentile = "high expression"
    "expression_percentile":     75,
    # Projection threshold: regions above this percentile = "significant projector"
    "projection_percentile":     50,
    # Minimum projection_energy to draw an edge in network visualization
    "network_edge_threshold":    0.001,
    # Cache directory for downloaded Allen data
    "cache_dir":                 "allen_connectivity_cache",
    # Output directory for saved figures
    "output_dir":                "outputs",
    # Figure resolution (DPI)
    "figure_dpi":                150,
    # Colormap for connectivity heatmap
    "connectivity_cmap":         "RdYlBu_r",
    # Colormap for expression heatmap
    "expression_cmap":           "YlOrRd",
}


# =============================================================================
# CATEGORY COLORS — used consistently across all figures
# =============================================================================
CATEGORY_COLORS = {
    "Prefrontal Cortex": "#3498db",   # Blue
    "Forebrain":         "#2ecc71",   # Green
    "Amygdala":          "#9b59b6",   # Purple
    "Thalamus":          "#1abc9c",   # Teal
    "Midbrain":          "#e74c3c",   # Red
    "Hindbrain":         "#e67e22",   # Orange
}
