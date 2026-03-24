"""
Centralized cell type definitions, colors, and patterns.

All cell type metadata used across validation and visualization scripts
is defined here to avoid duplication.
"""

# ── Hippocampal cell types (used in most validation figures) ─────────────────

HIPPO_CELL_TYPES = {
    'DG Granule': {
        'pattern': '037 DG Glut',
        'color': '#e74c3c',
        'expected': 'tightly clustered in dentate gyrus',
    },
    'CA1 Pyramidal': {
        'pattern': '016 CA1-ProS Glut',
        'color': '#2ecc71',
        'expected': 'clustered in CA1 layer',
    },
    'CA3 Pyramidal': {
        'pattern': '017 CA3 Glut',
        'color': '#3498db',
        'expected': 'clustered in CA3 layer',
    },
}

# ── Validation cell types (hippocampal + cortical + negative control) ────────

VALIDATION_CELL_TYPES = {
    **HIPPO_CELL_TYPES,
    'L2/3 IT CTX': {
        'pattern': 'L2/3 IT CTX',
        'color': '#f39c12',
        'expected': 'localized to cortical L2/3',
    },
    'Oligodendrocytes': {
        'pattern': 'Oligo',
        'color': '#9b59b6',
        'expected': 'distributed throughout tissue (negative control)',
    },
}

# ── IMN cell types ───────────────────────────────────────────────────────────

IMN_CELL_TYPES = {
    'DG-PIR Ex IMN': {
        'subclass': '038 DG-PIR Ex IMN',
        'color': '#ff00ff',
        'marker': 'D',
        'edgecolor': 'white',
        's': 100,
    },
    'OB-STR-CTX Inh IMN': {
        'subclass': '045 OB-STR-CTX Inh IMN',
        'color': '#00ffff',
        'marker': '^',
        'edgecolor': 'white',
        's': 100,
    },
}

# ── Neuronal class palettes ──────────────────────────────────────────────────

GLUT_PALETTE = [
    '#e74c3c', '#ff6b35', '#f39c12', '#27ae60', '#2ecc71',
    '#1abc9c', '#e67e22', '#d4ac0d', '#a04000', '#7d3c98',
    '#2980b9', '#16a085',
]

GABA_PALETTE = [
    '#3498db', '#9b59b6', '#8e44ad', '#2c3e50', '#1a5276',
    '#5dade2', '#af7ac5', '#85929e',
]


def build_neuronal_colormap(class_names):
    """Assign colors to neuronal class names: warm for Glut, cool for GABA.

    Parameters
    ----------
    class_names : list of str
        Sorted list of neuronal class names.

    Returns
    -------
    dict
        Mapping class_name → hex color string.
    """
    glut = [c for c in class_names if 'Glut' in c]
    gaba = [c for c in class_names if 'GABA' in c or 'Gaba' in c]

    colors = {}
    for i, c in enumerate(glut):
        colors[c] = GLUT_PALETTE[i % len(GLUT_PALETTE)]
    for i, c in enumerate(gaba):
        colors[c] = GABA_PALETTE[i % len(GABA_PALETTE)]
    return colors
