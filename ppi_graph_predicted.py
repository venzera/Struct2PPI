#!/usr/bin/env python3
"""
Struct2PPI -- Predicted-mode PPI graph (LIVIA / AFM-LIS integration)
====================================================================
Build an interactive protein-protein interaction graph from a *predicted*
complex (AlphaFold3, Boltz, Chai-1) and score every chain pair with the
Local Interaction Score framework (LIS / cLIS / iLIS, LIA / cLIA).

This is the structure-prediction counterpart to ``ppi_graph_3d_dg.py``:
instead of PRODIGY/FoldX binding energy from an experimental structure, edges
carry confidence-based interaction metrics from the predictor's PAE + pLDDT.

Clicking an edge opens a panel with, for that chain-chain interaction:
  * the 3D pair complex (two chains colored)
  * the inter-chain PAE interaction plot
  * the contact-filtered LIA (cLIA) map
  * LIS / cLIS / iLIS, LIA / cLIA / iLIA, ipTM and pLDDT scores

One prediction -> one graph. Uses the best model (highest mean iLIS) unless
--model is given. LIS scoring is provided by the bundled ``lis_lib`` (a trimmed
adaptation of AFM-LIS: no parallelism, no archive handling, no ChimeraX export).

Reference:
    Kim, Ah-Ram, and Norbert Perrimon. "LIVIA: a browser-based tool for
    assessing and visualizing predicted protein interactions."
    bioRxiv (2026). doi:10.64898/2026.05.01.721633
    iLIS metric: Kim A.-R. et al. "FlyPredictome." bioRxiv (2026).
    doi:10.64898/2026.04.14.718529

Usage:
    python ppi_graph_predicted.py /path/to/af3_output_folder
    python ppi_graph_predicted.py /path/to/boltz_folder --model 0
"""

import argparse
import io
import json
import os
import tempfile
import warnings

import networkx as nx
import numpy as np

from Bio.PDB import PDBParser, PDBIO, Select, is_aa
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import gemmi

import lis_lib

warnings.filterwarnings('ignore', category=PDBConstructionWarning)

EDGE_COLORS = [
    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
    '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
]


# ============================================================================
# Structure -> per-chain PDB strings (for 3Dmol rendering)
# ============================================================================

class _ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == ' ' and is_aa(residue, standard=True)


def load_biopython_structure(struct_path, fmt):
    """Return a BioPython structure, converting CIF -> PDB via gemmi if needed."""
    if fmt == 'cif':
        st = gemmi.read_structure(struct_path)
        st.setup_entities()
        tmp = tempfile.mktemp(suffix='.pdb')
        st.write_pdb(tmp)
        struct = PDBParser(QUIET=True).get_structure('pred', tmp)
        try:
            os.remove(tmp)
        except OSError:
            pass
        return struct
    return PDBParser(QUIET=True).get_structure('pred', struct_path)


def chain_pdb_string(structure, chain_id):
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    out = io.StringIO()
    pdb_io.save(out, _ChainSelect(chain_id))
    return out.getvalue()


# ============================================================================
# Model selection
# ============================================================================

def _mean_ilis(model):
    pairs = model['pairs']
    if not pairs:
        return 0.0
    return sum(p['iLIS'] for p in pairs) / len(pairs)


def select_model(models, requested=None):
    """Pick the model to visualize. requested matches rank string; else best iLIS."""
    if requested is not None:
        for m in models:
            if str(m['rank']) == str(requested):
                return m
        raise SystemExit(f"Model rank '{requested}' not found. "
                         f"Available: {[m['rank'] for m in models]}")
    return max(models, key=_mean_ilis)


# ============================================================================
# Graph + JSON payload
# ============================================================================

def get_chain_sequences(struct_path):
    """One-letter sequence per chain (positions align with LIS residue indices)."""
    st = gemmi.read_structure(struct_path)
    st.setup_entities()
    seqs = {}
    for chain in st[0]:
        poly = chain.get_polymer()
        if poly:
            seqs[chain.name] = gemmi.one_letter_code(poly.extract_sequence()).upper()
    return seqs


def build_payload(model, edge_min_ilis=0.0):
    """Build node/edge JSON for the HTML from a single analyzed model."""
    structure = load_biopython_structure(model['struct_path'], model['fmt'])
    bio_chains = {c.id for c in structure[0]}
    seqs = get_chain_sequences(model['struct_path'])

    chains = model['chains']
    plddt = model['chain_plddt']

    # Graph for layout
    G = nx.Graph()
    for c in chains:
        G.add_node(c)
    edge_pairs = [p for p in model['pairs'] if p['iLIS'] >= edge_min_ilis]
    for p in edge_pairs:
        G.add_edge(p['ci'], p['cj'], weight=p['iLIS'])

    if len(G.nodes()) <= 1:
        pos = {n: (0.0, 0.0) for n in G.nodes()}
    else:
        pos = nx.spring_layout(G, k=3, iterations=100, seed=42)

    xs = [p[0] for p in pos.values()] or [0]
    ys = [p[1] for p in pos.values()] or [0]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    xr = (x_max - x_min) or 1
    yr = (y_max - y_min) or 1
    npos = {n: (100 + (x - x_min) / xr * 1200, 100 + (y - y_min) / yr * 800)
            for n, (x, y) in pos.items()}

    nodes_js = []
    for c in chains:
        x, y = npos.get(c, (600, 450))
        pdb = ''
        if c in bio_chains:
            pdb = chain_pdb_string(structure, c)
        nodes_js.append({
            'id': c,
            'x': x, 'y': y,
            'label': c,
            'plddt': None if plddt.get(c) is None else round(plddt[c], 1),
            'degree': G.degree(c) if c in G else 0,
            'pdb': pdb,
        })

    # Sort edges by iLIS (strongest first) for stable coloring/legend
    edge_pairs = sorted(edge_pairs, key=lambda p: -p['iLIS'])
    edges_js = []
    for idx, p in enumerate(edge_pairs):
        x1, y1 = npos[p['ci']]
        x2, y2 = npos[p['cj']]
        positive = p['iLIS'] >= lis_lib.ILIS_THRESHOLD
        edges_js.append({
            'name': f"{p['ci']}-{p['cj']}",
            'chain1': p['ci'], 'chain2': p['cj'],
            'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2,
            'color': EDGE_COLORS[idx % len(EDGE_COLORS)],
            'width': float(min(2 + 8 * p['iLIS'], 10)),
            'positive': positive,
            'iLIS': round(p['iLIS'], 3),
            'LIS': round(p['LIS'], 3),
            'cLIS': round(p['cLIS'], 3),
            'iLIA': round(p['iLIA'], 1),
            'LIA': int(p['LIA']),
            'cLIA': int(p['cLIA']),
            'ipTM': round(p['ipTM'], 3) if p['ipTM'] is not None else None,
            'plddt_i': None if p.get('pLDDT_i') is None else round(p['pLDDT_i'], 1),
            'plddt_j': None if p.get('pLDDT_j') is None else round(p['pLDDT_j'], 1),
            'liplddt_i': _fmt(p.get('LIpLDDT_i')),
            'liplddt_j': _fmt(p.get('LIpLDDT_j')),
            'cliplddt_i': _fmt(p.get('cLIpLDDT_i')),
            'cliplddt_j': _fmt(p.get('cLIpLDDT_j')),
            'pae_map': p['pae_map'],
            'clia_map': p['clia_map'],
            'res_i': p['res_i'],
            'res_j': p['res_j'],
            'seq_i': seqs.get(p['ci'], ''),
            'seq_j': seqs.get(p['cj'], ''),
            'lir_i': [int(x) for x in sorted(p['lirI'])],
            'lir_j': [int(x) for x in sorted(p['lirJ'])],
            'clir_i': [int(x) for x in sorted(p['clirI'])],
            'clir_j': [int(x) for x in sorted(p['clirJ'])],
        })

    return nodes_js, edges_js, G


def _fmt(v):
    if v is None:
        return None
    try:
        if v != v:  # NaN
            return None
        return round(float(v), 1)
    except (TypeError, ValueError):
        return None


def _normalize_layout(G):
    """Spring layout -> pixel coordinates in the 1400x1000 canvas."""
    if len(G.nodes()) <= 1:
        pos = {n: (0.0, 0.0) for n in G.nodes()}
    else:
        pos = nx.spring_layout(G, k=3, iterations=100, seed=42)
    xs = [p[0] for p in pos.values()] or [0]
    ys = [p[1] for p in pos.values()] or [0]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    xr = (x_max - x_min) or 1
    yr = (y_max - y_min) or 1
    return {n: (100 + (x - x_min) / xr * 1200, 100 + (y - y_min) / yr * 800)
            for n, (x, y) in pos.items()}


def find_structure_file(path):
    """Return (struct_path, fmt) for a file or the best structure inside a folder."""
    if os.path.isfile(path):
        return path, ('cif' if path.endswith('.cif') else 'pdb')
    cands = []
    for root, _dirs, files in os.walk(path):
        for f in files:
            if f.endswith(('.cif', '.pdb')) and not f.startswith('._'):
                cands.append(os.path.join(root, f))
    if not cands:
        return None, None
    cands.sort(key=lambda p: (0 if ('model_0' in p or 'rank_0' in p) else 1, p))
    return cands[0], ('cif' if cands[0].endswith('.cif') else 'pdb')


def detect_contacts(structure, cutoff=5.0):
    """Distance-based residue contacts between protein chains (heavy atoms)."""
    from scipy.spatial.distance import cdist
    chains = {}
    for chain in structure[0]:
        coords, resids = [], []
        for res in chain:
            if res.id[0] != ' ' or not is_aa(res, standard=True):
                continue
            for atom in res:
                if atom.element == 'H':
                    continue
                coords.append(atom.coord)
                resids.append(res.id[1])
        if coords:
            chains[chain.id] = (np.array(coords), np.array(resids))
    ids = list(chains)
    contacts = {}
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            ca, ra = chains[ids[a]]
            cb, rb = chains[ids[b]]
            d = cdist(ca, cb)
            ai, bi = np.where(d < cutoff)
            if len(ai):
                pairs = set(zip(ra[ai].tolist(), rb[bi].tolist()))
                contacts[(ids[a], ids[b])] = len(pairs)
    return ids, contacts


def build_structure_only_payload(struct_path, fmt, cutoff=5.0):
    """Node/edge JSON when no PAE is available: contacts + pLDDT, no LIS."""
    structure = load_biopython_structure(struct_path, fmt)
    with open(struct_path) as fh:
        struct_text = fh.read()
    plddt = lis_lib.compute_chain_plddt(struct_text, fmt)
    chain_ids, contacts = detect_contacts(structure, cutoff)

    G = nx.Graph()
    for c in chain_ids:
        G.add_node(c)
    for (c1, c2), n in contacts.items():
        G.add_edge(c1, c2, weight=n)
    npos = _normalize_layout(G)

    nodes_js = []
    for c in chain_ids:
        x, y = npos.get(c, (600, 450))
        nodes_js.append({
            'id': c, 'x': x, 'y': y, 'label': c,
            'plddt': None if plddt.get(c) is None else round(plddt[c], 1),
            'degree': G.degree(c) if c in G else 0,
            'pdb': chain_pdb_string(structure, c),
        })

    items = sorted(contacts.items(), key=lambda kv: -kv[1])
    edges_js = []
    for idx, ((c1, c2), n) in enumerate(items):
        x1, y1 = npos[c1]
        x2, y2 = npos[c2]
        edges_js.append({
            'name': f"{c1}-{c2}", 'chain1': c1, 'chain2': c2,
            'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2,
            'color': EDGE_COLORS[idx % len(EDGE_COLORS)],
            'width': float(min(2 + np.log1p(n), 8)),
            'positive': True, 'contacts': int(n),
        })
    return nodes_js, edges_js


# ============================================================================
# HTML rendering
# ============================================================================

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html>
<head>
    <title>Predicted PPI Graph: __TITLE__</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }
        h1 { text-align: center; color: #333; margin-bottom: 4px; }
        #mode-banner {
            text-align: center; margin: 0 auto 12px auto; max-width: 1400px;
            background: #2c3e50; color: #fff; padding: 8px 14px; border-radius: 6px;
            font-size: 13px; letter-spacing: 0.3px;
        }
        #mode-banner b { color: #ffd54f; }
        #container {
            position: relative; width: 1400px; height: 1000px; margin: 0 auto;
            background: white; border: 1px solid #ddd; border-radius: 8px; overflow: hidden;
        }
        #graph-canvas { position: absolute; top: 0; left: 0; width: 100%; height: 100%; z-index: 1; }
        #graph-canvas.edge-hover { cursor: pointer; }
        .node-container {
            position: absolute; width: 100px; height: 100px; border-radius: 50%;
            border: 3px solid #333; background: white; cursor: move; z-index: 10;
            overflow: hidden; box-shadow: 0 2px 10px rgba(0,0,0,0.2);
        }
        .node-container:hover { border-color: #0066cc; box-shadow: 0 4px 20px rgba(0,102,204,0.4); }
        .node-viewer { width: 100%; height: 100%; border-radius: 50%; }
        .node-label {
            position: absolute; bottom: -25px; left: 50%; transform: translateX(-50%);
            font-size: 12px; font-weight: bold; white-space: nowrap;
            background: rgba(255,255,255,0.9); padding: 2px 6px; border-radius: 3px;
        }
        #legend {
            position: absolute; top: 10px; right: 10px; background: rgba(255,255,255,0.95);
            border: 1px solid #ddd; border-radius: 5px; padding: 10px; max-height: 900px;
            overflow-y: auto; z-index: 100; font-size: 11px; max-width: 280px;
        }
        #legend h3 { margin: 0 0 10px 0; font-size: 14px; }
        .legend-item { display: flex; align-items: center; margin: 3px 0; cursor: pointer; padding: 2px; }
        .legend-item:hover { background: #f0f0f0; }
        .legend-color { width: 20px; height: 4px; margin-right: 8px; flex-shrink: 0; border-radius: 2px; }
        .legend-text { display: flex; flex-direction: column; }
        .legend-name { font-weight: bold; }
        .legend-ilis { color: #0066cc; font-size: 10px; }
        .pos-dot { color: #2e7d32; font-weight: bold; }
        #info { text-align: center; margin-top: 10px; color: #666; }
        #tooltip {
            position: absolute; background: rgba(0,0,0,0.85); color: white; padding: 8px 12px;
            border-radius: 4px; font-size: 12px; pointer-events: none; z-index: 1000;
            display: none; max-width: 320px;
        }
        .edge-label {
            position: absolute; font-size: 10px; font-weight: bold; color: #0066cc;
            background: rgba(255,255,255,0.9); padding: 1px 4px; border-radius: 3px;
            z-index: 5; pointer-events: none;
        }
        #modal {
            display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%;
            background: rgba(0,0,0,0.6); z-index: 2000; justify-content: center; align-items: center;
        }
        #modal.visible { display: flex; }
        #modal-content {
            background: white; border-radius: 8px; width: 1080px; max-width: 96%;
            height: 760px; max-height: 96%; position: relative; display: flex;
            flex-direction: column; overflow: hidden; box-shadow: 0 10px 40px rgba(0,0,0,0.4);
        }
        #modal-header {
            padding: 12px 20px; border-bottom: 1px solid #ddd; background: #f9f9f9;
            display: flex; justify-content: space-between; align-items: center;
        }
        #modal-title { font-weight: bold; font-size: 16px; color: #333; }
        #modal-close {
            background: #e6194b; color: white; border: none; border-radius: 4px;
            width: 30px; height: 30px; cursor: pointer; font-size: 18px; font-weight: bold; line-height: 1;
        }
        #modal-close:hover { background: #c4003a; }
        #modal-body { flex: 1; display: flex; overflow: hidden; }
        #modal-left { width: 420px; display: flex; flex-direction: column; border-right: 1px solid #eee; }
        #viewer3d { flex: 1; position: relative; background: white; min-height: 200px; }
        #score-table { padding: 10px 14px; font-size: 12px; }
        #score-table table { width: 100%; border-collapse: collapse; }
        #score-table td { padding: 3px 6px; border-bottom: 1px solid #f0f0f0; }
        #score-table td.k { color: #666; }
        #score-table td.v { text-align: right; font-weight: bold; color: #222; }
        .verdict { margin: 6px 14px; padding: 6px 10px; border-radius: 5px; font-size: 12px; font-weight: bold; }
        .verdict.pos { background: #e8f5e9; color: #2e7d32; }
        .verdict.neg { background: #fdecea; color: #b71c1c; }
        #modal-right { flex: 1; display: flex; flex-direction: column; padding: 10px; overflow-y: auto; }
        .map-block { margin-bottom: 14px; }
        .map-block h4 { margin: 4px 0; font-size: 13px; color: #333; }
        .map-block .axis { font-size: 11px; color: #888; }
        .plot { width: 100%; height: 300px; }
        .colorbar { font-size: 10px; color: #666; margin-top: 2px; }
        .swatch { display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 4px; vertical-align: middle; }
        .seqtrack { font-family: 'Courier New', monospace; font-size: 13px; line-height: 1.6;
            word-break: break-all; margin: 4px 0; padding: 4px; background: #fafafa; border-radius: 4px; }
        .seqlabel { font-weight: bold; margin-right: 6px; }
        .res { padding: 0 0.5px; border-radius: 2px; }
        .res.lir { background: #fff3b0; }
        .res.clir { background: #fb8c00; color: #fff; font-weight: bold; }
        .seqlegend { font-size: 10px; color: #666; margin-bottom: 2px; }
    </style>
</head>
<body>
    <h1>Predicted PPI Graph: __TITLE__</h1>
    <div id="mode-banner">__BANNER__</div>
    <div id="container">
        <canvas id="graph-canvas"></canvas>
        <div id="legend">
            <h3>__LEGEND_TITLE__</h3>
            <div id="legend-items"></div>
        </div>
    </div>
    <div id="info">
        __NNODES__ chains | __NEDGES__ chain pairs | Drag nodes to rearrange |
        Double-click a node for its chain + pLDDT | __INFO_HINT__
    </div>
    <div id="tooltip"></div>

    <div id="modal">
        <div id="modal-content">
            <div id="modal-header">
                <div id="modal-title">Interaction</div>
                <button id="modal-close" title="Close">&times;</button>
            </div>
            <div id="modal-body">
                <div id="modal-left">
                    <div id="viewer3d"></div>
                    <div class="verdict" id="verdict"></div>
                    <div id="score-table"></div>
                </div>
                <div id="modal-right">
                    <div class="map-block">
                        <h4>PAE interaction plot (inter-chain)</h4>
                        <div class="axis" id="pae-axis"></div>
                        <div class="plot" id="pae-plot"></div>
                    </div>
                    <div class="map-block">
                        <h4>Contact-filtered LIA map (cLIA)</h4>
                        <div class="axis" id="clia-axis"></div>
                        <div class="plot" id="clia-plot"></div>
                        <div class="colorbar"><span class="swatch" style="background:#08306b"></span>confident contact (PAE&le;12 &amp; C&beta;&le;8&Aring;)
                            <span class="swatch" style="background:#f7fbff;border:1px solid #ccc"></span>none</div>
                    </div>
                    <div class="map-block" id="seq-block">
                        <h4>Interface sequences</h4>
                        <div class="seqlegend">
                            <span class="swatch" style="background:#fb8c00"></span>cLIR (confident contact residue)
                            <span class="swatch" style="background:#fff3b0;border:1px solid #ccc"></span>LIR (confident interface, PAE&le;12)
                        </div>
                        <div class="seqtrack" id="seq-i"></div>
                        <div class="seqtrack" id="seq-j"></div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <script>
        const nodes = __NODES_JSON__;
        const edges = __EDGES_JSON__;
        const ILIS_THRESHOLD = __THRESH__;
        const HAS_SCORES = __HAS_SCORES__;
        const sortedEdges = [...edges];

        const container = document.getElementById('container');
        const canvas = document.getElementById('graph-canvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');
        canvas.width = 1400; canvas.height = 1000;

        const nodeElements = {};
        const viewers = {};
        const edgeLabelElements = {};
        let draggedNode = null;
        let dragOffset = {x: 0, y: 0};

        function drawEdges() {
            ctx.clearRect(0, 0, canvas.width, canvas.height);
            edges.forEach(edge => {
                const n1 = nodes.find(n => n.id === edge.chain1);
                const n2 = nodes.find(n => n.id === edge.chain2);
                if (n1 && n2) {
                    ctx.beginPath();
                    ctx.moveTo(n1.x + 50, n1.y + 50);
                    ctx.lineTo(n2.x + 50, n2.y + 50);
                    ctx.strokeStyle = edge.color;
                    ctx.lineWidth = edge.width;
                    ctx.setLineDash(edge.positive ? [] : [6, 5]);
                    ctx.stroke();
                    ctx.setLineDash([]);
                    const labelEl = edgeLabelElements[edge.name];
                    if (labelEl) {
                        labelEl.style.left = ((n1.x + n2.x) / 2 + 50) + 'px';
                        labelEl.style.top = ((n1.y + n2.y) / 2 + 50) + 'px';
                    }
                }
            });
        }

        function createEdgeLabels() {
            edges.forEach(edge => {
                const n1 = nodes.find(n => n.id === edge.chain1);
                const n2 = nodes.find(n => n.id === edge.chain2);
                if (n1 && n2) {
                    const el = document.createElement('div');
                    el.className = 'edge-label';
                    el.textContent = HAS_SCORES ? ('iLIS ' + edge.iLIS.toFixed(3))
                                                : (edge.contacts + ' contacts');
                    el.style.left = ((n1.x + n2.x) / 2 + 50) + 'px';
                    el.style.top = ((n1.y + n2.y) / 2 + 50) + 'px';
                    container.appendChild(el);
                    edgeLabelElements[edge.name] = el;
                }
            });
        }

        function createNodes() {
            nodes.forEach(node => {
                const div = document.createElement('div');
                div.className = 'node-container';
                div.id = 'node-' + node.id;
                div.style.left = node.x + 'px';
                div.style.top = node.y + 'px';
                const viewer = document.createElement('div');
                viewer.className = 'node-viewer';
                viewer.id = 'viewer-' + node.id;
                div.appendChild(viewer);
                const label = document.createElement('div');
                label.className = 'node-label';
                label.textContent = node.id + (node.plddt != null ? ' (pLDDT ' + node.plddt + ')' : '');
                div.appendChild(label);

                div.addEventListener('mousedown', (e) => {
                    if (e.target === div || e.target === viewer) {
                        draggedNode = node;
                        const rect = div.getBoundingClientRect();
                        dragOffset.x = e.clientX - rect.left;
                        dragOffset.y = e.clientY - rect.top;
                        div.style.zIndex = 100;
                    }
                });
                div.addEventListener('dblclick', (e) => {
                    e.preventDefault(); e.stopPropagation();
                    showSingleChain(node);
                });
                div.addEventListener('mouseenter', () => {
                    tooltip.innerHTML = '<b>Chain ' + node.id + '</b><br>' +
                        (node.plddt != null ? 'avg pLDDT: ' + node.plddt + '<br>' : '') +
                        'Interactions: ' + node.degree;
                    tooltip.style.display = 'block';
                });
                div.addEventListener('mousemove', (e) => {
                    tooltip.style.left = (e.pageX + 15) + 'px';
                    tooltip.style.top = (e.pageY + 15) + 'px';
                });
                div.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
                container.appendChild(div);
                nodeElements[node.id] = div;
            });
        }

        function initViewers() {
            nodes.forEach(node => {
                const el = document.getElementById('viewer-' + node.id);
                if (el && node.pdb) {
                    const v = $3Dmol.createViewer(el, { backgroundColor: 'white' });
                    v.addModel(node.pdb, 'pdb');
                    v.setStyle({}, {cartoon: {color: 'spectrum'}});
                    v.zoomTo(); v.zoom(0.8); v.render();
                    viewers[node.id] = v;
                }
            });
        }

        function createLegend() {
            const box = document.getElementById('legend-items');
            sortedEdges.forEach(edge => {
                const item = document.createElement('div');
                item.className = 'legend-item';
                const detail = HAS_SCORES
                    ? ('iLIS ' + edge.iLIS.toFixed(3) + ' | LIS ' + edge.LIS.toFixed(3) +
                       ' | cLIS ' + edge.cLIS.toFixed(3))
                    : (edge.contacts + ' residue contacts');
                item.innerHTML = '<div class="legend-color" style="background:' + edge.color + '"></div>' +
                    '<div class="legend-text"><span class="legend-name">' + edge.name +
                    (HAS_SCORES && edge.positive ? ' <span class="pos-dot">&#10003;</span>' : '') + '</span>' +
                    '<span class="legend-ilis">' + detail + '</span></div>';
                item.addEventListener('click', () => showInteraction(edge));
                box.appendChild(item);
            });
        }

        document.addEventListener('mousemove', (e) => {
            if (draggedNode) {
                const r = container.getBoundingClientRect();
                let nx = e.clientX - r.left - dragOffset.x;
                let ny = e.clientY - r.top - dragOffset.y;
                nx = Math.max(0, Math.min(nx, container.offsetWidth - 100));
                ny = Math.max(0, Math.min(ny, container.offsetHeight - 100));
                draggedNode.x = nx; draggedNode.y = ny;
                const div = nodeElements[draggedNode.id];
                div.style.left = nx + 'px'; div.style.top = ny + 'px';
                drawEdges();
            }
        });
        document.addEventListener('mouseup', () => {
            if (draggedNode) { nodeElements[draggedNode.id].style.zIndex = 10; draggedNode = null; }
        });

        createNodes(); createEdgeLabels(); drawEdges(); createLegend();
        setTimeout(initViewers, 100);

        // ===== Plotly heatmaps (ColabFold style) =====
        // bwr colorscale: blue (confident, 0) -> white -> red (uncertain, 30)
        const BWR = [[0.0, 'rgb(0,0,255)'], [0.5, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,0)']];
        const CLIA_SCALE = [[0.0, '#f7fbff'], [1.0, '#08306b']];

        function plotPAE(divId, edge) {
            const div = document.getElementById(divId);
            if (!edge || !edge.pae_map || !edge.pae_map.length) { Plotly.purge(div); return; }
            const data = [{
                z: edge.pae_map, x: edge.res_j, y: edge.res_i, type: 'heatmap',
                colorscale: BWR, zmin: 0, zmax: 30, zsmooth: false,
                colorbar: { title: 'PAE (Å)', titleside: 'right', thickness: 12, len: 0.9 },
                hovertemplate: edge.chain1 + ' %{y} · ' + edge.chain2 + ' %{x}<br>PAE %{z:.1f} Å<extra></extra>'
            }];
            const layout = {
                margin: { l: 44, r: 10, t: 6, b: 36 },
                xaxis: { title: 'chain ' + edge.chain2, constrain: 'domain' },
                yaxis: { title: 'chain ' + edge.chain1, autorange: 'reversed', scaleanchor: 'x', constrain: 'domain' },
                font: { size: 10 }
            };
            Plotly.react(div, data, layout, { displayModeBar: false, responsive: true });
        }

        function plotCLIA(divId, edge) {
            const div = document.getElementById(divId);
            if (!edge || !edge.clia_map || !edge.clia_map.length) { Plotly.purge(div); return; }
            const data = [{
                z: edge.clia_map, x: edge.res_j, y: edge.res_i, type: 'heatmap',
                colorscale: CLIA_SCALE, zmin: 0, zmax: 1, showscale: false,
                hovertemplate: edge.chain1 + ' %{y} · ' + edge.chain2 + ' %{x}<br>%{z}<extra></extra>'
            }];
            const layout = {
                margin: { l: 44, r: 10, t: 6, b: 36 },
                xaxis: { title: 'chain ' + edge.chain2, constrain: 'domain' },
                yaxis: { title: 'chain ' + edge.chain1, autorange: 'reversed', scaleanchor: 'x', constrain: 'domain' },
                font: { size: 10 }
            };
            Plotly.react(div, data, layout, { displayModeBar: false, responsive: true });
        }

        function clearPlots() {
            Plotly.purge(document.getElementById('pae-plot'));
            Plotly.purge(document.getElementById('clia-plot'));
        }

        function renderSequence(elId, chainId, seq, lir, clir, color) {
            const el = document.getElementById(elId);
            if (!seq) { el.innerHTML = ''; return; }
            const L = new Set(lir), C = new Set(clir);
            let html = '<span class="seqlabel" style="color:' + color + '">' + chainId + '</span>';
            for (let i = 0; i < seq.length; i++) {
                const pos = i + 1;
                let cls = 'res';
                if (C.has(pos)) cls += ' clir';
                else if (L.has(pos)) cls += ' lir';
                html += '<span class="' + cls + '" title="' + chainId + ' ' + pos + ' ' + seq[i] + '">' + seq[i] + '</span>';
            }
            el.innerHTML = html;
        }

        // ===== Modal =====
        const modal = document.getElementById('modal');
        const modalTitle = document.getElementById('modal-title');
        const viewer3dEl = document.getElementById('viewer3d');
        const scoreTable = document.getElementById('score-table');
        const verdictEl = document.getElementById('verdict');
        const paeAxis = document.getElementById('pae-axis');
        const cliaAxis = document.getElementById('clia-axis');
        let modalViewer = null;

        function row(k, v) { return '<tr><td class="k">' + k + '</td><td class="v">' + v + '</td></tr>'; }
        function fmt(v, d) { return (v == null) ? '&ndash;' : (typeof v === 'number' ? v.toFixed(d) : v); }

        const modalRight = document.getElementById('modal-right');

        function open3DPair(edge) {
            const colorA = '#4363d8', colorB = '#e6194b';
            viewer3dEl.innerHTML = '';
            modalViewer = $3Dmol.createViewer(viewer3dEl, { backgroundColor: 'white' });
            const nA = nodes.find(n => n.id === edge.chain1);
            const nB = nodes.find(n => n.id === edge.chain2);
            if (nA && nA.pdb) modalViewer.addModel(nA.pdb, 'pdb');
            if (nB && nB.pdb) modalViewer.addModel(nB.pdb, 'pdb');
            modalViewer.setStyle({chain: edge.chain1}, {cartoon: {color: colorA}});
            modalViewer.setStyle({chain: edge.chain2}, {cartoon: {color: colorB}});
            modalViewer.zoomTo(); modalViewer.render();
            setTimeout(() => { if (modalViewer) modalViewer.resize(); }, 60);
        }

        function showInteractionNoScores(edge) {
            modalTitle.innerHTML = 'Interaction ' + edge.chain1 + ' &ndash; ' + edge.chain2;
            verdictEl.className = 'verdict neg';
            verdictEl.innerHTML = 'Score files (PAE) were not provided &mdash; LIS metrics unavailable.';
            scoreTable.innerHTML = '<table>' +
                row('Chain pair', edge.chain1 + ' &ndash; ' + edge.chain2) +
                row('Residue contacts', edge.contacts) + '</table>';
            modalRight.style.display = 'none';
            clearPlots();
            modal.classList.add('visible');
            open3DPair(edge);
        }

        function showInteraction(edge) {
            if (!HAS_SCORES) { showInteractionNoScores(edge); return; }
            modalRight.style.display = 'flex';
            const colorA = '#4363d8', colorB = '#e6194b';
            modalTitle.innerHTML = 'Interaction ' + edge.chain1 + ' &ndash; ' + edge.chain2 +
                ' &nbsp;<span style="color:#0066cc">iLIS ' + edge.iLIS.toFixed(3) + '</span>';

            verdictEl.className = 'verdict ' + (edge.positive ? 'pos' : 'neg');
            verdictEl.innerHTML = edge.positive
                ? '&#10003; Likely interaction (iLIS &ge; ' + ILIS_THRESHOLD + ')'
                : '&#10007; Below iLIS threshold (' + ILIS_THRESHOLD + ')';

            scoreTable.innerHTML = '<table>' +
                row('iLIS = &radic;(LIS&middot;cLIS)', '<b style="color:#0066cc">' + edge.iLIS.toFixed(3) + '</b>') +
                row('LIS (PAE&le;12)', edge.LIS.toFixed(3)) +
                row('cLIS (PAE&le;12 &amp; C&beta;&le;8)', edge.cLIS.toFixed(3)) +
                row('LIA / cLIA', edge.LIA + ' / ' + edge.cLIA) +
                row('iLIA', edge.iLIA.toFixed(1)) +
                row('ipTM', fmt(edge.ipTM, 3)) +
                row('avg pLDDT ' + edge.chain1 + ' / ' + edge.chain2,
                    fmt(edge.plddt_i, 1) + ' / ' + fmt(edge.plddt_j, 1)) +
                row('LI-pLDDT ' + edge.chain1 + ' / ' + edge.chain2,
                    fmt(edge.liplddt_i, 1) + ' / ' + fmt(edge.liplddt_j, 1)) +
                row('cLI-pLDDT ' + edge.chain1 + ' / ' + edge.chain2,
                    fmt(edge.cliplddt_i, 1) + ' / ' + fmt(edge.cliplddt_j, 1)) +
                '</table>';

            paeAxis.innerHTML = 'rows = chain <b>' + edge.chain1 + '</b> (' + edge.res_i.length +
                ' res), cols = chain <b>' + edge.chain2 + '</b> (' + edge.res_j.length + ' res)';
            cliaAxis.innerHTML = paeAxis.innerHTML;

            modal.classList.add('visible');
            plotPAE('pae-plot', edge);
            plotCLIA('clia-plot', edge);
            renderSequence('seq-i', edge.chain1, edge.seq_i, edge.lir_i, edge.clir_i, '#4363d8');
            renderSequence('seq-j', edge.chain2, edge.seq_j, edge.lir_j, edge.clir_j, '#e6194b');

            viewer3dEl.innerHTML = '';
            modalViewer = $3Dmol.createViewer(viewer3dEl, { backgroundColor: 'white' });
            const nA = nodes.find(n => n.id === edge.chain1);
            const nB = nodes.find(n => n.id === edge.chain2);
            if (nA && nA.pdb) modalViewer.addModel(nA.pdb, 'pdb');
            if (nB && nB.pdb) modalViewer.addModel(nB.pdb, 'pdb');
            modalViewer.setStyle({chain: edge.chain1}, {cartoon: {color: colorA}});
            modalViewer.setStyle({chain: edge.chain2}, {cartoon: {color: colorB}});
            modalViewer.zoomTo(); modalViewer.render();
            setTimeout(() => { if (modalViewer) modalViewer.resize(); }, 60);
        }

        function showSingleChain(node) {
            modalTitle.innerHTML = 'Chain ' + node.id;
            verdictEl.className = 'verdict';
            verdictEl.innerHTML = '';
            scoreTable.innerHTML = '<table>' +
                row('Chain', node.id) +
                row('avg pLDDT', fmt(node.plddt, 1)) +
                row('Interactions', node.degree) + '</table>';
            paeAxis.innerHTML = ''; cliaAxis.innerHTML = '';
            modalRight.style.display = 'none';
            modal.classList.add('visible');
            clearPlots();
            viewer3dEl.innerHTML = '';
            modalViewer = $3Dmol.createViewer(viewer3dEl, { backgroundColor: 'white' });
            if (node.pdb) {
                modalViewer.addModel(node.pdb, 'pdb');
                modalViewer.setStyle({}, {cartoon: {color: 'spectrum'}});
            }
            modalViewer.zoomTo(); modalViewer.render();
            setTimeout(() => { if (modalViewer) modalViewer.resize(); }, 60);
        }

        function closeModal() {
            modal.classList.remove('visible');
            if (modalViewer) { try { modalViewer.clear(); } catch (e) {} modalViewer = null; }
            viewer3dEl.innerHTML = '';
            clearPlots();
        }
        document.getElementById('modal-close').addEventListener('click', closeModal);
        modal.addEventListener('click', (e) => { if (e.target === modal) closeModal(); });
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape' && modal.classList.contains('visible')) closeModal();
        });

        // ===== Edge picking on canvas =====
        function pointToSeg(px, py, x1, y1, x2, y2) {
            const dx = x2 - x1, dy = y2 - y1, lenSq = dx*dx + dy*dy;
            if (lenSq === 0) return Math.hypot(px - x1, py - y1);
            let t = ((px - x1) * dx + (py - y1) * dy) / lenSq;
            t = Math.max(0, Math.min(1, t));
            return Math.hypot(px - (x1 + t*dx), py - (y1 + t*dy));
        }
        function findEdgeAt(mx, my) {
            let best = null, bestD = Infinity;
            edges.forEach(edge => {
                const n1 = nodes.find(n => n.id === edge.chain1);
                const n2 = nodes.find(n => n.id === edge.chain2);
                if (!n1 || !n2) return;
                const d = pointToSeg(mx, my, n1.x+50, n1.y+50, n2.x+50, n2.y+50);
                const thr = Math.max(edge.width, 4) + 4;
                if (d <= thr && d < bestD) { best = edge; bestD = d; }
            });
            return best;
        }
        canvas.addEventListener('click', (e) => {
            const r = canvas.getBoundingClientRect();
            const mx = (e.clientX - r.left) * (canvas.width / r.width);
            const my = (e.clientY - r.top) * (canvas.height / r.height);
            const edge = findEdgeAt(mx, my);
            if (edge) showInteraction(edge);
        });
        canvas.addEventListener('mousemove', (e) => {
            const r = canvas.getBoundingClientRect();
            const mx = (e.clientX - r.left) * (canvas.width / r.width);
            const my = (e.clientY - r.top) * (canvas.height / r.height);
            const edge = findEdgeAt(mx, my);
            if (edge) {
                canvas.classList.add('edge-hover');
                tooltip.innerHTML = HAS_SCORES
                    ? ('<b>' + edge.name + '</b><br>' +
                       'iLIS: ' + edge.iLIS.toFixed(3) + (edge.positive ? ' &#10003;' : '') + '<br>' +
                       'LIS: ' + edge.LIS.toFixed(3) + ' | cLIS: ' + edge.cLIS.toFixed(3) + '<br>' +
                       'ipTM: ' + fmt(edge.ipTM, 3) + '<br><i>Click for PAE / cLIA maps</i>')
                    : ('<b>' + edge.name + '</b><br>' + edge.contacts +
                       ' residue contacts<br><i>Click to view pair (no LIS: scores not provided)</i>');
                tooltip.style.display = 'block';
                tooltip.style.left = (e.pageX + 15) + 'px';
                tooltip.style.top = (e.pageY + 15) + 'px';
            } else {
                canvas.classList.remove('edge-hover');
                if (!draggedNode) tooltip.style.display = 'none';
            }
        });
        canvas.addEventListener('mouseleave', () => {
            canvas.classList.remove('edge-hover');
            tooltip.style.display = 'none';
        });
    </script>
</body>
</html>
"""


def render_html(nodes_js, edges_js, output_file, title, banner, legend_title, info_hint, has_scores):
    html = (HTML_TEMPLATE
            .replace('__TITLE__', title)
            .replace('__BANNER__', banner)
            .replace('__LEGEND_TITLE__', legend_title)
            .replace('__INFO_HINT__', info_hint)
            .replace('__HAS_SCORES__', 'true' if has_scores else 'false')
            .replace('__THRESH__', str(lis_lib.ILIS_THRESHOLD))
            .replace('__NNODES__', str(len(nodes_js)))
            .replace('__NEDGES__', str(len(edges_js)))
            .replace('__NODES_JSON__', json.dumps(nodes_js))
            .replace('__EDGES_JSON__', json.dumps(edges_js)))
    with open(output_file, 'w') as f:
        f.write(html)
    print(f"Saved predicted PPI graph: {output_file}")


def render_scored(model, nodes_js, edges_js, output_file, title):
    banner = (f"<b>PREDICTED MODE</b> &nbsp;|&nbsp; platform: <b>{model['platform']}</b>"
              f" &nbsp;|&nbsp; model: <b>{model['name']} (rank {model['rank']})</b>"
              f" &nbsp;|&nbsp; scoring: LIS / cLIS / iLIS (PAE + contacts)"
              f" &nbsp;|&nbsp; iLIS &ge; {lis_lib.ILIS_THRESHOLD} = likely interaction")
    render_html(nodes_js, edges_js, output_file, title, banner,
                legend_title='Chain pairs (by iLIS)',
                info_hint='Click an edge for PAE / cLIA maps + LIS scores',
                has_scores=True)


def render_unscored(struct_name, nodes_js, edges_js, output_file, title):
    banner = (f"<b>PREDICTED MODE</b> &nbsp;|&nbsp; "
              f"<b style=\"color:#ff8a80\">score files (PAE) not provided</b>"
              f" &nbsp;|&nbsp; structure-only PPI graph (no LIS metrics)"
              f" &nbsp;|&nbsp; source: <b>{struct_name}</b>")
    render_html(nodes_js, edges_js, output_file, title, banner,
                legend_title='Chain pairs (by contacts)',
                info_hint='Click an edge to view the pair complex (no LIS: scores not provided)',
                has_scores=False)


# ============================================================================
# Output text files
# ============================================================================

def save_lis_scores(model, output_file):
    """Write per-chain-pair LIS scores sorted by iLIS."""
    pairs = sorted(model['pairs'], key=lambda p: -p['iLIS'])
    with open(output_file, 'w') as f:
        f.write(f"# Predicted-mode LIS scores | platform={model['platform']} "
                f"| model={model['name']} rank={model['rank']}\n")
        f.write(f"# iLIS >= {lis_lib.ILIS_THRESHOLD} suggests a likely interaction (Kim et al., 2026)\n")
        f.write("chain_i\tchain_j\tiLIS\tLIS\tcLIS\tLIA\tcLIA\tiLIA\tipTM\tpLDDT_i\tpLDDT_j\tcall\n")
        for p in pairs:
            call = 'POSITIVE' if p['iLIS'] >= lis_lib.ILIS_THRESHOLD else 'negative'
            iptm = f"{p['ipTM']:.3f}" if p['ipTM'] is not None else 'NA'
            pi = f"{p['pLDDT_i']:.1f}" if p.get('pLDDT_i') is not None else 'NA'
            pj = f"{p['pLDDT_j']:.1f}" if p.get('pLDDT_j') is not None else 'NA'
            f.write(f"{p['ci']}\t{p['cj']}\t{p['iLIS']:.3f}\t{p['LIS']:.3f}\t{p['cLIS']:.3f}\t"
                    f"{int(p['LIA'])}\t{int(p['cLIA'])}\t{p['iLIA']:.1f}\t{iptm}\t{pi}\t{pj}\t{call}\n")
    print(f"Saved LIS scores: {output_file}")


# ============================================================================
# Modes
# ============================================================================

def run_predicted(path, output_dir, model_req, edge_min_ilis, pae_cutoff, cb_cutoff,
                  contact_cutoff=5.0):
    os.makedirs(output_dir, exist_ok=True)
    try:
        models = lis_lib.analyze_prediction(path, pae_cutoff=pae_cutoff, cb_cutoff=cb_cutoff)
    except ValueError:
        # No PAE/score files -> structure-only PPI graph (no LIS metrics)
        struct_path, fmt = find_structure_file(path)
        if not struct_path:
            raise SystemExit(f"Error: no structure (.cif/.pdb) found at {path}")
        base = os.path.splitext(os.path.basename(struct_path))[0]
        print("No PAE / score files found — generating a structure-only PPI graph "
              "(contacts + pLDDT, no LIS metrics).")
        nodes_js, edges_js = build_structure_only_payload(struct_path, fmt, contact_cutoff)
        html_out = os.path.join(output_dir, f"{base}_predicted_ppi.html")
        render_unscored(os.path.basename(struct_path), nodes_js, edges_js, html_out, base)
        return None

    model = select_model(models, model_req)
    print(f"Platform: {model['platform']} | models found: {len(models)} | "
          f"using model '{model['name']}' rank {model['rank']} "
          f"(mean iLIS {_mean_ilis(model):.3f})")

    nodes_js, edges_js, _G = build_payload(model, edge_min_ilis=edge_min_ilis)
    base = model['name'] or 'prediction'
    html_out = os.path.join(output_dir, f"{base}_predicted_ppi.html")
    txt_out = os.path.join(output_dir, f"{base}_lis_scores.txt")
    render_scored(model, nodes_js, edges_js, html_out, base)
    save_lis_scores(model, txt_out)
    return model


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Predicted-mode PPI graph with LIS/cLIS/iLIS scoring "
                    "(AlphaFold3 / Boltz / Chai-1).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument('input', help='Prediction folder (or single structure file with sibling PAE)')
    parser.add_argument('--predicted', action='store_true',
                        help='Predicted mode (default). One prediction -> one graph.')
    parser.add_argument('--model', default=None,
                        help='Model rank to use (default: best mean iLIS).')
    parser.add_argument('--edge-min-ilis', type=float, default=0.0,
                        help='Only draw edges with iLIS >= this value (default: 0, draw all).')
    parser.add_argument('--pae-cutoff', type=float, default=12.0, help='PAE cutoff (default 12).')
    parser.add_argument('--cb-cutoff', type=float, default=8.0, help='Cb contact cutoff Å (default 8).')
    parser.add_argument('--contact-cutoff', type=float, default=5.0,
                        help='Heavy-atom contact cutoff Å for the structure-only fallback '
                             'when no PAE/score files are provided (default 5).')
    parser.add_argument('--output-dir', default='.', help='Output directory (default: current).')
    args = parser.parse_args()

    try:
        run_predicted(args.input, args.output_dir, args.model, args.edge_min_ilis,
                      args.pae_cutoff, args.cb_cutoff, args.contact_cutoff)
    except ValueError as e:
        raise SystemExit(f"Error: {e}")


if __name__ == '__main__':
    main()
