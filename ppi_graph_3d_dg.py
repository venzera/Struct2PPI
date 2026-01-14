#!/usr/bin/env python3
"""
Protein-Protein Interaction Graph Generator with 3D Structure Nodes and Binding Energy

Similar to STRING DB visualization - shows protein structures inside nodes.
Uses py3Dmol for 3D structure visualization embedded in the network graph.
Calculates binding energy (ΔG) using PRODIGY for each interacting chain pair.

Usage:
    python ppi_graph_3d_dg.py structure.pdb
    python ppi_graph_3d_dg.py structure.pdb --cutoff 5.0
"""

import argparse
import os
import subprocess
import re
import warnings
import io
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.spatial.distance import cdist

from Bio.PDB import PDBParser, is_aa, PDBIO, Select
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import networkx as nx

warnings.filterwarnings('ignore', category=PDBConstructionWarning)


class ChainSelect(Select):
    """Select a specific chain for output."""
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == ' ' and is_aa(residue, standard=True)


class ChainPairSelect(Select):
    """Select two specific chains for output (for PRODIGY complex)."""
    def __init__(self, chain_id1, chain_id2):
        self.chain_ids = {chain_id1, chain_id2}

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        return residue.id[0] == ' ' and is_aa(residue, standard=True)


def parse_structure(filepath):
    """Parse PDB file and return structure object."""
    ext = os.path.splitext(filepath)[1].lower()
    structure_id = os.path.splitext(os.path.basename(filepath))[0]

    if ext != '.pdb':
        raise ValueError(f"Unsupported file format: {ext}. Only .pdb files are supported.")

    parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, filepath)


def get_chain_labels(filepath):
    """Extract chain labels from COMPND lines in PDB file."""
    chain_labels = {}
    current_mol = None
    current_chains = []

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('COMPND'):
                content = line[10:].strip()

                if 'MOL_ID:' in content:
                    if current_mol and current_chains:
                        for chain in current_chains:
                            chain_labels[chain] = current_mol
                    current_mol = None
                    current_chains = []

                if 'MOLECULE:' in content:
                    mol_name = content.replace('MOLECULE:', '').strip()
                    mol_name = mol_name.rstrip(';').strip()
                    current_mol = mol_name

                if 'CHAIN:' in content:
                    chains_str = content.replace('CHAIN:', '').strip()
                    chains_str = chains_str.rstrip(';').strip()
                    current_chains = [c.strip() for c in chains_str.split(',')]

        if current_mol and current_chains:
            for chain in current_chains:
                chain_labels[chain] = current_mol

    return chain_labels


def is_protein_chain(chain):
    """Check if chain contains protein residues."""
    aa_count = 0
    total_count = 0

    for residue in chain.get_residues():
        if residue.id[0] == ' ':
            total_count += 1
            if is_aa(residue, standard=True):
                aa_count += 1

    return aa_count > 0 and (aa_count / max(total_count, 1)) > 0.5


def get_heavy_atoms_coords(chain):
    """Get coordinates of all heavy (non-hydrogen) atoms in a chain."""
    coords = []
    atom_info = []

    for residue in chain.get_residues():
        if residue.id[0] != ' ':
            continue
        if not is_aa(residue, standard=True):
            continue

        res_name = residue.resname
        res_id = residue.id[1]

        for atom in residue.get_atoms():
            element = atom.element.strip().upper() if atom.element else ''
            if element == 'H' or element == '':
                if atom.name.startswith('H') or atom.name.startswith('1H') or atom.name.startswith('2H') or atom.name.startswith('3H'):
                    continue

            coords.append(atom.coord)
            atom_info.append((res_name, res_id, atom.name))

    return np.array(coords) if coords else None, atom_info


def find_interactions(structure, distance_cutoff=5.0):
    """Find protein-protein interactions based on distance cutoff."""
    interactions = defaultdict(int)
    residue_contacts = defaultdict(set)

    protein_chains = []
    for model in structure:
        for chain in model:
            if is_protein_chain(chain):
                protein_chains.append(chain)
        break

    print(f"Found {len(protein_chains)} protein chains")

    chain_data = {}
    for chain in protein_chains:
        coords, atom_info = get_heavy_atoms_coords(chain)
        if coords is not None and len(coords) > 0:
            chain_data[chain.id] = (coords, atom_info)

    chain_ids = list(chain_data.keys())
    total_pairs = len(list(combinations(chain_ids, 2)))

    print(f"Calculating interactions between {total_pairs} chain pairs...")

    for i, (chain1_id, chain2_id) in enumerate(combinations(chain_ids, 2)):
        coords1, info1 = chain_data[chain1_id]
        coords2, info2 = chain_data[chain2_id]

        distances = cdist(coords1, coords2)
        contact_mask = distances < distance_cutoff
        contact_indices = np.where(contact_mask)

        if len(contact_indices[0]) > 0:
            pair_key = tuple(sorted([chain1_id, chain2_id]))

            seen_residue_pairs = set()
            for idx1, idx2 in zip(contact_indices[0], contact_indices[1]):
                res1 = (info1[idx1][0], info1[idx1][1])
                res2 = (info2[idx2][0], info2[idx2][1])

                res_pair = (f"{res1[0]}{res1[1]}", f"{res2[0]}{res2[1]}")
                if res_pair not in seen_residue_pairs:
                    seen_residue_pairs.add(res_pair)

                    if pair_key[0] == chain1_id:
                        residue_contacts[pair_key].add((res_pair[0], res_pair[1]))
                    else:
                        residue_contacts[pair_key].add((res_pair[1], res_pair[0]))

            interactions[pair_key] = len(seen_residue_pairs)

        if (i + 1) % 100 == 0 or i + 1 == total_pairs:
            print(f"  Processed {i + 1}/{total_pairs} pairs")

    return dict(interactions), {k: list(v) for k, v in residue_contacts.items()}


def save_chain_pair_pdb(structure, chain1_id, chain2_id, output_path):
    """Save a PDB file containing only two chains for PRODIGY analysis."""
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(output_path, ChainPairSelect(chain1_id, chain2_id))


def run_prodigy(pdb_path, chain1_id, chain2_id):
    """Run PRODIGY on a chain pair and return binding affinity (ΔG) and Kd."""
    try:
        cmd = ['prodigy', pdb_path, '--selection', chain1_id, chain2_id]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        output = result.stdout + result.stderr

        # Parse binding affinity (ΔG in kcal/mol)
        dg_match = re.search(r'Predicted binding affinity.*?:\s*([-\d.]+)', output)
        kd_match = re.search(r'Predicted dissociation constant.*?:\s*([\d.eE+-]+)', output)

        dg = float(dg_match.group(1)) if dg_match else None
        kd = kd_match.group(1) if kd_match else None

        return dg, kd, output

    except subprocess.TimeoutExpired:
        print(f"  Warning: PRODIGY timeout for {chain1_id}-{chain2_id}")
        return None, None, "Timeout"
    except FileNotFoundError:
        print("  Error: PRODIGY not found. Please install with: pip install prodigy-prot")
        return None, None, "PRODIGY not installed"
    except Exception as e:
        print(f"  Warning: PRODIGY error for {chain1_id}-{chain2_id}: {e}")
        return None, None, str(e)


def calculate_binding_energies(structure, interactions, complexes_dir):
    """Calculate binding energies for all interacting chain pairs using PRODIGY."""
    binding_data = {}

    os.makedirs(complexes_dir, exist_ok=True)

    total_pairs = len(interactions)
    print(f"\nCalculating binding energies for {total_pairs} interacting pairs...")

    for i, (chain1, chain2) in enumerate(interactions.keys()):
        print(f"  [{i+1}/{total_pairs}] Processing {chain1}-{chain2}...", end=" ")

        # Save chain pair PDB
        pdb_filename = f"complex_{chain1}_{chain2}.pdb"
        pdb_path = os.path.join(complexes_dir, pdb_filename)
        save_chain_pair_pdb(structure, chain1, chain2, pdb_path)

        # Run PRODIGY
        dg, kd, raw_output = run_prodigy(pdb_path, chain1, chain2)

        binding_data[(chain1, chain2)] = {
            'dG': dg,
            'Kd': kd,
            'pdb_file': pdb_filename,
            'raw_output': raw_output
        }

        if dg is not None:
            print(f"ΔG = {dg:.1f} kcal/mol")
        else:
            print("Failed")

    return binding_data


def get_chain_pdb_string(structure, chain_id):
    """Extract a single chain as PDB string."""
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    output = io.StringIO()
    pdb_io.save(output, ChainSelect(chain_id))
    return output.getvalue()


def build_graph(chain_labels, interactions, protein_chain_ids, binding_data):
    """Build NetworkX graph from interactions with binding energy data."""
    G = nx.Graph()

    for chain_id in protein_chain_ids:
        label = chain_labels.get(chain_id, f"Chain {chain_id}")
        G.add_node(chain_id, label=label)

    for (chain1, chain2), contact_count in interactions.items():
        pair_key = tuple(sorted([chain1, chain2]))
        bd = binding_data.get(pair_key, {})
        dg = bd.get('dG')
        kd = bd.get('Kd')
        G.add_edge(chain1, chain2, weight=contact_count, dG=dg, Kd=kd)

    return G


def create_3d_visualization(G, chain_labels, structure, interactions, binding_data, output_file, structure_name):
    """Create interactive visualization with 3D structure nodes using py3Dmol."""

    if len(G.nodes()) == 0:
        print("Warning: No nodes in graph, skipping visualization")
        return

    if len(G.nodes()) == 1:
        pos = {list(G.nodes())[0]: (0, 0)}
    else:
        pos = nx.spring_layout(G, k=3, iterations=100, seed=42)

    x_coords = [p[0] for p in pos.values()]
    y_coords = [p[1] for p in pos.values()]
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    x_range = x_max - x_min if x_max != x_min else 1
    y_range = y_max - y_min if y_max != y_min else 1

    normalized_pos = {}
    for node, (x, y) in pos.items():
        nx_pos = 100 + (x - x_min) / x_range * 1200
        ny_pos = 100 + (y - y_min) / y_range * 800
        normalized_pos[node] = (nx_pos, ny_pos)

    edge_colors = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
        '#ff0000', '#00ff00', '#0000ff', '#ffff00', '#ff00ff',
        '#00ffff', '#ff8000', '#8000ff', '#0080ff', '#ff0080',
    ]

    chain_pdb_data = {}
    print("Extracting chain structures...")
    for node in G.nodes():
        pdb_str = get_chain_pdb_string(structure, node)
        chain_pdb_data[node] = pdb_str.replace('\n', '\\n').replace("'", "\\'")

    edges_js = []
    edge_list = list(G.edges(data=True))
    for i, (n1, n2, data) in enumerate(edge_list):
        x1, y1 = normalized_pos[n1]
        x2, y2 = normalized_pos[n2]
        weight = data.get('weight', 1)
        dg = data.get('dG')
        kd = data.get('Kd')
        color = edge_colors[i % len(edge_colors)]
        width = min(2 + np.log1p(weight), 8)
        label1 = chain_labels.get(n1, n1)
        label2 = chain_labels.get(n2, n2)

        # Build edge label with binding energy
        dg_str = f"{dg:.1f}" if dg is not None else "N/A"
        edge_label = f"{n1}-{n2}: {weight} contacts, ΔG={dg_str} kcal/mol"

        edges_js.append({
            'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2,
            'color': color, 'width': width,
            'name': f"{n1}-{n2}",
            'label': f"{n1} ({label1}) - {n2} ({label2}): {weight} contacts",
            'dG': dg_str,
            'Kd': kd if kd else "N/A"
        })

    nodes_js = []
    for node in G.nodes():
        x, y = normalized_pos[node]
        label = chain_labels.get(node, f"Chain {node}")
        degree = G.degree(node)
        total_contacts = sum(G[node][neighbor].get('weight', 0) for neighbor in G.neighbors(node))
        nodes_js.append({
            'id': node,
            'x': x, 'y': y,
            'label': label,
            'degree': degree,
            'contacts': total_contacts,
            'pdb': chain_pdb_data.get(node, '')
        })

    html_template = """<!DOCTYPE html>
<html>
<head>
    <title>PPI Graph with 3D Structures: {structure_name}</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        h1 {{
            text-align: center;
            color: #333;
        }}
        #container {{
            position: relative;
            width: 1400px;
            height: 1000px;
            margin: 0 auto;
            background: white;
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
        }}
        #graph-canvas {{
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            z-index: 1;
        }}
        .node-container {{
            position: absolute;
            width: 100px;
            height: 100px;
            border-radius: 50%;
            border: 3px solid #333;
            background: white;
            cursor: move;
            z-index: 10;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
        }}
        .node-container:hover {{
            border-color: #0066cc;
            box-shadow: 0 4px 20px rgba(0,102,204,0.4);
        }}
        .node-viewer {{
            width: 100%;
            height: 100%;
            border-radius: 50%;
        }}
        .node-label {{
            position: absolute;
            bottom: -25px;
            left: 50%;
            transform: translateX(-50%);
            font-size: 12px;
            font-weight: bold;
            white-space: nowrap;
            background: rgba(255,255,255,0.9);
            padding: 2px 6px;
            border-radius: 3px;
        }}
        #legend {{
            position: absolute;
            top: 10px;
            right: 10px;
            background: rgba(255,255,255,0.95);
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 10px;
            max-height: 900px;
            overflow-y: auto;
            z-index: 100;
            font-size: 11px;
            max-width: 280px;
            transition: all 0.3s ease;
        }}
        #legend.collapsed {{
            max-width: 40px;
            max-height: 40px;
            overflow: hidden;
            padding: 5px;
        }}
        #legend.collapsed #legend-items,
        #legend.collapsed h3 span {{
            display: none;
        }}
        #legend-toggle {{
            position: absolute;
            top: 5px;
            right: 5px;
            width: 24px;
            height: 24px;
            border: none;
            background: #0066cc;
            color: white;
            border-radius: 3px;
            cursor: pointer;
            font-size: 14px;
            line-height: 24px;
            text-align: center;
        }}
        #legend-toggle:hover {{
            background: #0055aa;
        }}
        #legend h3 {{
            margin: 0 30px 10px 0;
            font-size: 14px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin: 3px 0;
            cursor: pointer;
            padding: 2px;
        }}
        .legend-item:hover {{
            background: #f0f0f0;
        }}
        .legend-color {{
            width: 20px;
            height: 3px;
            margin-right: 8px;
            flex-shrink: 0;
        }}
        .legend-text {{
            display: flex;
            flex-direction: column;
        }}
        .legend-name {{
            font-weight: bold;
        }}
        .legend-dg {{
            color: #0066cc;
            font-size: 10px;
        }}
        #info {{
            text-align: center;
            margin-top: 10px;
            color: #666;
        }}
        #tooltip {{
            position: absolute;
            background: rgba(0,0,0,0.8);
            color: white;
            padding: 8px 12px;
            border-radius: 4px;
            font-size: 12px;
            pointer-events: none;
            z-index: 1000;
            display: none;
            max-width: 300px;
        }}
        .edge-label {{
            position: absolute;
            font-size: 10px;
            font-weight: bold;
            color: #0066cc;
            background: rgba(255,255,255,0.9);
            padding: 1px 4px;
            border-radius: 3px;
            z-index: 5;
            pointer-events: none;
        }}
    </style>
</head>
<body>
    <h1>Protein-Protein Interaction Graph: {structure_name}</h1>
    <div id="container">
        <canvas id="graph-canvas"></canvas>
        <div id="legend">
            <button id="legend-toggle" title="Toggle legend">≡</button>
            <h3><span>Interactions (by binding strength)</span></h3>
            <div id="legend-items"></div>
        </div>
    </div>
    <div id="info">
        Nodes: {num_nodes} chains | Edges: {num_edges} interactions | Drag nodes to rearrange
    </div>
    <div id="tooltip"></div>

    <script>
        const nodes = {nodes_json};
        const edges = {edges_json};

        // Sort edges by binding strength (more negative = stronger)
        const sortedEdges = [...edges].sort((a, b) => {{
            const dgA = a.dG === 'N/A' ? 0 : parseFloat(a.dG);
            const dgB = b.dG === 'N/A' ? 0 : parseFloat(b.dG);
            return dgA - dgB;  // More negative first
        }});

        const container = document.getElementById('container');
        const canvas = document.getElementById('graph-canvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');

        canvas.width = 1400;
        canvas.height = 1000;

        const nodeElements = {{}};
        const viewers = {{}};
        const edgeLabelElements = {{}};
        let draggedNode = null;
        let dragOffset = {{x: 0, y: 0}};

        function drawEdges() {{
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            edges.forEach(edge => {{
                const n1 = nodes.find(n => n.id === edge.name.split('-')[0]);
                const n2 = nodes.find(n => n.id === edge.name.split('-')[1]);

                if (n1 && n2) {{
                    ctx.beginPath();
                    ctx.moveTo(n1.x + 50, n1.y + 50);
                    ctx.lineTo(n2.x + 50, n2.y + 50);
                    ctx.strokeStyle = edge.color;
                    ctx.lineWidth = edge.width;
                    ctx.stroke();

                    // Update edge label position
                    const labelEl = edgeLabelElements[edge.name];
                    if (labelEl) {{
                        const midX = (n1.x + n2.x) / 2 + 50;
                        const midY = (n1.y + n2.y) / 2 + 50;
                        labelEl.style.left = midX + 'px';
                        labelEl.style.top = midY + 'px';
                    }}
                }}
            }});
        }}

        function createEdgeLabels() {{
            edges.forEach(edge => {{
                const n1 = nodes.find(n => n.id === edge.name.split('-')[0]);
                const n2 = nodes.find(n => n.id === edge.name.split('-')[1]);

                if (n1 && n2 && edge.dG !== 'N/A') {{
                    const labelEl = document.createElement('div');
                    labelEl.className = 'edge-label';
                    labelEl.textContent = edge.dG;
                    labelEl.title = 'ΔG (kcal/mol)';

                    const midX = (n1.x + n2.x) / 2 + 50;
                    const midY = (n1.y + n2.y) / 2 + 50;
                    labelEl.style.left = midX + 'px';
                    labelEl.style.top = midY + 'px';

                    container.appendChild(labelEl);
                    edgeLabelElements[edge.name] = labelEl;
                }}
            }});
        }}

        function createNodes() {{
            nodes.forEach(node => {{
                const div = document.createElement('div');
                div.className = 'node-container';
                div.id = 'node-' + node.id;
                div.style.left = (node.x) + 'px';
                div.style.top = (node.y) + 'px';

                const viewer = document.createElement('div');
                viewer.className = 'node-viewer';
                viewer.id = 'viewer-' + node.id;
                div.appendChild(viewer);

                const label = document.createElement('div');
                label.className = 'node-label';
                label.textContent = node.id;
                div.appendChild(label);

                div.addEventListener('mousedown', (e) => {{
                    if (e.target === div || e.target === viewer) {{
                        draggedNode = node;
                        const rect = div.getBoundingClientRect();
                        const containerRect = container.getBoundingClientRect();
                        dragOffset.x = e.clientX - rect.left;
                        dragOffset.y = e.clientY - rect.top;
                        div.style.zIndex = 100;
                    }}
                }});

                div.addEventListener('mouseenter', () => {{
                    tooltip.innerHTML = '<b>Chain ' + node.id + '</b><br>' +
                        node.label + '<br>' +
                        'Interactions: ' + node.degree + '<br>' +
                        'Total contacts: ' + node.contacts;
                    tooltip.style.display = 'block';
                }});

                div.addEventListener('mousemove', (e) => {{
                    tooltip.style.left = (e.pageX + 15) + 'px';
                    tooltip.style.top = (e.pageY + 15) + 'px';
                }});

                div.addEventListener('mouseleave', () => {{
                    tooltip.style.display = 'none';
                }});

                container.appendChild(div);
                nodeElements[node.id] = div;
            }});
        }}

        function initViewers() {{
            nodes.forEach(node => {{
                const element = document.getElementById('viewer-' + node.id);
                if (element && node.pdb) {{
                    const viewer = $3Dmol.createViewer(element, {{
                        backgroundColor: 'white'
                    }});

                    const pdbData = node.pdb.replace(/\\\\n/g, '\\n');
                    viewer.addModel(pdbData, 'pdb');
                    viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                    viewer.zoomTo();
                    viewer.zoom(0.8);
                    viewer.render();
                    viewers[node.id] = viewer;
                }}
            }});
        }}

        function createLegend() {{
            const legendItems = document.getElementById('legend-items');
            sortedEdges.forEach(edge => {{
                const item = document.createElement('div');
                item.className = 'legend-item';
                item.innerHTML = '<div class="legend-color" style="background:' + edge.color + '"></div>' +
                    '<div class="legend-text">' +
                    '<span class="legend-name">' + edge.name + '</span>' +
                    '<span class="legend-dg">ΔG: ' + edge.dG + ' kcal/mol</span>' +
                    '</div>';
                item.title = edge.label + '\\nKd: ' + edge.Kd;
                legendItems.appendChild(item);
            }});
        }}

        document.addEventListener('mousemove', (e) => {{
            if (draggedNode) {{
                const containerRect = container.getBoundingClientRect();
                let newX = e.clientX - containerRect.left - dragOffset.x;
                let newY = e.clientY - containerRect.top - dragOffset.y;

                newX = Math.max(0, Math.min(newX, container.offsetWidth - 100));
                newY = Math.max(0, Math.min(newY, container.offsetHeight - 100));

                draggedNode.x = newX;
                draggedNode.y = newY;

                const div = nodeElements[draggedNode.id];
                div.style.left = newX + 'px';
                div.style.top = newY + 'px';

                drawEdges();
            }}
        }});

        document.addEventListener('mouseup', () => {{
            if (draggedNode) {{
                nodeElements[draggedNode.id].style.zIndex = 10;
                draggedNode = null;
            }}
        }});

        createNodes();
        createEdgeLabels();
        drawEdges();
        createLegend();

        setTimeout(initViewers, 100);

        // Legend toggle functionality
        document.getElementById('legend-toggle').addEventListener('click', function() {{
            const legend = document.getElementById('legend');
            legend.classList.toggle('collapsed');
            this.textContent = legend.classList.contains('collapsed') ? '+' : '≡';
        }});
    </script>
</body>
</html>
"""

    import json

    html_content = html_template.format(
        structure_name=structure_name,
        num_nodes=len(G.nodes()),
        num_edges=len(G.edges()),
        nodes_json=json.dumps(nodes_js),
        edges_json=json.dumps(edges_js)
    )

    with open(output_file, 'w') as f:
        f.write(html_content)

    print(f"Saved 3D interactive graph: {output_file}")


def save_chain_info(chain_labels, protein_chain_ids, output_file):
    """Save chain information to text file."""
    with open(output_file, 'w') as f:
        f.write("Chain ID\tDescription\n")
        f.write("-" * 60 + "\n")

        for chain_id in sorted(protein_chain_ids):
            label = chain_labels.get(chain_id, "Unknown")
            f.write(f"{chain_id}\t{label}\n")

    print(f"Saved chain info: {output_file}")


def save_residue_contacts(residue_contacts, chain_labels, binding_data, output_file):
    """Save residue contacts with binding energy data to text file."""
    with open(output_file, 'w') as f:
        f.write("Chain Pair\tChain 1 Description\tChain 2 Description\tContacts\tΔG (kcal/mol)\tKd (M)\tInteracting Residues\n")
        f.write("=" * 140 + "\n")

        for (chain1, chain2), contacts in sorted(residue_contacts.items()):
            label1 = chain_labels.get(chain1, "Unknown")
            label2 = chain_labels.get(chain2, "Unknown")

            # Get binding data
            pair_key = tuple(sorted([chain1, chain2]))
            bd = binding_data.get(pair_key, {})
            dg = bd.get('dG')
            kd = bd.get('Kd')
            dg_str = f"{dg:.2f}" if dg is not None else "N/A"
            kd_str = kd if kd else "N/A"

            contacts_sorted = sorted(contacts, key=lambda x: (x[0], x[1]))

            f.write(f"\n{chain1}-{chain2}\t{label1}\t{label2}\t{len(contacts)} contacts\tΔG: {dg_str}\tKd: {kd_str}\n")
            f.write("-" * 100 + "\n")

            for res1, res2 in contacts_sorted:
                f.write(f"  {chain1}:{res1} <-> {chain2}:{res2}\n")

    print(f"Saved residue contacts: {output_file}")


def save_binding_strength_ranking(binding_data, chain_labels, output_file):
    """Save chain pairs sorted by binding strength (more negative ΔG = stronger binding)."""
    # Filter pairs with valid binding data
    valid_pairs = [(k, v) for k, v in binding_data.items() if v.get('dG') is not None]

    # Sort by ΔG (ascending - more negative = stronger)
    sorted_pairs = sorted(valid_pairs, key=lambda x: x[1]['dG'])

    with open(output_file, 'w') as f:
        f.write("# Chain pairs sorted by binding strength (strongest first)\n")
        f.write("# More negative ΔG = stronger binding\n")
        f.write("=" * 100 + "\n\n")
        f.write(f"{'Rank':<6}{'Chain Pair':<15}{'ΔG (kcal/mol)':<18}{'Kd (M)':<15}{'Chain 1':<30}{'Chain 2':<30}\n")
        f.write("-" * 100 + "\n")

        for rank, ((chain1, chain2), data) in enumerate(sorted_pairs, 1):
            dg = data['dG']
            kd = data.get('Kd', 'N/A')
            label1 = chain_labels.get(chain1, "Unknown")
            label2 = chain_labels.get(chain2, "Unknown")

            # Truncate long labels
            label1_short = label1[:28] + ".." if len(label1) > 30 else label1
            label2_short = label2[:28] + ".." if len(label2) > 30 else label2

            f.write(f"{rank:<6}{chain1}-{chain2:<13}{dg:<18.2f}{kd:<15}{label1_short:<30}{label2_short:<30}\n")

        f.write("\n" + "=" * 100 + "\n")
        f.write(f"\nTotal pairs with binding data: {len(sorted_pairs)}\n")

        if sorted_pairs:
            strongest = sorted_pairs[0]
            weakest = sorted_pairs[-1]
            f.write(f"Strongest interaction: {strongest[0][0]}-{strongest[0][1]} (ΔG = {strongest[1]['dG']:.2f} kcal/mol)\n")
            f.write(f"Weakest interaction: {weakest[0][0]}-{weakest[0][1]} (ΔG = {weakest[1]['dG']:.2f} kcal/mol)\n")

    print(f"Saved binding strength ranking: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate PPI graph with 3D structure nodes and binding energy (PRODIGY)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python ppi_graph_3d_dg.py structure.pdb
    python ppi_graph_3d_dg.py structure.pdb --cutoff 4.0
        """
    )

    parser.add_argument('input_file', help='Input PDB file')
    parser.add_argument('--cutoff', type=float, default=5.0,
                        help='Distance cutoff for interactions in Angstroms (default: 5.0)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--skip-prodigy', action='store_true',
                        help='Skip PRODIGY calculation (for testing)')

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Error: File not found: {args.input_file}")
        return 1

    ext = os.path.splitext(args.input_file)[1].lower()
    if ext != '.pdb':
        print(f"Error: Only PDB files are supported. Got: {ext}")
        return 1

    os.makedirs(args.output_dir, exist_ok=True)

    basename = os.path.splitext(os.path.basename(args.input_file))[0]

    # Create complexes directory for PRODIGY PDB files
    complexes_dir = os.path.join(args.output_dir, f"{basename}_complexes")

    print(f"Processing: {args.input_file}")
    print(f"Distance cutoff: {args.cutoff} A")
    print()

    print("Parsing structure...")
    structure = parse_structure(args.input_file)

    print("Extracting chain labels...")
    chain_labels = get_chain_labels(args.input_file)

    if chain_labels:
        print(f"Found labels for {len(chain_labels)} chains")
    else:
        print("No chain labels found in file")

    print("\nFinding protein-protein interactions...")
    interactions, residue_contacts = find_interactions(structure, args.cutoff)

    protein_chain_ids = set()
    for model in structure:
        for chain in model:
            if is_protein_chain(chain):
                protein_chain_ids.add(chain.id)
        break

    protein_chain_labels = {k: v for k, v in chain_labels.items() if k in protein_chain_ids}
    chain_labels = protein_chain_labels

    print(f"\nFound {len(interactions)} interacting chain pairs")

    # Calculate binding energies using PRODIGY
    binding_data = {}
    if not args.skip_prodigy and interactions:
        binding_data = calculate_binding_energies(structure, interactions, complexes_dir)
    elif args.skip_prodigy:
        print("\nSkipping PRODIGY calculations (--skip-prodigy flag)")

    print("\nBuilding interaction graph...")
    G = build_graph(chain_labels, interactions, protein_chain_ids, binding_data)

    # Output files
    html_file = os.path.join(args.output_dir, f"{basename}_ppi_3d_dg.html")
    chain_info_file = os.path.join(args.output_dir, f"{basename}_chain_info.txt")
    contacts_file = os.path.join(args.output_dir, f"{basename}_residue_contacts.txt")
    binding_file = os.path.join(args.output_dir, f"{basename}_binding_strength.txt")

    print("\nGenerating outputs...")
    create_3d_visualization(G, chain_labels, structure, interactions, binding_data, html_file, basename)
    save_chain_info(chain_labels, protein_chain_ids, chain_info_file)
    save_residue_contacts(residue_contacts, chain_labels, binding_data, contacts_file)

    if binding_data:
        save_binding_strength_ranking(binding_data, chain_labels, binding_file)

    print("\nDone!")
    print(f"\nSummary:")
    print(f"  Protein chains: {len(protein_chain_ids)}")
    print(f"  Interacting pairs: {len(interactions)}")
    print(f"  Total contacts: {sum(interactions.values())}")
    print(f"  Pairs with binding data: {len([v for v in binding_data.values() if v.get('dG') is not None])}")
    print(f"\nOutput files:")
    print(f"  HTML visualization: {html_file}")
    print(f"  Chain info: {chain_info_file}")
    print(f"  Residue contacts: {contacts_file}")
    if binding_data:
        print(f"  Binding strength: {binding_file}")
        print(f"  Complex PDBs: {complexes_dir}/")

    return 0


if __name__ == '__main__':
    exit(main())
