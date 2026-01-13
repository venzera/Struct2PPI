#!/usr/bin/env python3
"""
Protein-Protein Interaction Graph Generator with 3D Structure Nodes

Similar to STRING DB visualization - shows protein structures inside nodes.
Uses py3Dmol for 3D structure visualization embedded in the network graph.

Usage:
    python ppi_graph_3d.py structure.cif
    python ppi_graph_3d.py structure.pdb --cutoff 5.0
"""

import argparse
import os
import warnings
import base64
import io
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.spatial.distance import cdist

from Bio.PDB import PDBParser, MMCIFParser, is_aa, PDBIO, Select
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

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


def parse_structure(filepath):
    """Parse PDB or CIF file and return structure object."""
    ext = os.path.splitext(filepath)[1].lower()
    structure_id = os.path.splitext(os.path.basename(filepath))[0]

    if ext == '.cif':
        parser = MMCIFParser(QUIET=True)
    elif ext == '.pdb':
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .pdb or .cif")

    return parser.get_structure(structure_id, filepath)


def get_chain_labels_from_pdb(filepath):
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


def get_chain_labels_from_cif(filepath):
    """Extract chain labels from CIF file entity information (polymers only)."""
    chain_labels = {}

    try:
        mmcif_dict = MMCIF2Dict(filepath)

        entity_ids = mmcif_dict.get('_entity.id', [])
        entity_descriptions = mmcif_dict.get('_entity.pdbx_description', [])
        entity_types = mmcif_dict.get('_entity.type', [])

        if isinstance(entity_ids, str):
            entity_ids = [entity_ids]
        if isinstance(entity_descriptions, str):
            entity_descriptions = [entity_descriptions]
        if isinstance(entity_types, str):
            entity_types = [entity_types]

        polymer_entities = set()
        entity_to_desc = {}
        for i, eid in enumerate(entity_ids):
            etype = entity_types[i] if i < len(entity_types) else 'polymer'
            if etype == 'polymer':
                polymer_entities.add(eid)
                if i < len(entity_descriptions):
                    desc = entity_descriptions[i]
                    if desc and desc != '?':
                        entity_to_desc[eid] = desc

        asym_ids = mmcif_dict.get('_struct_asym.id', [])
        entity_refs = mmcif_dict.get('_struct_asym.entity_id', [])

        if isinstance(asym_ids, str):
            asym_ids = [asym_ids]
        if isinstance(entity_refs, str):
            entity_refs = [entity_refs]

        polymer_label_ids = set()
        for i, asym_id in enumerate(asym_ids):
            if i < len(entity_refs):
                entity_id = entity_refs[i]
                if entity_id in polymer_entities:
                    polymer_label_ids.add(asym_id)
                    if entity_id in entity_to_desc:
                        chain_labels[asym_id] = entity_to_desc[entity_id]

        auth_asym_ids = mmcif_dict.get('_atom_site.auth_asym_id', [])
        label_asym_ids = mmcif_dict.get('_atom_site.label_asym_id', [])

        if isinstance(auth_asym_ids, str):
            auth_asym_ids = [auth_asym_ids]
        if isinstance(label_asym_ids, str):
            label_asym_ids = [label_asym_ids]

        label_to_auth = {}
        for i, label_id in enumerate(label_asym_ids):
            if label_id in polymer_label_ids:
                if i < len(auth_asym_ids):
                    auth_id = auth_asym_ids[i]
                    if label_id not in label_to_auth:
                        label_to_auth[label_id] = auth_id

        updated_labels = {}
        for label_id, desc in chain_labels.items():
            if label_id in label_to_auth:
                auth_id = label_to_auth[label_id]
                if auth_id not in updated_labels:
                    updated_labels[auth_id] = desc
            else:
                updated_labels[label_id] = desc

        chain_labels = updated_labels

    except Exception as e:
        print(f"Warning: Could not parse CIF entity info: {e}")

    return chain_labels


def get_chain_labels(filepath):
    """Get chain labels from PDB or CIF file."""
    ext = os.path.splitext(filepath)[1].lower()

    if ext == '.cif':
        return get_chain_labels_from_cif(filepath)
    elif ext == '.pdb':
        return get_chain_labels_from_pdb(filepath)
    else:
        return {}


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


def get_chain_pdb_string(structure, chain_id):
    """Extract a single chain as PDB string."""
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    output = io.StringIO()
    pdb_io.save(output, ChainSelect(chain_id))
    return output.getvalue()


def build_graph(chain_labels, interactions, protein_chain_ids):
    """Build NetworkX graph from interactions."""
    G = nx.Graph()

    for chain_id in protein_chain_ids:
        label = chain_labels.get(chain_id, f"Chain {chain_id}")
        G.add_node(chain_id, label=label)

    for (chain1, chain2), contact_count in interactions.items():
        G.add_edge(chain1, chain2, weight=contact_count)

    return G


def create_3d_visualization(G, chain_labels, structure, interactions, output_file, structure_name):
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
        nx_pos = 100 + (x - x_min) / x_range * 800
        ny_pos = 100 + (y - y_min) / y_range * 600
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
        color = edge_colors[i % len(edge_colors)]
        width = min(2 + np.log1p(weight), 8)
        label1 = chain_labels.get(n1, n1)
        label2 = chain_labels.get(n2, n2)
        edges_js.append({
            'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2,
            'color': color, 'width': width,
            'name': f"{n1}-{n2}",
            'label': f"{n1} ({label1}) - {n2} ({label2}): {weight} contacts"
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
            width: 1000px;
            height: 800px;
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
            max-height: 700px;
            overflow-y: auto;
            z-index: 100;
            font-size: 11px;
        }}
        #legend h3 {{
            margin: 0 0 10px 0;
            font-size: 14px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin: 3px 0;
            cursor: pointer;
        }}
        .legend-item:hover {{
            background: #f0f0f0;
        }}
        .legend-color {{
            width: 20px;
            height: 3px;
            margin-right: 8px;
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
        }}
    </style>
</head>
<body>
    <h1>Protein-Protein Interaction Graph: {structure_name}</h1>
    <div id="container">
        <canvas id="graph-canvas"></canvas>
        <div id="legend">
            <h3>Interactions</h3>
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

        const container = document.getElementById('container');
        const canvas = document.getElementById('graph-canvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');

        canvas.width = 1000;
        canvas.height = 800;

        const nodeElements = {{}};
        const viewers = {{}};
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
            edges.forEach(edge => {{
                const item = document.createElement('div');
                item.className = 'legend-item';
                item.innerHTML = '<div class="legend-color" style="background:' + edge.color + '"></div>' +
                    '<span>' + edge.name + '</span>';
                item.title = edge.label;
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
        drawEdges();
        createLegend();

        setTimeout(initViewers, 100);
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


def save_residue_contacts(residue_contacts, chain_labels, output_file):
    """Save residue contacts to text file."""
    with open(output_file, 'w') as f:
        f.write("Chain Pair\tChain 1 Description\tChain 2 Description\tContacts\tInteracting Residues\n")
        f.write("=" * 120 + "\n")

        for (chain1, chain2), contacts in sorted(residue_contacts.items()):
            label1 = chain_labels.get(chain1, "Unknown")
            label2 = chain_labels.get(chain2, "Unknown")

            contacts_sorted = sorted(contacts, key=lambda x: (x[0], x[1]))

            f.write(f"\n{chain1}-{chain2}\t{label1}\t{label2}\t{len(contacts)} contacts\n")
            f.write("-" * 80 + "\n")

            for res1, res2 in contacts_sorted:
                f.write(f"  {chain1}:{res1} <-> {chain2}:{res2}\n")

    print(f"Saved residue contacts: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate PPI graph with 3D structure nodes (STRING DB style)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python ppi_graph_3d.py 8xks.cif
    python ppi_graph_3d.py structure.pdb --cutoff 4.0
        """
    )

    parser.add_argument('input_file', help='Input PDB or CIF file')
    parser.add_argument('--cutoff', type=float, default=5.0,
                        help='Distance cutoff for interactions in Angstroms (default: 5.0)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory (default: current directory)')

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Error: File not found: {args.input_file}")
        return 1

    os.makedirs(args.output_dir, exist_ok=True)

    basename = os.path.splitext(os.path.basename(args.input_file))[0]

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

    print("\nBuilding interaction graph...")
    G = build_graph(chain_labels, interactions, protein_chain_ids)

    html_file = os.path.join(args.output_dir, f"{basename}_ppi_3d.html")
    chain_info_file = os.path.join(args.output_dir, f"{basename}_chain_info.txt")
    contacts_file = os.path.join(args.output_dir, f"{basename}_residue_contacts.txt")

    print("\nGenerating outputs...")
    create_3d_visualization(G, chain_labels, structure, interactions, html_file, basename)
    save_chain_info(chain_labels, protein_chain_ids, chain_info_file)
    save_residue_contacts(residue_contacts, chain_labels, contacts_file)

    print("\nDone!")
    print(f"\nSummary:")
    print(f"  Protein chains: {len(protein_chain_ids)}")
    print(f"  Interacting pairs: {len(interactions)}")
    print(f"  Total contacts: {sum(interactions.values())}")

    return 0


if __name__ == '__main__':
    exit(main())
