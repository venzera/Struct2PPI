#!/usr/bin/env python3
"""
Protein-Protein Interaction Graph Generator

Generates interactive protein-protein interaction graphs from PDB/CIF files.
Extracts chain labels from COMPND lines (PDB) or entity info (CIF).
Outputs: HTML visualization, chain info, and residue contacts.

Usage:
    python ppi_graph.py structure.cif
    python ppi_graph.py structure.pdb --cutoff 5.0
"""

import argparse
import os
import warnings
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.spatial.distance import cdist

from Bio.PDB import PDBParser, MMCIFParser, is_aa
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import networkx as nx
import plotly.graph_objects as go

warnings.filterwarnings('ignore', category=PDBConstructionWarning)


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
    """
    Find protein-protein interactions based on distance cutoff.

    Returns:
        interactions: dict of (chain1, chain2) -> number of contacts
        residue_contacts: dict of (chain1, chain2) -> list of (res1, res2) pairs
    """
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


def build_graph(chain_labels, interactions, protein_chain_ids):
    """Build NetworkX graph from interactions."""
    G = nx.Graph()

    for chain_id in protein_chain_ids:
        label = chain_labels.get(chain_id, f"Chain {chain_id}")
        G.add_node(chain_id, label=label)

    for (chain1, chain2), contact_count in interactions.items():
        G.add_edge(chain1, chain2, weight=contact_count)

    return G


def create_plotly_visualization(G, chain_labels, output_file, structure_name):
    """Create interactive Plotly visualization and save as HTML."""
    if len(G.nodes()) == 0:
        print("Warning: No nodes in graph, skipping visualization")
        return

    if len(G.nodes()) == 1:
        pos = {list(G.nodes())[0]: (0, 0)}
    else:
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    edge_colors = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
        '#ff0000', '#00ff00', '#0000ff', '#ffff00', '#ff00ff',
        '#00ffff', '#ff8000', '#8000ff', '#0080ff', '#ff0080',
        '#80ff00', '#00ff80', '#800080', '#008000', '#000080',
        '#808000', '#800000', '#008080', '#c0c0c0', '#400000',
        '#004000', '#000040', '#404000', '#400040', '#004040',
        '#404040', '#200000', '#002000', '#000020', '#202000',
        '#e60000', '#00e600', '#0000e6', '#e6e600', '#e600e6',
        '#00e6e6', '#e68000', '#8000e6', '#0080e6', '#e60080',
        '#b30000', '#00b300', '#0000b3', '#b3b300', '#b300b3',
        '#00b3b3', '#b36600', '#6600b3', '#0066b3', '#b30066',
        '#990000', '#009900', '#000099', '#999900', '#990099',
        '#009999', '#994d00', '#4d0099', '#004d99', '#990066',
        '#660000', '#006600', '#000066', '#666600', '#660066',
        '#006666', '#663300', '#330066', '#003366', '#660033',
        '#cc0000', '#00cc00', '#0000cc', '#cccc00', '#cc00cc',
        '#00cccc', '#cc6600', '#6600cc', '#0066cc', '#cc0066',
    ]

    edge_traces = []
    edge_list = list(G.edges(data=True))

    for i, edge in enumerate(edge_list):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        weight = edge[2].get('weight', 1)

        width = min(1 + np.log1p(weight) * 2, 10)
        color = edge_colors[i % len(edge_colors)]

        label1 = chain_labels.get(edge[0], edge[0])
        label2 = chain_labels.get(edge[1], edge[1])
        legend_name = f"{edge[0]}-{edge[1]}"

        edge_trace = go.Scatter(
            x=[x0, x1, None],
            y=[y0, y1, None],
            mode='lines',
            line=dict(width=width, color=color),
            hoverinfo='text',
            hovertext=f"{edge[0]} ({label1}) - {edge[1]} ({label2}): {weight} contacts",
            name=legend_name,
            showlegend=True,
            legendgroup=legend_name
        )
        edge_traces.append(edge_trace)

    node_x = []
    node_y = []
    node_text = []
    node_hover = []
    node_sizes = []

    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

        label = chain_labels.get(node, f"Chain {node}")
        node_text.append(node)

        degree = G.degree(node)
        total_contacts = sum(G[node][neighbor].get('weight', 0) for neighbor in G.neighbors(node))
        node_hover.append(f"<b>Chain {node}</b><br>{label}<br>Interactions: {degree}<br>Total contacts: {total_contacts}")

        node_sizes.append(20 + min(total_contacts / 10, 30))

    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5'
    ]

    node_colors = [colors[i % len(colors)] for i in range(len(G.nodes()))]

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=node_text,
        textposition='top center',
        hovertext=node_hover,
        marker=dict(
            size=node_sizes,
            color=node_colors,
            line=dict(width=2, color='white')
        ),
        showlegend=False
    )

    fig = go.Figure(
        data=edge_traces + [node_trace],
        layout=go.Layout(
            title=dict(
                text=f'Protein-Protein Interaction Graph: {structure_name}',
                x=0.5,
                font=dict(size=20)
            ),
            showlegend=True,
            hovermode='closest',
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            margin=dict(l=40, r=250, t=60, b=40),
            legend=dict(
                title=dict(text='Interactions', font=dict(size=14)),
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="lightgray",
                borderwidth=1,
                font=dict(size=10),
                itemsizing='constant'
            ),
            annotations=[
                dict(
                    text=f"Nodes: {len(G.nodes())} chains | Edges: {len(G.edges())} interactions",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.5, y=-0.05,
                    font=dict(size=12, color='gray')
                )
            ]
        )
    )

    html_content = fig.to_html(include_plotlyjs=True, full_html=True)

    drag_script = """
<script>
document.addEventListener('DOMContentLoaded', function() {
    var gd = document.querySelector('.plotly-graph-div');
    if (!gd) return;

    var draggedNode = null;
    var nodeTraceIndex = null;
    var edgeData = {};

    // Wait for Plotly to fully render
    setTimeout(function() {
        // Find the node trace (last trace with markers)
        var traces = gd.data;
        for (var i = traces.length - 1; i >= 0; i--) {
            if (traces[i].mode && traces[i].mode.includes('markers')) {
                nodeTraceIndex = i;
                break;
            }
        }

        // Build edge connectivity map
        for (var i = 0; i < traces.length; i++) {
            if (traces[i].mode === 'lines' && traces[i].name) {
                var parts = traces[i].name.split('-');
                if (parts.length === 2) {
                    var chain1 = parts[0], chain2 = parts[1];
                    if (!edgeData[chain1]) edgeData[chain1] = [];
                    if (!edgeData[chain2]) edgeData[chain2] = [];
                    edgeData[chain1].push({traceIdx: i, isStart: true});
                    edgeData[chain2].push({traceIdx: i, isStart: false});
                }
            }
        }

        // Get node labels
        var nodeLabels = traces[nodeTraceIndex].text;

        gd.on('plotly_click', function(data) {
            if (data.points[0].curveNumber === nodeTraceIndex) {
                draggedNode = data.points[0].pointIndex;
                gd.style.cursor = 'grabbing';
            }
        });

        gd.addEventListener('mousemove', function(evt) {
            if (draggedNode === null) return;

            var bb = gd.getBoundingClientRect();
            var xaxis = gd._fullLayout.xaxis;
            var yaxis = gd._fullLayout.yaxis;

            var newX = xaxis.p2d(evt.clientX - bb.left - gd._fullLayout.margin.l);
            var newY = yaxis.p2d(evt.clientY - bb.top - gd._fullLayout.margin.t);

            // Update node position
            var nodeX = traces[nodeTraceIndex].x.slice();
            var nodeY = traces[nodeTraceIndex].y.slice();
            nodeX[draggedNode] = newX;
            nodeY[draggedNode] = newY;

            var update = {'x': [nodeX], 'y': [nodeY]};
            Plotly.restyle(gd, update, [nodeTraceIndex]);

            // Update connected edges
            var nodeLabel = nodeLabels[draggedNode];
            if (edgeData[nodeLabel]) {
                edgeData[nodeLabel].forEach(function(edge) {
                    var edgeX = traces[edge.traceIdx].x.slice();
                    var edgeY = traces[edge.traceIdx].y.slice();
                    if (edge.isStart) {
                        edgeX[0] = newX;
                        edgeY[0] = newY;
                    } else {
                        edgeX[1] = newX;
                        edgeY[1] = newY;
                    }
                    Plotly.restyle(gd, {'x': [edgeX], 'y': [edgeY]}, [edge.traceIdx]);
                });
            }
        });

        gd.addEventListener('mouseup', function() {
            if (draggedNode !== null) {
                draggedNode = null;
                gd.style.cursor = 'default';
            }
        });

        gd.addEventListener('mouseleave', function() {
            if (draggedNode !== null) {
                draggedNode = null;
                gd.style.cursor = 'default';
            }
        });

    }, 500);
});
</script>
"""

    html_content = html_content.replace('</body>', drag_script + '</body>')

    with open(output_file, 'w') as f:
        f.write(html_content)

    print(f"Saved interactive graph: {output_file}")


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
        description='Generate protein-protein interaction graph from PDB/CIF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python ppi_graph.py 8xks.cif
    python ppi_graph.py structure.pdb --cutoff 4.0
    python ppi_graph.py complex.cif --output-dir ./results
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

    html_file = os.path.join(args.output_dir, f"{basename}_ppi_graph.html")
    chain_info_file = os.path.join(args.output_dir, f"{basename}_chain_info.txt")
    contacts_file = os.path.join(args.output_dir, f"{basename}_residue_contacts.txt")

    print("\nGenerating outputs...")
    create_plotly_visualization(G, chain_labels, html_file, basename)
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
