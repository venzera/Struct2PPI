#!/usr/bin/env python3
"""
Protein-Protein Interaction Graph Generator with 3D Structure Nodes and Binding Energy

Similar to STRING DB visualization - shows protein structures inside nodes.
Uses py3Dmol for 3D structure visualization embedded in the network graph.
Calculates binding energy (ΔG) using PRODIGY for each interacting chain pair.

Usage:
    python ppi_graph_3d_dg.py structure.pdb
    python ppi_graph_3d_dg.py structure.cif
    python ppi_graph_3d_dg.py structure.pdb --cutoff 5.0
"""

import argparse
import os
import subprocess
import re
import shutil
import warnings
import io
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.spatial.distance import cdist

from Bio.PDB import PDBParser, is_aa, PDBIO, Select
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import gemmi
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
    structure_id = os.path.splitext(os.path.basename(filepath))[0]
    parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, filepath)


def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert CIF to PDB using gemmi (AlphaFold3 & general CIF support)."""
    structure = gemmi.read_structure(cif_path)
    structure.setup_entities()
    structure.write_pdb(pdb_path)
    return pdb_path


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
    """Extract chain labels from CIF using gemmi.

    Handles missing pdbx_description gracefully (e.g. AlphaFold3 output)
    by returning {} for unnamed chains — the caller falls back to the
    chain id.
    """
    chain_labels = {}
    try:
        doc = gemmi.cif.read(filepath)
        block = doc.sole_block() if len(doc) == 1 else doc[0]

        entity_to_desc = {}
        polymer_entities = set()
        entity_loop = block.find('_entity.', ['id', 'pdbx_description', 'type'])
        for row in entity_loop:
            eid = row[0]
            desc = row[1] if len(row) > 1 else ''
            etype = row[2] if len(row) > 2 else ''
            if etype == 'polymer':
                polymer_entities.add(eid)
                if desc and desc not in ('?', '.'):
                    entity_to_desc[eid] = desc.strip('"').strip("'")

        label_to_entity = {}
        asym_loop = block.find('_struct_asym.', ['id', 'entity_id'])
        for row in asym_loop:
            label_to_entity[row[0]] = row[1]

        label_to_auth = {}
        atom_loop = block.find('_atom_site.', ['label_asym_id', 'auth_asym_id'])
        for row in atom_loop:
            label_id = row[0]
            auth_id = row[1]
            if label_id not in label_to_auth:
                label_to_auth[label_id] = auth_id

        for label_id, entity_id in label_to_entity.items():
            if entity_id in polymer_entities and entity_id in entity_to_desc:
                auth_id = label_to_auth.get(label_id, label_id)
                if auth_id not in chain_labels:
                    chain_labels[auth_id] = entity_to_desc[entity_id]

    except Exception as e:
        print(f"  Warning: Could not parse CIF entity info: {e}")

    return chain_labels


def get_chain_labels(filepath):
    """Get chain labels from PDB or CIF; may return {} when unavailable."""
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.cif':
        return get_chain_labels_from_cif(filepath)
    if ext == '.pdb':
        return get_chain_labels_from_pdb(filepath)
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


def run_foldx_repair(fx_path, pdb_path, work_dir):
    """Run FoldX RepairPDB command."""
    pdb_name = os.path.basename(pdb_path)
    shutil.copy(pdb_path, os.path.join(work_dir, pdb_name))

    cmd = [fx_path, '--command=RepairPDB', f'--pdb={pdb_name}']
    print(f"Running FoldX RepairPDB on {pdb_name}...")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir, timeout=3600)

    basename = os.path.splitext(pdb_name)[0]
    repaired = os.path.join(work_dir, f"{basename}_Repair.pdb")

    if os.path.exists(repaired):
        print(f"  Repaired structure: {repaired}")
        return repaired
    else:
        print(f"  Error: RepairPDB failed")
        if result.stderr:
            print(f"  {result.stderr[:500]}")
        return None


def fake_foldx_repair(pdb_path, work_dir):
    """Fake FoldX repair by copying with _Repair.pdb suffix."""
    basename = os.path.splitext(os.path.basename(pdb_path))[0]
    repaired = os.path.join(work_dir, f"{basename}_Repair.pdb")
    shutil.copy(pdb_path, repaired)
    print(f"Faked repair: {repaired}")
    return repaired


def parse_foldx_interactions(fxout_path):
    """Parse ALL interaction energies from FoldX AnalyseComplex output.

    Returns dict of (chain1, chain2) -> interaction_energy.
    """
    results = {}
    with open(fxout_path, 'r') as f:
        header_found = False
        for line in f:
            if line.strip().startswith('Pdb'):
                header_found = True
                continue
            if header_found:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    try:
                        chain1 = parts[1]
                        chain2 = parts[2]
                        energy = float(parts[5])
                        pair = tuple(sorted([chain1, chain2]))
                        results[pair] = energy
                    except (ValueError, IndexError):
                        continue
    return results


def run_foldx_analyse_complex(fx_path, pdb_path, work_dir):
    """Run FoldX AnalyseComplex on full structure (all chain pairs at once)."""
    pdb_name = os.path.basename(pdb_path)
    work_pdb = os.path.join(work_dir, pdb_name)
    if not os.path.exists(work_pdb):
        shutil.copy(pdb_path, work_pdb)

    cmd = [fx_path, '--command=AnalyseComplex', f'--pdb={pdb_name}']
    print(f"Running FoldX AnalyseComplex on {pdb_name}...")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir, timeout=3600)

    # Find the Interaction output file (named Interaction_{basename}_*.fxout)
    basename = os.path.splitext(pdb_name)[0]
    import glob
    fxout_files = glob.glob(os.path.join(work_dir, f"Interaction_{basename}*.fxout"))

    if fxout_files:
        fxout = fxout_files[0]
        return parse_foldx_interactions(fxout)

    print(f"  Warning: No FoldX AnalyseComplex output found")
    if result.stderr:
        print(f"  {result.stderr[:500]}")
    return {}


def calculate_foldx_scores(fx_path, repaired_pdb, interactions, work_dir):
    """Run FoldX AnalyseComplex once and extract energies for interacting pairs."""
    all_energies = run_foldx_analyse_complex(fx_path, repaired_pdb, work_dir)

    foldx_data = {}
    for (chain1, chain2) in interactions.keys():
        pair = tuple(sorted([chain1, chain2]))
        energy = all_energies.get(pair)
        foldx_data[pair] = energy
        if energy is not None:
            print(f"  {chain1}-{chain2}: IE = {energy:.2f} kcal/mol")
        else:
            print(f"  {chain1}-{chain2}: No FoldX data")

    print(f"\n  Scored {len([v for v in foldx_data.values() if v is not None])}/{len(interactions)} pairs")
    return foldx_data


def run_foldx_buildmodel(fx_path, pdb_path, mutant_file, work_dir):
    """Run FoldX BuildModel with mutations."""
    pdb_name = os.path.basename(pdb_path)
    work_pdb = os.path.join(work_dir, pdb_name)
    if not os.path.exists(work_pdb):
        shutil.copy(pdb_path, work_pdb)

    shutil.copy(mutant_file, os.path.join(work_dir, 'individual_list.txt'))

    cmd = [fx_path, '--command=BuildModel', f'--pdb={pdb_name}',
           '--mutant-file=individual_list.txt']
    print(f"Running FoldX BuildModel...")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir, timeout=3600)

    basename = os.path.splitext(pdb_name)[0]
    wt_pdb = os.path.join(work_dir, f"WT_{basename}_1.pdb")
    mut_pdb = os.path.join(work_dir, f"{basename}_1.pdb")
    dif_file = os.path.join(work_dir, f"Dif_{basename}.fxout")

    if os.path.exists(wt_pdb) and os.path.exists(mut_pdb):
        print(f"  WT structure: {wt_pdb}")
        print(f"  Mutant structure: {mut_pdb}")

        total_ddg = None
        if os.path.exists(dif_file):
            with open(dif_file, 'r') as f:
                for line in f:
                    if line.startswith('Pdb'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        try:
                            total_ddg = float(parts[1])
                        except ValueError:
                            pass

        return wt_pdb, mut_pdb, total_ddg
    else:
        print(f"  Error: BuildModel failed")
        if result.stderr:
            print(f"  {result.stderr[:500]}")
        return None, None, None


def calculate_foldx_ddg(fx_path, wt_pdb, mut_pdb, interactions, work_dir):
    """Calculate binding ddG by running AnalyseComplex on WT and mutant."""
    wt_dir = os.path.join(work_dir, 'wt')
    mut_dir = os.path.join(work_dir, 'mut')
    os.makedirs(wt_dir, exist_ok=True)
    os.makedirs(mut_dir, exist_ok=True)

    print("\nRunning AnalyseComplex on WT structure...")
    wt_energies = run_foldx_analyse_complex(fx_path, wt_pdb, wt_dir)
    print("Running AnalyseComplex on mutant structure...")
    mut_energies = run_foldx_analyse_complex(fx_path, mut_pdb, mut_dir)

    ddg_data = {}
    print(f"\nCalculating binding ddG for {len(interactions)} pairs...")

    for (chain1, chain2) in interactions.keys():
        pair = tuple(sorted([chain1, chain2]))
        wt_energy = wt_energies.get(pair)
        mut_energy = mut_energies.get(pair)

        if wt_energy is not None and mut_energy is not None:
            ddg = mut_energy - wt_energy
            ddg_data[pair] = {
                'wt_energy': wt_energy,
                'mut_energy': mut_energy,
                'ddG': ddg
            }
            print(f"  {chain1}-{chain2}: ddG = {ddg:.2f} (WT: {wt_energy:.2f}, Mut: {mut_energy:.2f})")
        else:
            ddg_data[pair] = {
                'wt_energy': wt_energy,
                'mut_energy': mut_energy,
                'ddG': None
            }
            print(f"  {chain1}-{chain2}: Failed")

    return ddg_data


def run_prodigy_one_vs_all(pdb_path, chain_id, other_chain_ids):
    """Run PRODIGY for one chain vs all others combined."""
    others_str = ','.join(sorted(other_chain_ids))
    try:
        cmd = ['prodigy', pdb_path, '--selection', chain_id, others_str]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        output = result.stdout + result.stderr

        dg_match = re.search(r'Predicted binding affinity.*?:\s*([-\d.]+)', output)
        kd_match = re.search(r'Predicted dissociation constant.*?:\s*([\d.eE+-]+)', output)

        dg = float(dg_match.group(1)) if dg_match else None
        kd = kd_match.group(1) if kd_match else None

        return dg, kd
    except subprocess.TimeoutExpired:
        print(f"  Warning: PRODIGY timeout for {chain_id} vs all")
        return None, None
    except FileNotFoundError:
        print("  Error: PRODIGY not found. Install with: pip install prodigy-prot")
        return None, None
    except Exception as e:
        print(f"  Warning: PRODIGY error for {chain_id} vs all: {e}")
        return None, None


def calculate_one_vs_all(pdb_path, protein_chain_ids):
    """Calculate PRODIGY one-vs-all for each chain."""
    chain_list = sorted(protein_chain_ids)
    results = {}
    total = len(chain_list)
    print(f"\nCalculating one-vs-all binding energies for {total} chains...")

    for i, chain_id in enumerate(chain_list):
        other_chains = [c for c in chain_list if c != chain_id]
        if not other_chains:
            continue
        others_str = ','.join(other_chains)
        print(f"  [{i+1}/{total}] {chain_id} vs {others_str}...", end=" ")

        dg, kd = run_prodigy_one_vs_all(pdb_path, chain_id, other_chains)
        results[chain_id] = {'dG': dg, 'Kd': kd, 'vs': others_str}

        if dg is not None:
            print(f"ΔG = {dg:.1f} kcal/mol")
        else:
            print("Failed")

    return results


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
            'chain1': n1,
            'chain2': n2,
            'chain1_label': label1,
            'chain2_label': label2,
            'contacts': weight,
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
        #graph-canvas {{
            cursor: default;
        }}
        #graph-canvas.edge-hover {{
            cursor: pointer;
        }}
        #complex-modal {{
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0,0,0,0.6);
            z-index: 2000;
            justify-content: center;
            align-items: center;
        }}
        #complex-modal.visible {{
            display: flex;
        }}
        #complex-modal-content {{
            background: white;
            border-radius: 8px;
            width: 900px;
            max-width: 95%;
            height: 700px;
            max-height: 95%;
            position: relative;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            box-shadow: 0 10px 40px rgba(0,0,0,0.4);
        }}
        #complex-modal-header {{
            padding: 12px 20px;
            border-bottom: 1px solid #ddd;
            background: #f9f9f9;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        #complex-modal-title {{
            font-weight: bold;
            font-size: 15px;
            color: #333;
        }}
        #complex-modal-close {{
            background: #e6194b;
            color: white;
            border: none;
            border-radius: 4px;
            width: 30px;
            height: 30px;
            cursor: pointer;
            font-size: 18px;
            font-weight: bold;
            line-height: 1;
        }}
        #complex-modal-close:hover {{
            background: #c4003a;
        }}
        #complex-viewer {{
            flex: 1;
            position: relative;
            background: white;
        }}
        #complex-modal-footer {{
            padding: 10px 20px;
            border-top: 1px solid #ddd;
            background: #f9f9f9;
            font-size: 12px;
            color: #555;
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
        }}
        #complex-modal-footer span b {{
            color: #0066cc;
        }}
        .chain-swatch {{
            display: inline-block;
            width: 12px;
            height: 12px;
            border-radius: 2px;
            margin-right: 4px;
            vertical-align: middle;
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
        Nodes: {num_nodes} chains | Edges: {num_edges} interactions | Drag nodes to rearrange | Double-click a node to view its chain | Click an edge to view the pair complex
    </div>
    <div id="tooltip"></div>

    <div id="complex-modal">
        <div id="complex-modal-content">
            <div id="complex-modal-header">
                <div id="complex-modal-title">Complex</div>
                <button id="complex-modal-close" title="Close">&times;</button>
            </div>
            <div id="complex-viewer"></div>
            <div id="complex-modal-footer"></div>
        </div>
    </div>

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

                div.addEventListener('dblclick', (e) => {{
                    e.preventDefault();
                    e.stopPropagation();
                    showSingleChain(node);
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

        // ===== Edge click → show pair complex in modal =====
        const modal = document.getElementById('complex-modal');
        const modalTitle = document.getElementById('complex-modal-title');
        const modalFooter = document.getElementById('complex-modal-footer');
        const modalClose = document.getElementById('complex-modal-close');
        const complexViewerEl = document.getElementById('complex-viewer');
        let complexViewer = null;

        function pointToSegmentDistance(px, py, x1, y1, x2, y2) {{
            const dx = x2 - x1;
            const dy = y2 - y1;
            const lenSq = dx * dx + dy * dy;
            if (lenSq === 0) {{
                const ex = px - x1, ey = py - y1;
                return Math.sqrt(ex * ex + ey * ey);
            }}
            let t = ((px - x1) * dx + (py - y1) * dy) / lenSq;
            t = Math.max(0, Math.min(1, t));
            const projX = x1 + t * dx;
            const projY = y1 + t * dy;
            const ex = px - projX, ey = py - projY;
            return Math.sqrt(ex * ex + ey * ey);
        }}

        function findEdgeAt(mx, my) {{
            let closest = null;
            let closestDist = Infinity;
            edges.forEach(edge => {{
                const n1 = nodes.find(n => n.id === edge.chain1);
                const n2 = nodes.find(n => n.id === edge.chain2);
                if (!n1 || !n2) return;
                const x1 = n1.x + 50, y1 = n1.y + 50;
                const x2 = n2.x + 50, y2 = n2.y + 50;
                const threshold = Math.max(edge.width, 4) + 4;
                const d = pointToSegmentDistance(mx, my, x1, y1, x2, y2);
                if (d <= threshold && d < closestDist) {{
                    closest = edge;
                    closestDist = d;
                }}
            }});
            return closest;
        }}

        function showSingleChain(node) {{
            const colorA = '#4363d8';
            modalTitle.textContent = 'Chain ' + node.id;
            modalFooter.innerHTML =
                '<span><span class="chain-swatch" style="background:' + colorA + '"></span>' +
                '<b>' + node.id + '</b>: ' + (node.label || node.id) + '</span>' +
                '<span>Interactions: <b>' + node.degree + '</b></span>' +
                '<span>Total contacts: <b>' + node.contacts + '</b></span>';

            modal.classList.add('visible');

            complexViewerEl.innerHTML = '';
            complexViewer = $3Dmol.createViewer(complexViewerEl, {{ backgroundColor: 'white' }});
            if (node.pdb) {{
                complexViewer.addModel(node.pdb.replace(/\\\\n/g, '\\n'), 'pdb');
                complexViewer.setStyle({{chain: node.id}}, {{cartoon: {{color: colorA}}}});
            }}
            complexViewer.zoomTo();
            complexViewer.render();
            setTimeout(() => {{ if (complexViewer) complexViewer.resize(); }}, 50);
        }}

        function showComplex(edge) {{
            const colorA = '#4363d8';
            const colorB = '#e6194b';
            modalTitle.textContent = 'Complex ' + edge.chain1 + ' – ' + edge.chain2;
            modalFooter.innerHTML =
                '<span><span class="chain-swatch" style="background:' + colorA + '"></span>' +
                '<b>' + edge.chain1 + '</b>: ' + edge.chain1_label + '</span>' +
                '<span><span class="chain-swatch" style="background:' + colorB + '"></span>' +
                '<b>' + edge.chain2 + '</b>: ' + edge.chain2_label + '</span>' +
                '<span>Contacts: <b>' + edge.contacts + '</b></span>' +
                '<span>ΔG: <b>' + edge.dG + '</b> kcal/mol</span>' +
                '<span>Kd: <b>' + edge.Kd + '</b></span>';

            modal.classList.add('visible');

            complexViewerEl.innerHTML = '';
            complexViewer = $3Dmol.createViewer(complexViewerEl, {{ backgroundColor: 'white' }});

            const nodeA = nodes.find(n => n.id === edge.chain1);
            const nodeB = nodes.find(n => n.id === edge.chain2);
            if (nodeA && nodeA.pdb) {{
                complexViewer.addModel(nodeA.pdb.replace(/\\\\n/g, '\\n'), 'pdb');
            }}
            if (nodeB && nodeB.pdb) {{
                complexViewer.addModel(nodeB.pdb.replace(/\\\\n/g, '\\n'), 'pdb');
            }}
            complexViewer.setStyle({{chain: edge.chain1}}, {{cartoon: {{color: colorA}}}});
            complexViewer.setStyle({{chain: edge.chain2}}, {{cartoon: {{color: colorB}}}});
            complexViewer.zoomTo();
            complexViewer.render();
            setTimeout(() => {{ if (complexViewer) complexViewer.resize(); }}, 50);
        }}

        function closeModal() {{
            modal.classList.remove('visible');
            if (complexViewer) {{
                try {{ complexViewer.clear(); }} catch (e) {{}}
                complexViewer = null;
            }}
            complexViewerEl.innerHTML = '';
        }}

        modalClose.addEventListener('click', closeModal);
        modal.addEventListener('click', (e) => {{ if (e.target === modal) closeModal(); }});
        document.addEventListener('keydown', (e) => {{
            if (e.key === 'Escape' && modal.classList.contains('visible')) closeModal();
        }});

        canvas.addEventListener('click', (e) => {{
            const rect = canvas.getBoundingClientRect();
            const scaleX = canvas.width / rect.width;
            const scaleY = canvas.height / rect.height;
            const mx = (e.clientX - rect.left) * scaleX;
            const my = (e.clientY - rect.top) * scaleY;
            const edge = findEdgeAt(mx, my);
            if (edge) showComplex(edge);
        }});

        canvas.addEventListener('mousemove', (e) => {{
            const rect = canvas.getBoundingClientRect();
            const scaleX = canvas.width / rect.width;
            const scaleY = canvas.height / rect.height;
            const mx = (e.clientX - rect.left) * scaleX;
            const my = (e.clientY - rect.top) * scaleY;
            const edge = findEdgeAt(mx, my);
            if (edge) {{
                canvas.classList.add('edge-hover');
                tooltip.innerHTML = '<b>' + edge.name + '</b><br>' +
                    edge.chain1_label + ' &harr; ' + edge.chain2_label + '<br>' +
                    'Contacts: ' + edge.contacts + '<br>' +
                    'ΔG: ' + edge.dG + ' kcal/mol<br>' +
                    '<i>Click to view complex</i>';
                tooltip.style.display = 'block';
                tooltip.style.left = (e.pageX + 15) + 'px';
                tooltip.style.top = (e.pageY + 15) + 'px';
            }} else {{
                canvas.classList.remove('edge-hover');
                if (!draggedNode) tooltip.style.display = 'none';
            }}
        }});

        canvas.addEventListener('mouseleave', () => {{
            canvas.classList.remove('edge-hover');
            tooltip.style.display = 'none';
        }});

        // Also allow legend items to open the complex
        setTimeout(() => {{
            const items = document.querySelectorAll('.legend-item');
            items.forEach((item, idx) => {{
                const edge = sortedEdges[idx];
                if (!edge) return;
                item.addEventListener('click', () => showComplex(edge));
            }});
        }}, 200);
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


def save_foldx_ranking(foldx_data, chain_labels, output_file):
    """Save FoldX interaction energy ranking."""
    valid = [(k, v) for k, v in foldx_data.items() if v is not None]
    sorted_pairs = sorted(valid, key=lambda x: x[1])

    with open(output_file, 'w') as f:
        f.write("# FoldX Interaction Energies (AnalyseComplex)\n")
        f.write("# More negative = stronger binding\n")
        f.write("=" * 90 + "\n\n")
        f.write(f"{'Rank':<6}{'Chain Pair':<15}{'IE (kcal/mol)':<18}{'Chain 1':<30}{'Chain 2':<30}\n")
        f.write("-" * 90 + "\n")

        for rank, ((c1, c2), energy) in enumerate(sorted_pairs, 1):
            l1 = chain_labels.get(c1, "Unknown")[:28]
            l2 = chain_labels.get(c2, "Unknown")[:28]
            f.write(f"{rank:<6}{c1}-{c2:<13}{energy:<18.2f}{l1:<30}{l2:<30}\n")

        f.write(f"\nTotal pairs: {len(sorted_pairs)}\n")

    print(f"Saved FoldX ranking: {output_file}")


def save_foldx_ddg_ranking(ddg_data, chain_labels, output_file):
    """Save FoldX ddG (mutation effect on binding) ranking."""
    valid = [(k, v) for k, v in ddg_data.items() if v.get('ddG') is not None]
    sorted_pairs = sorted(valid, key=lambda x: x[1]['ddG'])

    with open(output_file, 'w') as f:
        f.write("# FoldX Binding ddG (BuildModel + AnalyseComplex)\n")
        f.write("# Negative ddG = mutation stabilizes binding\n")
        f.write("# Positive ddG = mutation destabilizes binding\n")
        f.write("=" * 120 + "\n\n")
        f.write(f"{'Rank':<6}{'Pair':<12}{'ddG':<12}{'WT IE':<14}{'Mut IE':<14}{'Chain 1':<25}{'Chain 2':<25}\n")
        f.write("-" * 120 + "\n")

        for rank, ((c1, c2), data) in enumerate(sorted_pairs, 1):
            l1 = chain_labels.get(c1, "Unknown")[:23]
            l2 = chain_labels.get(c2, "Unknown")[:23]
            f.write(f"{rank:<6}{c1}-{c2:<10}{data['ddG']:<12.2f}{data['wt_energy']:<14.2f}{data['mut_energy']:<14.2f}{l1:<25}{l2:<25}\n")

        f.write(f"\nTotal pairs: {len(sorted_pairs)}\n")

    print(f"Saved FoldX ddG ranking: {output_file}")


def save_one_vs_all_ranking(ova_data, chain_labels, output_file):
    """Save one-vs-all PRODIGY ranking."""
    valid = [(k, v) for k, v in ova_data.items() if v.get('dG') is not None]
    sorted_chains = sorted(valid, key=lambda x: x[1]['dG'])

    with open(output_file, 'w') as f:
        f.write("# PRODIGY One-vs-All Binding Energies\n")
        f.write("# Each chain scored against all other chains combined\n")
        f.write("=" * 90 + "\n\n")
        f.write(f"{'Rank':<6}{'Chain':<10}{'ΔG (kcal/mol)':<18}{'Kd (M)':<18}{'Description':<35}\n")
        f.write("-" * 90 + "\n")

        for rank, (chain_id, data) in enumerate(sorted_chains, 1):
            label = chain_labels.get(chain_id, "Unknown")[:33]
            kd = data.get('Kd', 'N/A')
            f.write(f"{rank:<6}{chain_id:<10}{data['dG']:<18.2f}{kd:<18}{label:<35}\n")

        f.write(f"\nTotal chains: {len(sorted_chains)}\n")

    print(f"Saved one-vs-all ranking: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate PPI graph with 3D structure nodes and binding energy (PRODIGY)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python ppi_graph_3d_dg.py structure.pdb
    python ppi_graph_3d_dg.py structure.cif
    python ppi_graph_3d_dg.py structure.pdb --cutoff 4.0
    python ppi_graph_3d_dg.py structure.pdb --fx_score --fx_path /path/to/foldx
    python ppi_graph_3d_dg.py structure.pdb --fx_mut mutations.txt --fx_path /path/to/foldx
    python ppi_graph_3d_dg.py structure.pdb --one_vs_all
        """
    )

    parser.add_argument('input_file', help='Input PDB or CIF file (CIF is converted to PDB via gemmi)')
    parser.add_argument('--cutoff', type=float, default=5.0,
                        help='Distance cutoff for interactions in Angstroms (default: 5.0)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--skip-prodigy', action='store_true',
                        help='Skip PRODIGY calculation (for testing)')
    parser.add_argument('--fx_score', action='store_true',
                        help='Run FoldX RepairPDB + AnalyseComplex to score chain interactions')
    parser.add_argument('--fx_path', type=str, default=None,
                        help='Path to FoldX executable')
    parser.add_argument('--fx_mut', type=str, default=None,
                        help='Path to individual_list.txt for FoldX BuildModel (skips repair, calculates binding ddG)')
    parser.add_argument('--one_vs_all', action='store_true',
                        help='Calculate PRODIGY binding energy for each chain vs all others combined')

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Error: File not found: {args.input_file}")
        return 1

    ext = os.path.splitext(args.input_file)[1].lower()
    if ext not in ('.pdb', '.cif'):
        print(f"Error: Only .pdb and .cif files are supported. Got: {ext}")
        return 1

    if (args.fx_score or args.fx_mut) and not args.fx_path:
        print("Error: --fx_path is required when using --fx_score or --fx_mut")
        return 1

    if args.fx_mut and not os.path.exists(args.fx_mut):
        print(f"Error: Mutation file not found: {args.fx_mut}")
        return 1

    os.makedirs(args.output_dir, exist_ok=True)

    basename = os.path.splitext(os.path.basename(args.input_file))[0]

    # Create complexes directory for PRODIGY PDB files
    complexes_dir = os.path.join(args.output_dir, f"{basename}_complexes")

    print(f"Processing: {args.input_file}")
    print(f"Distance cutoff: {args.cutoff} A")
    print()

    # Extract chain labels from the original file (CIF entity info, if any)
    print("Extracting chain labels...")
    chain_labels = get_chain_labels(args.input_file)

    # Convert CIF to PDB so PRODIGY / FoldX / BioPython PDBIO work downstream
    pdb_path = args.input_file
    if ext == '.cif':
        pdb_path = os.path.join(args.output_dir, f"{basename}_from_cif.pdb")
        print(f"Converting CIF to PDB via gemmi: {pdb_path}")
        convert_cif_to_pdb(args.input_file, pdb_path)

    print("Parsing structure...")
    structure = parse_structure(pdb_path)

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

    # FoldX scoring
    foldx_data = {}
    foldx_ddg_data = {}
    foldx_dir = os.path.join(args.output_dir, f"{basename}_foldx")

    if args.fx_score and interactions:
        os.makedirs(foldx_dir, exist_ok=True)
        repaired_pdb = run_foldx_repair(args.fx_path, pdb_path, foldx_dir)
        if repaired_pdb:
            foldx_data = calculate_foldx_scores(args.fx_path, repaired_pdb, interactions, foldx_dir)

    if args.fx_mut and interactions:
        os.makedirs(foldx_dir, exist_ok=True)
        repaired_pdb = fake_foldx_repair(pdb_path, foldx_dir)
        wt_pdb, mut_pdb, total_ddg = run_foldx_buildmodel(args.fx_path, repaired_pdb, args.fx_mut, foldx_dir)
        if total_ddg is not None:
            print(f"  Total stability ddG: {total_ddg:.2f} kcal/mol")
        if wt_pdb and mut_pdb:
            foldx_ddg_data = calculate_foldx_ddg(args.fx_path, wt_pdb, mut_pdb, interactions, foldx_dir)

    # One-vs-all PRODIGY
    ova_data = {}
    if args.one_vs_all:
        ova_data = calculate_one_vs_all(pdb_path, protein_chain_ids)

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

    if foldx_data:
        foldx_file = os.path.join(args.output_dir, f"{basename}_foldx_scores.txt")
        save_foldx_ranking(foldx_data, chain_labels, foldx_file)

    if foldx_ddg_data:
        ddg_file = os.path.join(args.output_dir, f"{basename}_foldx_ddg.txt")
        save_foldx_ddg_ranking(foldx_ddg_data, chain_labels, ddg_file)

    if ova_data:
        ova_file = os.path.join(args.output_dir, f"{basename}_one_vs_all.txt")
        save_one_vs_all_ranking(ova_data, chain_labels, ova_file)

    print("\nDone!")
    print(f"\nSummary:")
    print(f"  Protein chains: {len(protein_chain_ids)}")
    print(f"  Interacting pairs: {len(interactions)}")
    print(f"  Total contacts: {sum(interactions.values())}")
    print(f"  Pairs with binding data: {len([v for v in binding_data.values() if v.get('dG') is not None])}")
    if foldx_data:
        print(f"  FoldX scored pairs: {len([v for v in foldx_data.values() if v is not None])}")
    if foldx_ddg_data:
        print(f"  FoldX ddG pairs: {len([v for v in foldx_ddg_data.values() if v.get('ddG') is not None])}")
    if ova_data:
        print(f"  One-vs-all chains: {len([v for v in ova_data.values() if v.get('dG') is not None])}")
    print(f"\nOutput files:")
    print(f"  HTML visualization: {html_file}")
    print(f"  Chain info: {chain_info_file}")
    print(f"  Residue contacts: {contacts_file}")
    if binding_data:
        print(f"  Binding strength: {binding_file}")
        print(f"  Complex PDBs: {complexes_dir}/")
    if foldx_data:
        print(f"  FoldX scores: {os.path.join(args.output_dir, f'{basename}_foldx_scores.txt')}")
    if foldx_ddg_data:
        print(f"  FoldX ddG: {os.path.join(args.output_dir, f'{basename}_foldx_ddg.txt')}")
    if ova_data:
        print(f"  One-vs-all: {os.path.join(args.output_dir, f'{basename}_one_vs_all.txt')}")

    return 0


if __name__ == '__main__':
    exit(main())
