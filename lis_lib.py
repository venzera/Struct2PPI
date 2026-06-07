#!/usr/bin/env python3
"""
lis_lib -- Local Interaction Score (LIS / cLIS / iLIS) computation for Struct2PPI
================================================================================
Trimmed, single-model adaptation of the AFM-LIS reference implementation
(Kim & Perrimon, LIVIA, bioRxiv 2026) for use inside Struct2PPI's
``--predicted`` mode.

Compared to the upstream ``lis.py`` this module deliberately drops:
  * parallel / multiprocessing
  * zip / tar / gz / xz archive handling (plain folders / files only)
  * ChimeraX (.cxc) script generation

It keeps the core scoring (LIS, cLIS, iLIS, LIA, cLIA, iLIA, ipTM, pLDDT) and
adds per-chain-pair PAE / contact-filtered-LIA maps for visualization.

Supported platforms (auto-detected): AlphaFold3, Boltz, Chai-1 (+ generic).

Reference:
    Kim, Ah-Ram, and Norbert Perrimon. "LIVIA: a browser-based tool for
    assessing and visualizing predicted protein interactions."
    bioRxiv (2026). doi:10.64898/2026.05.01.721633
    iLIS metric: Kim A.-R. et al. "FlyPredictome." bioRxiv (2026).
    doi:10.64898/2026.04.14.718529
"""

import io
import json
import math
import os
import re
from collections import OrderedDict, defaultdict

import numpy as np
from scipy.spatial.distance import pdist, squareform

# ============================================================================
# Constants
# ============================================================================

ION_NAMES = {'ZN', 'CA', 'MG', 'MN', 'FE', 'CU', 'NA', 'K', 'CL', 'NI', 'CO', 'CD'}

# iLIS positive-interaction threshold (Kim et al., 2026)
ILIS_THRESHOLD = 0.223


# ============================================================================
# File scanning (folders / single files only -- no archives)
# ============================================================================

def scan_files(path):
    """Return (filenames, read_fn). filenames are relative-ish paths.

    read_fn(name) returns str for text files, bytes for binary (.npy/.npz).
    """
    filenames = []
    if os.path.isdir(path):
        root = path
        for dirpath, _dirs, files in os.walk(path):
            for f in files:
                if f.startswith('._') or '__MACOSX' in dirpath:
                    continue
                full = os.path.join(dirpath, f)
                filenames.append(os.path.relpath(full, root))
        base = root
    else:
        root = os.path.dirname(path) or '.'
        base = root
        filenames = [os.path.basename(path)]

    def read_fn(name):
        full = os.path.join(base, name)
        binary = name.endswith(('.npy', '.npz'))
        try:
            if binary:
                with open(full, 'rb') as fh:
                    return fh.read()
            with open(full, 'r') as fh:
                return fh.read()
        except (OSError, UnicodeDecodeError):
            return None

    return filenames, read_fn


# ============================================================================
# Platform detection
# ============================================================================

def detect_platform(filenames):
    """Auto-detect AlphaFold3 / Boltz / Chai from filenames (generic fallback)."""
    basenames = [os.path.basename(f) for f in filenames]

    # AlphaFold3 (server + standard): *_model_N.cif + *_full_data_N.json
    has_af3_model = any(re.search(r'_model_\d+\.cif$', b) for b in basenames)
    has_af3_full = any(re.search(r'_full_data_\d+\.json$', b) for b in basenames)
    if has_af3_model and has_af3_full:
        return 'alphafold3'

    # Boltz
    has_boltz_conf = any(b.startswith('confidence_') and b.endswith('.json') for b in basenames)
    has_boltz_struct = any(b.endswith('.cif') or b.endswith('.pdb') for b in basenames)
    if has_boltz_conf and has_boltz_struct:
        return 'boltz'

    # Chai-1
    has_chai_model = any(re.match(r'^pred\.(rank_\d+|model_idx_\d+)\.cif$', b) for b in basenames)
    has_chai_scores = any(re.match(r'^scores\.(rank_\d+|model_idx_\d+)\.(json|npz)$', b) for b in basenames)
    if has_chai_model or (has_chai_scores and any(b.startswith('pred.') for b in basenames)):
        return 'chai'

    return 'generic'


# ============================================================================
# Model discovery -- yields (name, rank, struct_path, pae_path, scores_path, fmt)
# ============================================================================

def find_models(filenames, platform):
    filenames = [f for f in filenames if '__MACOSX' not in f and not os.path.basename(f).startswith('._')]
    basenames_map = {os.path.basename(f): f for f in filenames}

    if platform == 'alphafold3':
        yield from _find_af3(filenames, basenames_map)
    elif platform == 'boltz':
        yield from _find_boltz(filenames, basenames_map)
    elif platform == 'chai':
        yield from _find_chai(filenames, basenames_map)
    else:
        yield from _find_generic(filenames, basenames_map)


def _toplevel_name(filepath):
    parts = filepath.split(os.sep)
    if len(parts) >= 2 and parts[0] and not parts[0].startswith('__'):
        return parts[0]
    return None


def _find_af3(filenames, basenames_map):
    """AlphaFold3: *_model_N.cif + *_full_data_N.json (+ *_summary_confidences_N.json)."""
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^(.+)_model_(\d+)\.cif$', base)
        if not m:
            continue
        prefix, idx = m.group(1), m.group(2)
        full_data = basenames_map.get(f'{prefix}_full_data_{idx}.json')
        summary = basenames_map.get(f'{prefix}_summary_confidences_{idx}.json')
        if not full_data:
            continue
        yield (prefix, idx, name, full_data, summary or full_data, 'cif')


def _find_boltz(filenames, basenames_map):
    """Boltz: *_model_N.pdb/.cif + confidence_*_model_N.json + pae_*_model_N.npz."""
    dir_groups = defaultdict(list)
    for f in filenames:
        dir_groups[os.path.dirname(f)].append(f)
    if len(dir_groups) <= 1:
        dir_groups = {'': filenames}

    for dirpath, dir_files in sorted(dir_groups.items()):
        struct_files = sorted([f for f in dir_files if f.endswith(('.cif', '.pdb'))],
                              key=lambda f: os.path.basename(f))
        conf_files = [f for f in dir_files if os.path.basename(f).startswith('confidence') and f.endswith('.json')]
        pae_files = [f for f in dir_files if os.path.basename(f).startswith('pae') and f.endswith('.npz')]
        if not struct_files:
            continue
        pred_name = os.path.basename(dirpath) if dirpath else ''
        if not pred_name or pred_name in ('Boltz1_outputs', 'boltz_outputs', 'predictions', 'result', 'results'):
            pred_name = _toplevel_name(struct_files[0]) or 'boltz'

        for i, sf in enumerate(struct_files):
            fmt = 'cif' if sf.endswith('.cif') else 'pdb'
            sb = os.path.basename(sf)
            model_idx = str(i)
            m = re.search(r'model_(\d+)', sb)
            if m:
                model_idx = m.group(1)
            conf_path = next((cf for cf in conf_files if f'model_{model_idx}' in os.path.basename(cf)), None)
            if not conf_path and i < len(conf_files):
                conf_path = conf_files[i]
            pae_path = next((pf for pf in pae_files if f'model_{model_idx}' in os.path.basename(pf)), None)
            if not pae_path and i < len(pae_files):
                pae_path = pae_files[i]
            yield (pred_name, model_idx, sf, pae_path or conf_path, conf_path, fmt)


def _find_chai(filenames, basenames_map):
    """Chai-1: pred.rank_N.cif + scores.rank_N.json + pae.rank_N.npy/.npz."""
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^pred\.(rank_(\d+)|model_idx_(\d+))\.cif$', base)
        if not m:
            continue
        rank = m.group(2) or m.group(3)
        score_path = (basenames_map.get(f'scores.rank_{rank}.json')
                      or basenames_map.get(f'scores.model_idx_{rank}.json'))
        if not score_path:
            score_path = next((fn for fn in filenames
                               if 'scores' in os.path.basename(fn) and fn.endswith('.json')), None)
        pae_path = (basenames_map.get(f'pae.rank_{rank}.npy')
                    or basenames_map.get(f'pae.rank_{rank}.npz')
                    or basenames_map.get(f'pae.model_idx_{rank}.npy')
                    or basenames_map.get(f'pae.model_idx_{rank}.npz'))
        pred_name = _toplevel_name(name) or 'prediction'
        yield (pred_name, rank, name, pae_path or score_path, score_path, 'cif')


def _find_generic(filenames, basenames_map):
    """Fallback: pair any structure with a sibling .npy/.npz/.json PAE source."""
    by_dir = defaultdict(list)
    for f in filenames:
        by_dir[os.path.dirname(f)].append(f)
    for dirpath, files in sorted(by_dir.items()):
        structs = sorted(f for f in files if f.endswith(('.cif', '.pdb')))
        paes = [f for f in files if f.endswith(('.npy', '.npz'))]
        jsons = [f for f in files if f.endswith('.json')]
        if not structs:
            continue
        pae = paes[0] if paes else (jsons[0] if jsons else None)
        if not pae:
            continue
        scores = jsons[0] if jsons else pae
        name = os.path.basename(dirpath) or os.path.splitext(os.path.basename(structs[0]))[0]
        fmt = 'cif' if structs[0].endswith('.cif') else 'pdb'
        yield (name, '0', structs[0], pae, scores, fmt)


# ============================================================================
# PAE extraction
# ============================================================================

def extract_pae(pae_source, read_fn):
    """Extract a 2D PAE matrix (numpy array) or None."""
    if not pae_source:
        return None
    content = read_fn(pae_source)
    if content is None:
        return None
    basename = os.path.basename(pae_source)

    if basename.endswith('.npy'):
        if isinstance(content, str):
            content = content.encode('utf-8')
        arr = np.load(io.BytesIO(content))
        if arr.ndim == 3:
            arr = arr[0]
        return arr.astype(np.float32)

    if basename.endswith('.npz'):
        if isinstance(content, str):
            content = content.encode('utf-8')
        npz = np.load(io.BytesIO(content))
        arr = npz['pae'] if 'pae' in npz else None
        if arr is None:
            for key in npz.files:
                if npz[key].ndim >= 2:
                    arr = npz[key]
                    break
        if arr is None:
            return None
        if arr.ndim == 3:
            arr = arr[0]
        return arr.astype(np.float32)

    if not isinstance(content, str):
        return None
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        return None

    # AF3 full_data: {pae: [[...]]}
    if 'pae' in data and isinstance(data['pae'], list):
        pae = data['pae']
        if isinstance(pae[0], list):
            return np.array(pae, dtype=np.float32)
        n = int(round(math.sqrt(len(pae))))
        if n * n == len(pae):
            return np.array(pae, dtype=np.float32).reshape(n, n)

    if isinstance(data, list) and data and 'predicted_aligned_error' in (data[0] if isinstance(data[0], dict) else {}):
        return np.array(data[0]['predicted_aligned_error'], dtype=np.float32)

    if 'predicted_aligned_error' in data:
        pae = data['predicted_aligned_error']
        if isinstance(pae, list) and isinstance(pae[0], list):
            return np.array(pae, dtype=np.float32)

    for key in ('pae_matrix', 'predicted_aligned_error_matrix'):
        if key in data:
            return np.array(data[key], dtype=np.float32)

    return None


# ============================================================================
# Confidence score extraction
# ============================================================================

def _unwrap(v):
    if isinstance(v, (list, np.ndarray)) and len(v) == 1:
        return float(v[0])
    return v


def extract_confidence_scores(confidence_path, read_fn):
    """Extract pTM / ipTM / chain_pair_iptm / pLDDT from a confidence JSON."""
    if not confidence_path:
        return {}
    content = read_fn(confidence_path)
    if not content or not isinstance(content, str):
        return {}
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        return {}

    scores = {}
    if 'ptm' in data:
        scores['pTM'] = _unwrap(data['ptm'])
    if 'iptm' in data:
        scores['ipTM'] = _unwrap(data['iptm'])
    if 'chain_pair_iptm' in data:
        scores['chainPairIptm'] = data['chain_pair_iptm']
    if 'atom_plddts' in data:
        plddts = data['atom_plddts']
        scores['pLDDT'] = sum(plddts) / len(plddts) if plddts else 0

    # ColabFold-style plddt array
    if 'plddt' in data:
        plddts = data['plddt']
        if isinstance(plddts, list) and plddts:
            scores.setdefault('pLDDT', sum(plddts) / len(plddts))

    # Chai-1
    if 'per_chain_pair_iptm' in data and 'chainPairIptm' not in scores:
        raw = data['per_chain_pair_iptm']
        if (isinstance(raw, list) and raw and isinstance(raw[0], list)
                and raw[0] and isinstance(raw[0][0], list)):
            scores['chainPairIptm'] = raw[0]
        else:
            scores['chainPairIptm'] = raw

    # Boltz
    if 'ptm_score' in data:
        scores['pTM'] = _unwrap(data['ptm_score'])
    if 'iptm_score' in data:
        scores['ipTM'] = _unwrap(data['iptm_score'])
    if 'pair_chains_iptm' in data and 'chainPairIptm' not in scores:
        raw = data['pair_chains_iptm']
        if isinstance(raw, list):
            scores['chainPairIptm'] = raw
        elif isinstance(raw, dict):
            keys = sorted(raw.keys(), key=lambda x: int(x))
            matrix = []
            for rk in keys:
                row = raw[rk]
                ck = sorted(row.keys(), key=lambda x: int(x))
                matrix.append([row[c] for c in ck])
            scores['chainPairIptm'] = matrix
    if 'complex_plddt' in data:
        scores.setdefault('pLDDT', data['complex_plddt'])

    return scores


# ============================================================================
# Structure parsing -- coordinates, chains, B-factors
# ============================================================================

def parse_pdb_coords(pdb_text):
    """One Cb (Ca for GLY, P for nucleic) coordinate per residue from PDB text."""
    residues = OrderedDict()
    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        comp_id = line[17:20].strip()
        chain = line[21:22].strip() or 'A'
        try:
            resnum = int(line[22:26].strip())
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        key = f'{chain}:{resnum}'
        if atom_name == 'CB':
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
        elif atom_name == 'CA' and comp_id == 'GLY' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
        elif atom_name == 'P' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': True}
    return list(residues.values())


def _iter_cif_atom_site(cif_text):
    """Yield dict-per-atom for the _atom_site loop of an mmCIF."""
    in_atom_site = False
    col_names = []
    for line in cif_text.split('\n'):
        if line.startswith('_atom_site.'):
            in_atom_site = True
            col_names.append(line.strip().split('.')[1])
            continue
        if in_atom_site and not line.startswith('_atom_site.') and not line.startswith('#') and line.strip():
            if line.startswith('loop_') or line.startswith('_'):
                in_atom_site = False
                continue
            parts = line.strip().split()
            if len(parts) < len(col_names):
                continue
            yield {c: parts[i] for i, c in enumerate(col_names)}


def parse_cif_coords(cif_text):
    """One Cb (Ca for GLY, P for nucleic) coordinate per residue from mmCIF text."""
    residues = OrderedDict()
    for row in _iter_cif_atom_site(cif_text):
        group_pdb = row.get('group_PDB', '')
        if group_pdb not in ('ATOM', 'HETATM'):
            continue
        atom_name = row.get('label_atom_id', '')
        comp_id = row.get('label_comp_id', '')
        chain = row.get('label_asym_id', '')
        try:
            x = float(row.get('Cartn_x')); y = float(row.get('Cartn_y')); z = float(row.get('Cartn_z'))
        except (ValueError, TypeError):
            continue
        if group_pdb == 'HETATM' and comp_id in ION_NAMES:
            residues.setdefault(f'{chain}:1', {'chain': chain, 'resnum': 1, 'x': x, 'y': y, 'z': z, 'has_p': False})
            continue
        try:
            resnum = int(row.get('label_seq_id'))
        except (ValueError, TypeError):
            continue
        key = f'{chain}:{resnum}'
        if atom_name == 'CB':
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
        elif atom_name == 'CA' and comp_id == 'GLY' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
        elif atom_name == 'P' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': True}
    return list(residues.values())


def parse_structure_coords(text, fmt):
    return parse_pdb_coords(text) if fmt == 'pdb' else parse_cif_coords(text)


def get_chains_from_pdb(pdb_text):
    chain_order = []
    chain_counts = OrderedDict()
    seen = set()
    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM'):
            continue
        if line[12:16].strip() not in ('CA', 'P'):
            continue
        chain = line[21:22].strip() or 'A'
        rkey = f'{chain}:{line[22:26].strip()}'
        if rkey in seen:
            continue
        seen.add(rkey)
        if chain not in chain_counts:
            chain_order.append(chain)
            chain_counts[chain] = 0
        chain_counts[chain] += 1
    return {'names': chain_order, 'sizes': [chain_counts[c] for c in chain_order],
            'types': ['protein'] * len(chain_order)}


def get_chains_from_cif(cif_text):
    chain_order = []
    chain_counts = OrderedDict()
    chain_types = {}
    seen = set()
    for row in _iter_cif_atom_site(cif_text):
        group_pdb = row.get('group_PDB', '')
        atom_name = row.get('label_atom_id', '')
        comp_id = row.get('label_comp_id', '')
        chain = row.get('label_asym_id', '')
        res_seq = row.get('label_seq_id', '')
        counted = False
        if group_pdb == 'ATOM' and atom_name in ('CA', 'P'):
            rkey = f'{chain}:{res_seq}'
            if rkey not in seen:
                seen.add(rkey)
                counted = True
                if atom_name == 'P':
                    chain_types[chain] = 'dna' if comp_id.startswith('D') else 'rna'
                elif chain not in chain_types:
                    chain_types[chain] = 'protein'
        elif group_pdb == 'HETATM' and comp_id in ION_NAMES:
            rkey = f'{chain}:ion'
            if rkey not in seen:
                seen.add(rkey)
                counted = True
                chain_types[chain] = 'ion'
        if counted:
            if chain not in chain_counts:
                chain_order.append(chain)
                chain_counts[chain] = 0
            chain_counts[chain] += 1
    return {'names': chain_order, 'sizes': [chain_counts[c] for c in chain_order],
            'types': [chain_types.get(c, 'protein') for c in chain_order]}


def get_chains_from_structure(text, fmt):
    return get_chains_from_pdb(text) if fmt == 'pdb' else get_chains_from_cif(text)


def parse_bfactors_per_residue(text, fmt):
    """Per-residue pLDDT (CA B-factor) as {chain:resnum: value}."""
    bfactors = {}
    if fmt == 'pdb':
        for line in text.split('\n'):
            if not line.startswith('ATOM') or len(line) < 66:
                continue
            if line[12:16].strip() != 'CA':
                continue
            chain = line[21:22].strip() or 'A'
            try:
                rn = int(line[22:26].strip())
                bf = float(line[60:66].strip())
            except ValueError:
                continue
            bfactors[f'{chain}:{rn}'] = bf
    else:
        for row in _iter_cif_atom_site(text):
            if row.get('group_PDB') != 'ATOM' or row.get('label_atom_id') != 'CA':
                continue
            chain = row.get('label_asym_id', '')
            try:
                rn = int(row.get('label_seq_id'))
                bf = float(row.get('B_iso_or_equiv'))
            except (ValueError, TypeError):
                continue
            bfactors[f'{chain}:{rn}'] = bf

    if bfactors:
        mx = max(bfactors.values())
        if 0 < mx <= 1.0:
            bfactors = {k: v * 100.0 for k, v in bfactors.items()}
    return bfactors


def compute_chain_plddt(struct_text, fmt):
    bfs = parse_bfactors_per_residue(struct_text, fmt)
    chain_vals = {}
    for key, val in bfs.items():
        chain_vals.setdefault(key.split(':')[0], []).append(val)
    return {c: sum(v) / len(v) for c, v in chain_vals.items() if v}


# ============================================================================
# Contact map + PAE transform
# ============================================================================

def compute_contact_map(coords, threshold=8):
    """NxN binary Cb-Cb contact map (4A adjustment for phosphorus atoms)."""
    n = len(coords)
    if n == 0:
        return np.zeros((0, 0), dtype=np.uint8), 0
    xyz = np.array([[c['x'], c['y'], c['z']] for c in coords])
    has_p = np.array([c['has_p'] for c in coords])
    distances = squareform(pdist(xyz))
    p_adjustment = np.zeros_like(distances)
    p_mask = has_p[:, np.newaxis] | has_p[np.newaxis, :]
    p_adjustment[p_mask] = -4.0
    contact = ((distances + p_adjustment) < threshold).astype(np.uint8)
    return contact, n


def transform_pae_matrix(pae, pae_cutoff=12):
    """PAE -> per-direction confidence: 1 - PAE/cutoff where PAE < cutoff, else 0."""
    pae = np.asarray(pae, dtype=np.float64)
    transformed = np.zeros_like(pae)
    mask = pae < pae_cutoff
    transformed[mask] = 1.0 - pae[mask] / pae_cutoff
    return transformed


def _avg_bfactor(res_set, chain, bfactors):
    vals = [bfactors[f'{chain}:{r}'] for r in res_set if f'{chain}:{r}' in bfactors]
    return sum(vals) / len(vals) if vals else float('nan')


# ============================================================================
# Single-model analysis
# ============================================================================

def analyze_single_model(struct_text, pae_matrix, scores, fmt, platform,
                         pae_path, read_fn, pae_cutoff=12, cb_cutoff=8,
                         max_map_dim=400):
    """Compute LIS metrics + per-pair PAE/cLIA maps for one model.

    Returns dict:
        {'chains': [...], 'chain_plddt': {...}, 'pairs': [pair_dict, ...]}
    Each pair_dict carries scores plus 'pae_map', 'clia_map', 'res_i', 'res_j'
    for visualization.
    """
    chain_info = get_chains_from_structure(struct_text, fmt)

    # AF3 token_chain_ids give authoritative per-token chain assignment
    chain_ids = None
    if platform in ('alphafold3', 'openfold3') and pae_path:
        pae_content = read_fn(pae_path)
        if isinstance(pae_content, str):
            try:
                fd = json.loads(pae_content)
                if 'token_chain_ids' in fd:
                    chain_ids = fd['token_chain_ids']
            except json.JSONDecodeError:
                pass

    if chain_ids:
        chain_map = OrderedDict()
        for c in chain_ids:
            chain_map[c] = chain_map.get(c, 0) + 1
        chain_names = list(chain_map.keys())
        sizes = list(chain_map.values())
        n_total = len(chain_ids)
    else:
        chain_names = chain_info['names']
        sizes = chain_info['sizes']
        n_total = sum(sizes)

    pae = np.asarray(pae_matrix, dtype=np.float32)
    if pae.shape[0] != n_total:
        n_total = pae.shape[0]

    cum_sum = np.cumsum(sizes)
    starts = np.concatenate(([0], cum_sum[:-1]))

    transformed = transform_pae_matrix(pae[:n_total, :n_total], pae_cutoff)

    coords = parse_structure_coords(struct_text, fmt)
    contact, n_coords = compute_contact_map(coords, cb_cutoff)
    n_use = min(n_total, n_coords)

    iptm_matrix = scores.get('chainPairIptm')
    global_iptm = scores.get('ipTM', 0)
    bfs = parse_bfactors_per_residue(struct_text, fmt)
    chain_plddt = compute_chain_plddt(struct_text, fmt)

    nc = len(chain_names)
    pairs = {}

    for i in range(nc):
        for j in range(nc):
            if i == j:
                continue
            si, ei = int(starts[i]), int(min(cum_sum[i], n_total))
            sj, ej = int(starts[j]), int(min(cum_sum[j], n_total))

            t_block = transformed[si:ei, sj:ej]
            t_pos = t_block > 0
            lis_sum = float(t_block[t_pos].sum())
            lis_count_avg = int(t_pos.sum())

            t_block_rev = transformed[sj:ej, si:ei]
            either_pos = t_pos | (t_block_rev.T > 0)
            lir_i = set(np.where(either_pos.any(axis=1))[0] + 1)
            lir_j = set(np.where(either_pos.any(axis=0))[0] + 1)

            c_ei, c_ej = min(ei, n_use), min(ej, n_use)
            c_si, c_sj = si, sj
            if c_ei > c_si and c_ej > c_sj:
                contact_block = contact[c_si:c_ei, c_sj:c_ej].astype(bool)
                t_contact = t_block[:c_ei - c_si, :c_ej - c_sj]
                ct_pos = t_pos[:c_ei - c_si, :c_ej - c_sj] & contact_block
                clis_sum = float(t_contact[ct_pos].sum())
                clis_count_avg = int(ct_pos.sum())
                either_contact = either_pos[:c_ei - c_si, :c_ej - c_sj] & contact_block
                clir_i = set(np.where(either_contact.any(axis=1))[0] + 1)
                clir_j = set(np.where(either_contact.any(axis=0))[0] + 1)
            else:
                contact_block = None
                clis_sum, clis_count_avg = 0.0, 0
                clir_i, clir_j = set(), set()

            pae_ij = pae[si:ei, sj:ej]
            pae_ji = pae[sj:ej, si:ei]
            lis_count_ab = int((pae_ij < pae_cutoff).sum())
            lis_count_ba = int((pae_ji < pae_cutoff).sum())

            if contact_block is not None:
                pae_ij_c = pae_ij[:c_ei - c_si, :c_ej - c_sj]
                pae_ji_c = pae_ji[:c_ej - c_sj, :c_ei - c_si]
                clis_count_ab = int(((pae_ij_c < pae_cutoff) & contact_block).sum())
                clis_count_ba = int(((pae_ji_c < pae_cutoff) & contact_block.T).sum())
            else:
                clis_count_ab, clis_count_ba = 0, 0

            lis_val = lis_sum / lis_count_avg if lis_count_avg > 0 else 0.0
            clis_val = clis_sum / clis_count_avg if clis_count_avg > 0 else 0.0
            ilis_val = math.sqrt(lis_val * clis_val)

            iptm_val = global_iptm
            if iptm_matrix:
                try:
                    if isinstance(iptm_matrix, dict):
                        key = f'({chain_names[i]}, {chain_names[j]})'
                        if key in iptm_matrix:
                            iptm_val = float(iptm_matrix[key])
                    elif isinstance(iptm_matrix, (list, np.ndarray)):
                        if i < len(iptm_matrix) and j < len(iptm_matrix[i]):
                            v = iptm_matrix[i][j]
                            if v is not None:
                                iptm_val = float(v)
                except (KeyError, IndexError, TypeError):
                    pass

            lia_count = lis_count_ab + lis_count_ba
            clia_count = clis_count_ab + clis_count_ba

            liplddt_i = _avg_bfactor(lir_i, chain_names[i], bfs)
            liplddt_j = _avg_bfactor(lir_j, chain_names[j], bfs)
            cliplddt_i = _avg_bfactor(clir_i, chain_names[i], bfs)
            cliplddt_j = _avg_bfactor(clir_j, chain_names[j], bfs)

            pairs[f'{chain_names[i]},{chain_names[j]}'] = {
                'ci': chain_names[i], 'cj': chain_names[j],
                'LIS': lis_val, 'cLIS': clis_val, 'iLIS': ilis_val,
                'ipTM': iptm_val, 'LIA': lia_count, 'cLIA': clia_count,
                'LIpLDDT_i': liplddt_i, 'LIpLDDT_j': liplddt_j,
                'cLIpLDDT_i': cliplddt_i, 'cLIpLDDT_j': cliplddt_j,
                'lirI': lir_i, 'lirJ': lir_j, 'clirI': clir_i, 'clirJ': clir_j,
                'lenI': ei - si, 'lenJ': ej - sj,
                '_si': si, '_ei': ei, '_sj': sj, '_ej': ej,
            }

    # Symmetrize (i,j)+(j,i) into one canonical pair
    symmetric = {}
    seen = set()
    for key, val in pairs.items():
        ci, cj = key.split(',')
        canon = ','.join(sorted([ci, cj]))
        if canon in seen:
            continue
        seen.add(canon)
        rev = pairs.get(f'{cj},{ci}')
        if rev:
            s = {
                'ci': ci, 'cj': cj,
                'LIS': (val['LIS'] + rev['LIS']) / 2,
                'cLIS': (val['cLIS'] + rev['cLIS']) / 2,
                'ipTM': (val['ipTM'] + rev['ipTM']) / 2,
                'LIA': val['LIA'], 'cLIA': val['cLIA'],
                'lenI': val['lenI'], 'lenJ': val['lenJ'],
                'lirI': val['lirI'] | rev['lirJ'], 'lirJ': val['lirJ'] | rev['lirI'],
                'clirI': val['clirI'] | rev['clirJ'], 'clirJ': val['clirJ'] | rev['clirI'],
                'LIpLDDT_i': val['LIpLDDT_i'], 'LIpLDDT_j': rev['LIpLDDT_i'],
                'cLIpLDDT_i': val['cLIpLDDT_i'], 'cLIpLDDT_j': rev['cLIpLDDT_i'],
                '_si': val['_si'], '_ei': val['_ei'], '_sj': val['_sj'], '_ej': val['_ej'],
            }
        else:
            s = dict(val)
        s['iLIS'] = math.sqrt(s['LIS'] * s['cLIS'])
        s['iLIA'] = math.sqrt(s['LIA'] * s['cLIA'])
        symmetric[f'{ci}-{cj}'] = s

    # Build maps + finalize
    results = []
    for _key, v in sorted(symmetric.items()):
        si, ei, sj, ej = v['_si'], v['_ei'], v['_sj'], v['_ej']
        pae_map, clia_map, res_i, res_j = _build_pair_maps(
            pae, transformed, contact, n_use, si, ei, sj, ej, pae_cutoff, max_map_dim)
        v['pae_map'] = pae_map
        v['clia_map'] = clia_map
        v['res_i'] = res_i
        v['res_j'] = res_j
        v['pLDDT_i'] = chain_plddt.get(v['ci'])
        v['pLDDT_j'] = chain_plddt.get(v['cj'])
        v['pTM'] = scores.get('pTM')
        for k in ('_si', '_ei', '_sj', '_ej'):
            v.pop(k, None)
        results.append(v)

    return {'chains': chain_names, 'sizes': sizes, 'chain_plddt': chain_plddt, 'pairs': results}


def _downsample_indices(length, max_dim):
    """Return evenly-spaced row indices covering [0, length) capped at max_dim."""
    if length <= max_dim:
        return list(range(length))
    return [int(round(k * (length - 1) / (max_dim - 1))) for k in range(max_dim)]


def _build_pair_maps(pae, transformed, contact, n_use, si, ei, sj, ej, pae_cutoff, max_dim):
    """Inter-chain PAE block + contact-filtered LIA (cLIA) binary map.

    Both are downsampled to <= max_dim per axis for compact HTML embedding.
    Returns (pae_map, clia_map, res_i_labels, res_j_labels).
    """
    li, lj = ei - si, ej - sj
    if li <= 0 or lj <= 0:
        return [], [], [], []
    ridx = _downsample_indices(li, max_dim)
    cidx = _downsample_indices(lj, max_dim)

    pae_block = pae[si:ei, sj:ej]
    pae_map = [[round(float(pae_block[r, c]), 1) for c in cidx] for r in ridx]

    # cLIA: confident (PAE<cutoff) AND in direct contact (within n_use)
    c_ei, c_ej = min(ei, n_use), min(ej, n_use)
    clia_map = [[0 for _ in cidx] for _ in ridx]
    if c_ei > si and c_ej > sj:
        contact_block = contact[si:c_ei, sj:c_ej].astype(bool)
        conf_block = (pae_block[:c_ei - si, :c_ej - sj] < pae_cutoff)
        cmask = conf_block & contact_block
        for ri, r in enumerate(ridx):
            if r >= cmask.shape[0]:
                continue
            for cj_, c in enumerate(cidx):
                if c < cmask.shape[1] and cmask[r, c]:
                    clia_map[ri][cj_] = 1

    res_i = [r + 1 for r in ridx]
    res_j = [c + 1 for c in cidx]
    return pae_map, clia_map, res_i, res_j


# ============================================================================
# High-level entry point
# ============================================================================

def analyze_prediction(path, pae_cutoff=12, cb_cutoff=8, max_map_dim=400):
    """Analyze a single prediction (folder or file).

    Returns a list of model dicts. Each model dict:
        {'name', 'rank', 'platform', 'struct_path', 'struct_text', 'fmt',
         'chains', 'sizes', 'chain_plddt', 'pairs': [...]}
    Raises ValueError if no model with usable PAE is found.
    """
    filenames, read_fn = scan_files(path)
    if not filenames:
        raise ValueError(f"No files found at {path}")
    platform = detect_platform(filenames)

    models = []
    for name, rank, struct_path, pae_path, scores_path, fmt in find_models(filenames, platform):
        pae = extract_pae(pae_path, read_fn)
        if pae is None:
            continue
        scores = extract_confidence_scores(scores_path, read_fn)
        struct_text = read_fn(struct_path)
        if not isinstance(struct_text, str):
            continue
        analysis = analyze_single_model(
            struct_text, pae, scores, fmt, platform, pae_path, read_fn,
            pae_cutoff=pae_cutoff, cb_cutoff=cb_cutoff, max_map_dim=max_map_dim)
        models.append({
            'name': name, 'rank': rank, 'platform': platform,
            'struct_path': os.path.join(path, struct_path) if os.path.isdir(path) else path,
            'struct_text': struct_text, 'fmt': fmt,
            **analysis,
        })

    if not models:
        raise ValueError(
            f"No model with usable PAE found at {path}.\n"
            "  PAE is required for LIS and is NOT stored in the structure file.\n"
            "  AlphaFold3: need *_full_data_N.json next to *_model_N.cif\n"
            "  Boltz: need pae_*.npz / confidence_*.json\n"
            "  Chai-1: need pae.*.npy/.npz + scores.*.json")
    return models
