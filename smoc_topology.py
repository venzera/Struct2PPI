#!/usr/bin/env python3
"""
SMOC topology module -- graph-only filament-continuity check.

Call definition (also stated in results/README.md as the method definition):
    A SMOC-like filament is a single connected component whose spectral
    ordering yields a banded adjacency matrix with a flat interior
    cut-crossing profile.

Depends on the graph only -- no coordinates -- so it applies to any assembly
Struct2PPI can turn into a graph (edge list of chain pairs, optionally
weighted by contact count). Centroid geometry (PCA axial order, collinearity
R^2) is computed separately, as a baseline for comparison, never as part of
the topology call itself.

Method notes:
  - Spectral ordering: Fiedler vector (eigenvector of the second-smallest
    eigenvalue) of the SYMMETRIC NORMALIZED graph Laplacian
    L_sym = I - D^-1/2 A D^-1/2 (networkx.normalized_laplacian_matrix),
    dense eigendecomposition (all graphs here are small, N<=54). Nodes are
    ordered by increasing Fiedler-vector value. Computed per connected
    component (undefined across components).
  - Banding: for the spectral-ordered adjacency matrix, "modal offset set" =
    offsets |pos_i - pos_j| whose edge count is >= 50% of the single most
    frequent offset's count (a half-max band definition, stated explicitly
    since it is a judgement call). Only edges of >= 5 residue pairs enter
    this histogram, so a grazing chain-pair contact cannot widen the band
    (see banding_analysis). "Fraction of edges outside the modal offset set"
    is reported as the primary banding-failure diagnostic.
  - Cut-crossing profile: for each cut position k (between spectral ranks k
    and k+1, k=0..n-2), the number and summed contact-weight of edges with
    one endpoint at rank <=k and the other at rank >k. The end taper window
    is defined as the max offset in the modal offset set (the characteristic
    bandwidth) -- cuts within that many positions of either end are excluded
    from the "interior" plateau/CV/minimum calculations, since a filament's
    ends legitimately have fewer crossing edges even when continuous.
"""
import itertools
from collections import defaultdict

import networkx as nx
import numpy as np
from scipy import stats


def build_graph(edges_with_weight):
    """edges_with_weight: iterable of (chain_i, chain_j, weight)."""
    G = nx.Graph()
    for ci, cj, w in edges_with_weight:
        if G.has_edge(ci, cj):
            G[ci][cj]['weight'] += w
        else:
            G.add_edge(ci, cj, weight=w)
    return G


def fiedler_order(G_component):
    """Spectral order of a connected component's nodes by Fiedler vector of
    the symmetric normalized Laplacian. Returns (ordered_node_list, fiedler_values_dict)."""
    nodes = list(G_component.nodes())
    if len(nodes) < 3:
        return nodes, {n: 0.0 for n in nodes}
    L = nx.normalized_laplacian_matrix(G_component, nodelist=nodes, weight='weight').toarray()
    eigvals, eigvecs = np.linalg.eigh(L)
    # eigvals sorted ascending by eigh; index 0 is ~0 (connected -> single zero eigenvalue)
    fiedler = eigvecs[:, 1]
    order = np.argsort(fiedler)
    ordered_nodes = [nodes[i] for i in order]
    fiedler_values = {nodes[i]: float(fiedler[i]) for i in range(len(nodes))}
    return ordered_nodes, fiedler_values


def banding_analysis(G_component, spectral_order, min_interface_weight=5):
    """Offset histogram under spectral order; modal offset set (half-max
    rule); fraction of edges outside it.

    Only edges representing a real interface (weight >= min_interface_weight
    residue pairs) contribute to the offset histogram. Edges below that stay
    in the graph -- they still carry connectivity and still appear in the
    crossing profile -- but they do not define the band. The half-max rule
    counts edges rather than contact mass, so without this floor a few
    grazing contacts weigh as much as a full interface: in the RIPK1 fibril
    9HR6 each chain touches its stacking neighbour across 72 residue pairs
    at 2.7 Å, while the three second-neighbour pairs share 3 residue pairs
    at 3.9 Å and vanish entirely below a 3.5 Å cutoff. Those three edges
    were enough to put offset 2 in the modal set, which set bandwidth=2,
    left a 5-chain stack with no interior, and reported borderline for what
    is plainly a continuous stack.
    """
    pos = {n: i for i, n in enumerate(spectral_order)}
    offset_counts = defaultdict(int)
    offset_weight = defaultdict(float)
    total_edges = 0
    n_below_floor = 0
    for u, v, data in G_component.edges(data=True):
        off = abs(pos[u] - pos[v])
        w = data.get('weight', 1)
        if w < min_interface_weight:
            n_below_floor += 1
            continue
        offset_counts[off] += 1
        offset_weight[off] += w
        total_edges += 1
    if not offset_counts:
        return {'offset_counts': {}, 'modal_offsets': set(), 'frac_outside_modal': 1.0,
                'bandwidth': 0, 'n_edges_below_interface_floor': n_below_floor,
                'min_interface_weight': min_interface_weight}
    max_count = max(offset_counts.values())
    modal_offsets = {off for off, c in offset_counts.items() if c >= 0.5 * max_count}
    n_outside = sum(c for off, c in offset_counts.items() if off not in modal_offsets)
    frac_outside = n_outside / total_edges if total_edges else 0.0
    bandwidth = max(modal_offsets) if modal_offsets else 0
    return {'offset_counts': dict(offset_counts), 'offset_weight': dict(offset_weight),
            'modal_offsets': modal_offsets, 'frac_outside_modal': frac_outside,
            'bandwidth': bandwidth, 'total_edges': total_edges,
            'n_edges_below_interface_floor': n_below_floor,
            'min_interface_weight': min_interface_weight}


def crossing_profile(G_component, spectral_order):
    """For each cut k=0..n-2: (n_edges_crossing, weight_crossing)."""
    pos = {n: i for i, n in enumerate(spectral_order)}
    n = len(spectral_order)
    edge_positions = [(pos[u], pos[v], data.get('weight', 1)) for u, v, data in G_component.edges(data=True)]
    counts = np.zeros(n - 1, dtype=int)
    weights = np.zeros(n - 1, dtype=float)
    for pu, pv, w in edge_positions:
        lo, hi = min(pu, pv), max(pu, pv)
        counts[lo:hi] += 1
        weights[lo:hi] += w
    return counts, weights


def interior_stats(counts, weights, end_window):
    n_cuts = len(counts)
    lo = end_window
    hi = n_cuts - end_window
    # degenerate: structure too short for its own bandwidth to leave an
    # interior distinct from the end taper. Falling back to the full range
    # here would mix in the legitimately-lower end-taper cuts and inflate
    # CV for a reason that has nothing to do with interruption -- flag this
    # explicitly rather than silently computing a misleading CV over it.
    #
    # A single interior cut counts as degenerate too: the CV of one value is
    # zero by construction and its min equals its plateau, so both gates in
    # classify_topology pass without testing anything. That is how the PurE
    # tetramer 3RGG, a closed 4-chain clique, read as a continuous filament.
    # The reference filaments 3J63 and 2N1F have two interior cuts, so two
    # is the smallest width that can carry evidence here.
    degenerate = (hi - lo) < 2
    if degenerate:
        lo, hi = 0, n_cuts
    interior_counts = counts[lo:hi]
    interior_weights = weights[lo:hi]
    if len(interior_weights) == 0:
        return {'plateau_weight': None, 'cv_weight': None, 'min_interior_weight': None,
                'min_interior_pos': None, 'end_window_used': end_window, 'interior_range': (lo, hi),
                'degenerate_interior': True}
    plateau = float(np.mean(interior_weights))
    cv = float(np.std(interior_weights) / plateau) if plateau > 0 else float('inf')
    min_idx_local = int(np.argmin(interior_weights))
    min_pos = lo + min_idx_local
    return {'plateau_weight': plateau, 'cv_weight': cv,
            'min_interior_weight': float(interior_weights[min_idx_local]),
            'min_interior_count': int(interior_counts[min_idx_local]),
            'min_interior_pos': min_pos, 'end_window_used': end_window, 'interior_range': (lo, hi),
            'degenerate_interior': degenerate}


def degree_distribution(G_component, spectral_order, end_window):
    n = len(spectral_order)
    lo, hi = end_window, n - end_window
    if hi <= lo:
        lo, hi = 0, n
    interior_nodes = spectral_order[lo:hi]
    degs = [G_component.degree(node, weight=None) for node in interior_nodes]
    return {'interior_nodes': interior_nodes, 'degrees': degs,
            'mean_degree': float(np.mean(degs)) if degs else None}


def order_correlation(spectral_order, reference_order_map):
    """Spearman correlation between spectral rank and a reference rank
    (e.g. AF3 chain-letter index, or PCA axial projection rank), for nodes
    present in both."""
    common = [n for n in spectral_order if n in reference_order_map]
    if len(common) < 3:
        return None
    spectral_rank = {n: i for i, n in enumerate(spectral_order)}
    xs = [spectral_rank[n] for n in common]
    ys = [reference_order_map[n] for n in common]
    if len(set(ys)) < 2:
        return None
    rho = stats.spearmanr(xs, ys).correlation
    return float(rho) if rho == rho else None


def analyze_component(G_component, chain_index_map=None, pca_rank_map=None):
    spectral_order, fiedler_values = fiedler_order(G_component)
    band = banding_analysis(G_component, spectral_order)
    counts, weights = crossing_profile(G_component, spectral_order)
    end_window = max(1, band['bandwidth'])
    istats = interior_stats(counts, weights, end_window)
    degdist = degree_distribution(G_component, spectral_order, end_window)
    diameter = nx.diameter(G_component) if nx.is_connected(G_component) else None
    mean_path = nx.average_shortest_path_length(G_component) if nx.is_connected(G_component) else None

    corr_chain_index = order_correlation(spectral_order, chain_index_map) if chain_index_map else None
    corr_pca = order_correlation(spectral_order, pca_rank_map) if pca_rank_map else None

    return {
        'nodes': list(G_component.nodes()), 'n_nodes': G_component.number_of_nodes(),
        'n_edges': G_component.number_of_edges(),
        'spectral_order': spectral_order, 'fiedler_values': fiedler_values,
        'banding': band, 'crossing_counts': counts.tolist(), 'crossing_weights': weights.tolist(),
        'interior_stats': istats, 'degree_distribution': degdist,
        'diameter': diameter, 'mean_shortest_path': mean_path,
        'spearman_vs_chain_index': corr_chain_index, 'spearman_vs_pca_order': corr_pca,
    }


def classify_topology(components_analysis, n_total_chains,
                       banding_frac_outside_thresh=0.3,
                       dip_ratio_thresh=0.5,
                       cv_thresh=0.15):
    """Category from a list of per-component analyze_component() outputs.

    REVISED (after validating on a 64-structure death-fold PDB survey
    spanning DED/PYD/DD/RHIM families): banding (frac_outside_modal) is
    reported but no longer a hard gate on the call. Evidence: comparing two
    same-family, same-size structures (9U6E vs 9U7A, both N=24, 60 edges)
    showed near-identical, cleanly FLAT crossing profiles (CV 0.082 vs
    0.004) but opposite banding fractions (0.567 vs 0.233) -- banding is
    sensitive to a filament's specific helical parameters (start-number,
    rise/rotation) and a threshold fit to the ASC-PYD 3-start system does
    not transfer to other families. The crossing-profile SHAPE is the more
    architecture-agnostic signal and still discriminates correctly (e.g.
    correctly flat for the genuine MyD88 DD filament 6I3N, correctly
    non-flat/ramping for the known false case, human ASC-PYD N=42).

    Primary gate is now the interior crossing-weight coefficient of
    variation (CV): cv_thresh=0.15 sits between the highest CV among cases
    that should read as continuous (6I3N=0.119, a real filament) and the
    lowest CV among cases that should not (human N=42=0.227, two filaments
    lying alongside each other -- its dip_ratio alone does NOT catch this,
    since the failure mode is a broad asymmetric ramp, not a sharp
    localized dip, which is why CV rather than dip_ratio is now primary).
    dip_ratio remains a SEPARATE check for a sharp localized dip (a genuine
    interruption at one position) even within an otherwise low-CV profile.

    Categories: continuous filament / interrupted filament /
    multiple filaments / non-filamentous / borderline (with reason).
    """
    large_components = [c for c in components_analysis if c['n_nodes'] >= 3]
    if len(large_components) > 1:
        return {'call': 'multiple filaments', 'reason': f'{len(large_components)} components with >=3 chains',
                'n_components': len(large_components)}

    if not large_components:
        return {'call': 'non-filamentous', 'reason': 'no component with >=3 chains', 'n_components': 0}

    comp = large_components[0]
    frac_outside = comp['banding']['frac_outside_modal']

    istats = comp['interior_stats']
    plateau = istats['plateau_weight']
    min_w = istats['min_interior_weight']
    cv = istats['cv_weight']
    if plateau is None or plateau == 0:
        return {'call': 'borderline', 'reason': 'no defined interior (too short / bandwidth too large)',
                'n_components': 1}

    dip_ratio = min_w / plateau

    # too short for its own bandwidth to have a true interior distinct from
    # the end taper -- CV over the full (taper-included) profile isn't a
    # meaningful flatness measure here. Report borderline rather than guess.
    if istats.get('degenerate_interior'):
        return {'call': 'borderline',
                'reason': f'fewer than 2 interior cuts distinct from end-taper (structure too short for '
                          f'its own bandwidth={comp["banding"]["bandwidth"]}); CV={cv:.3f} not reliable '
                          f'here (banding frac_outside_modal={frac_outside:.2f}, reported not gating)',
                'n_components': 1}

    # sharp localized dip -> interrupted, regardless of overall CV
    if dip_ratio < dip_ratio_thresh:
        return {'call': 'interrupted filament',
                'reason': f'interior crossing dip: min/plateau = {dip_ratio:.2f} < {dip_ratio_thresh} '
                          f'(banding frac_outside_modal={frac_outside:.2f}, reported not gating)',
                'dip_position': istats['min_interior_pos'],
                'flanking_chains': (comp['spectral_order'][istats['min_interior_pos']],
                                     comp['spectral_order'][istats['min_interior_pos'] + 1])
                                    if istats['min_interior_pos'] is not None else None,
                'n_components': 1}

    # broad non-flatness (ramp/scatter rather than a sharp dip) -> non-filamentous,
    # with a borderline band just above the threshold checked FIRST so a value
    # like 0.158 against a 0.15 threshold reads as borderline, not a hard fail.
    if cv > cv_thresh * 1.2:
        return {'call': 'non-filamentous',
                'reason': f'interior crossing profile not flat (CV={cv:.3f} > {cv_thresh*1.2:.3f}; '
                          f'banding frac_outside_modal={frac_outside:.2f}, reported not gating)',
                'n_components': 1}

    if cv > cv_thresh:
        return {'call': 'borderline',
                'reason': f'CV={cv:.3f} in the borderline band above threshold {cv_thresh} '
                          f'(banding frac_outside_modal={frac_outside:.2f}, reported not gating)',
                'n_components': 1}

    # borderline zone just above dip_ratio threshold
    if dip_ratio < dip_ratio_thresh * 1.3:
        return {'call': 'borderline', 'reason': f'dip_ratio={dip_ratio:.2f} close to threshold {dip_ratio_thresh}',
                'n_components': 1}

    return {'call': 'continuous filament',
            'reason': f'dip_ratio={dip_ratio:.2f}, CV={cv:.3f}, frac_outside_modal={frac_outside:.2f} (not gating)',
            'n_components': 1}


def full_topology_report(edges_with_weight, chain_index_map=None, pca_rank_map=None):
    """Top-level entry point: edges_with_weight -> full report dict."""
    G = build_graph(edges_with_weight)
    n_total = G.number_of_nodes()
    components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    components_analysis = [analyze_component(c, chain_index_map, pca_rank_map) for c in components]
    call = classify_topology(components_analysis, n_total)
    return {'n_total_chains': n_total, 'n_components': len(components),
            'component_sizes': sorted([c.number_of_nodes() for c in components], reverse=True),
            'components': components_analysis, 'call': call}
