#!/usr/bin/env python3
"""
Generate correlation plots comparing FoldX and PRODIGY binding energies.
"""

import re
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from scipy import stats

def parse_foldx(filepath):
    """Parse FoldX output file and return dict of chain_pair -> dG."""
    foldx_data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('./') or line.startswith('8xks'):
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chain1 = parts[1]
                    chain2 = parts[2]
                    try:
                        dg = float(parts[5])  # Interaction Energy column
                        pair = tuple(sorted([chain1, chain2]))
                        # Only keep if there's actual interaction (non-zero or significant)
                        if abs(dg) > 0.1:
                            foldx_data[pair] = dg
                    except (ValueError, IndexError):
                        continue
    return foldx_data

def parse_prodigy(filepath):
    """Parse PRODIGY binding strength file and return dict of chain_pair -> dG."""
    prodigy_data = {}
    with open(filepath, 'r') as f:
        for line in f:
            # Match lines like "1     E-H            -43.90"
            match = re.match(r'\d+\s+([A-Z])-([A-Z])\s+([-\d.]+)', line.strip())
            if match:
                chain1, chain2, dg = match.groups()
                pair = tuple(sorted([chain1, chain2]))
                prodigy_data[pair] = float(dg)
    return prodigy_data

def calculate_correlation(x, y):
    """Calculate Pearson correlation coefficient and p-value."""
    r, p = stats.pearsonr(x, y)
    return r, p

# Parse all data files
foldx_data = parse_foldx('Interaction_8xks_Repair_AC.fxout')
prodigy_raw = parse_prodigy('../test/8xks_binding_strength.txt')
prodigy_repaired = parse_prodigy('8xks_Repair_binding_strength.txt')

# Find intersecting pairs for each comparison
pairs_raw = set(foldx_data.keys()) & set(prodigy_raw.keys())
pairs_repaired = set(foldx_data.keys()) & set(prodigy_repaired.keys())

print(f"FoldX pairs: {len(foldx_data)}")
print(f"PRODIGY Raw pairs: {len(prodigy_raw)}")
print(f"PRODIGY Repaired pairs: {len(prodigy_repaired)}")
print(f"Intersecting pairs (Raw): {len(pairs_raw)}")
print(f"Intersecting pairs (Repaired): {len(pairs_repaired)}")

# Prepare data for plotting
foldx_raw_values = [foldx_data[p] for p in pairs_raw]
prodigy_raw_values = [prodigy_raw[p] for p in pairs_raw]
labels_raw = [f"{p[0]}-{p[1]}" for p in pairs_raw]

foldx_repaired_values = [foldx_data[p] for p in pairs_repaired]
prodigy_repaired_values = [prodigy_repaired[p] for p in pairs_repaired]
labels_repaired = [f"{p[0]}-{p[1]}" for p in pairs_repaired]

# Calculate correlations
r_raw, p_raw = calculate_correlation(foldx_raw_values, prodigy_raw_values)
r_repaired, p_repaired = calculate_correlation(foldx_repaired_values, prodigy_repaired_values)

print(f"\nCorrelation (FoldX vs PRODIGY Raw): r={r_raw:.3f}, p={p_raw:.2e}")
print(f"Correlation (FoldX vs PRODIGY Repaired): r={r_repaired:.3f}, p={p_repaired:.2e}")

# Create subplots
fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=(
        f'FoldX vs PRODIGY (Raw)<br>r={r_raw:.3f}, n={len(pairs_raw)}',
        f'FoldX vs PRODIGY (Repaired)<br>r={r_repaired:.3f}, n={len(pairs_repaired)}'
    ),
    horizontal_spacing=0.12
)

# Subplot 1: FoldX vs PRODIGY Raw
fig.add_trace(
    go.Scatter(
        x=foldx_raw_values,
        y=prodigy_raw_values,
        mode='markers',
        marker=dict(size=8, color='#1f77b4', opacity=0.7),
        text=labels_raw,
        hovertemplate='<b>%{text}</b><br>FoldX: %{x:.1f}<br>PRODIGY: %{y:.1f}<extra></extra>',
        name='Raw'
    ),
    row=1, col=1
)

# Add trendline for subplot 1
z_raw = np.polyfit(foldx_raw_values, prodigy_raw_values, 1)
p_raw_line = np.poly1d(z_raw)
x_range_raw = np.linspace(min(foldx_raw_values), max(foldx_raw_values), 100)
fig.add_trace(
    go.Scatter(
        x=x_range_raw,
        y=p_raw_line(x_range_raw),
        mode='lines',
        line=dict(color='red', dash='dash'),
        name='Trendline',
        showlegend=False
    ),
    row=1, col=1
)

# Subplot 2: FoldX vs PRODIGY Repaired
fig.add_trace(
    go.Scatter(
        x=foldx_repaired_values,
        y=prodigy_repaired_values,
        mode='markers',
        marker=dict(size=8, color='#2ca02c', opacity=0.7),
        text=labels_repaired,
        hovertemplate='<b>%{text}</b><br>FoldX: %{x:.1f}<br>PRODIGY: %{y:.1f}<extra></extra>',
        name='Repaired'
    ),
    row=1, col=2
)

# Add trendline for subplot 2
z_repaired = np.polyfit(foldx_repaired_values, prodigy_repaired_values, 1)
p_repaired_line = np.poly1d(z_repaired)
x_range_repaired = np.linspace(min(foldx_repaired_values), max(foldx_repaired_values), 100)
fig.add_trace(
    go.Scatter(
        x=x_range_repaired,
        y=p_repaired_line(x_range_repaired),
        mode='lines',
        line=dict(color='red', dash='dash'),
        name='Trendline',
        showlegend=False
    ),
    row=1, col=2
)

# Update layout
fig.update_layout(
    title=dict(
        text='Binding Energy Correlation: FoldX vs PRODIGY (8XKS)',
        x=0.5,
        font=dict(size=16)
    ),
    showlegend=False,
    width=1000,
    height=500,
    template='plotly_white'
)

# Update axes
fig.update_xaxes(title_text='FoldX ΔG (kcal/mol)', row=1, col=1)
fig.update_yaxes(title_text='PRODIGY ΔG (kcal/mol)', row=1, col=1)
fig.update_xaxes(title_text='FoldX ΔG (kcal/mol)', row=1, col=2)
fig.update_yaxes(title_text='PRODIGY ΔG (kcal/mol)', row=1, col=2)

# Save as PNG
fig.write_image('evaluation_8xks.png', scale=2)
print(f"\nSaved: evaluation_8xks.png")

# Also save as HTML for interactive viewing
fig.write_html('evaluation_8xks.html')
print(f"Saved: evaluation_8xks.html")
