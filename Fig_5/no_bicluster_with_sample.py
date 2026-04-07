# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:51:15 2026

@author: Arkin Lab
"""

# -*- coding: utf-8 -*-
"""
2D Hierarchical clustering on env×species median-per-run PHB matrix,
then:
  1) randomly sample 20×20 (env×species) to plot a sampled heatmap;
  2) plot the full heatmap.
All figures saved to .../exports/heatmap.

Requires:
    pip install scipy scikit-learn
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, optimal_leaf_ordering
from scipy.spatial.distance import pdist

# ============== CONFIG ==============
RESULTS_DIR  = r"F:\branchwater\results\local"
OUT_DIR      = os.path.join(RESULTS_DIR, "exports")
HEATMAP_DIR  = os.path.join(OUT_DIR, "heatmap")   # 新的输出文件夹
METADATA_DIR = r"F:\branchwater\metadata"         # 这里这次不会写入任何 CSV

# Must match thresholds/TAG you used before
CANI_THR = 0.93
CONT_THR = 0.20

# Hierarchical clustering params
ROW_METHOD  = 'ward'
COL_METHOD  = 'ward'
ROW_METRIC  = 'euclidean'
COL_METRIC  = 'euclidean'

N_ENV_CLUSTERS     = 2
N_SPECIES_CLUSTERS = 2

LOG_TRANSFORM = True

# 抽样参数
RANDOM_SEED     = 0
SAMPLED_N_ENV   = 20  # y 轴抽样个数
SAMPLED_N_SPEC  = 20  # x 轴抽样个数

# “完整版”热图显示标签的最大数量（抽样显示刻度，避免挤）
MAX_X_LABELS_FULL   = 60
MAX_Y_LABELS_FULL   = 50
LABEL_FONTSIZE_FULL = 9
ROTATE_X_FULL       = 60

# “抽样版”热图字体（可大些）
LABEL_FONTSIZE_SAMPLED = 11
ROTATE_X_SAMPLED       = 45

def threshold_tag(cani: float, cont: float) -> str:
    return f"cani{cani:.2f}_cont{cont:.2f}".replace(".", "p")

TAG = threshold_tag(CANI_THR, CONT_THR)

matrix_csv = os.path.join(
    OUT_DIR, f"env_species_median_per_run_matrix_{TAG}.csv"
)

if not os.path.exists(matrix_csv):
    raise FileNotFoundError(
        f"Missing {matrix_csv} – run the median-per-run heatmap script first."
    )

os.makedirs(HEATMAP_DIR, exist_ok=True)

# ============== 1) Load matrix ==============
mat = pd.read_csv(matrix_csv, index_col=0)   # rows = env, cols = species
env_names = mat.index.to_list()
species_names = mat.columns.to_list()
X = mat.values.astype(float)

if LOG_TRANSFORM:
    X_work = np.log10(X + 1.0)
else:
    X_work = X

# ============== 2) Hierarchical clustering ==============
# rows
row_dist = pdist(X_work, metric=ROW_METRIC)
row_link = linkage(row_dist, method=ROW_METHOD, optimal_ordering=False)
row_link = optimal_leaf_ordering(row_link, row_dist)
row_order = dendrogram(row_link, no_plot=True)['leaves']
row_labels = fcluster(row_link, t=N_ENV_CLUSTERS, criterion='maxclust')

# cols
col_dist = pdist(X_work.T, metric=COL_METRIC)
col_link = linkage(col_dist, method=COL_METHOD, optimal_ordering=False)
col_link = optimal_leaf_ordering(col_link, col_dist)
col_order = dendrogram(col_link, no_plot=True)['leaves']
col_labels = fcluster(col_link, t=N_SPECIES_CLUSTERS, criterion='maxclust')

# ============== 3) Reorder matrix by leaves ==============
X_reordered = X[np.ix_(row_order, col_order)]
env_reordered = [env_names[i] for i in row_order]
species_reordered = [species_names[j] for j in col_order]

# 供绘图的值（颜色）
X_plot_full = np.log10(X_reordered + 1.0) if LOG_TRANSFORM else X_reordered

# ============== 4) Sampled heatmap (20×20) ==============
np.random.seed(RANDOM_SEED)

ny, nx = X_reordered.shape
n_env_pick  = min(SAMPLED_N_ENV, ny)
n_spec_pick = min(SAMPLED_N_SPEC, nx)

# 在“重排后的顺序”上随机抽样，保证可读+与聚类一致
sampled_row_idx = np.sort(np.random.choice(ny, size=n_env_pick, replace=False))
sampled_col_idx = np.sort(np.random.choice(nx, size=n_spec_pick, replace=False))

X_sampled = X_reordered[np.ix_(sampled_row_idx, sampled_col_idx)]
env_sampled = [env_reordered[i] for i in sampled_row_idx]
spec_sampled = [species_reordered[j] for j in sampled_col_idx]
X_plot_sampled = np.log10(X_sampled + 1.0) if LOG_TRANSFORM else X_sampled

# 画“抽样版”热图
fig_h = 0.42 * n_env_pick + 2.5
fig_w = 0.32 * n_spec_pick + 4.5
fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=200, constrained_layout=True)

im = ax.imshow(X_plot_sampled, aspect="auto")
ax.set_yticks(np.arange(n_env_pick))
ax.set_yticklabels(env_sampled, fontsize=LABEL_FONTSIZE_SAMPLED)
ax.set_xticks(np.arange(n_spec_pick))
ax.set_xticklabels(spec_sampled, rotation=ROTATE_X_SAMPLED, ha="right",
                   fontsize=LABEL_FONTSIZE_SAMPLED)

ax.tick_params(axis='both', length=0)
ax.set_xlabel("Species (genome IDs)")
ax.set_ylabel("Environment")
ax.set_title(
    "Sampled 20×20 heatmap after 2D hierarchical ordering\n"
    f"{TAG}, rows: {ROW_METHOD}/{ROW_METRIC}, cols: {COL_METHOD}/{COL_METRIC}"
)

cbar = fig.colorbar(im, ax=ax, shrink=0.9)
cbar.set_label("log10(median per-run PHB mass proxy + 1)" if LOG_TRANSFORM
               else "median per-run PHB mass proxy")

out_base_sampled = os.path.join(HEATMAP_DIR, f"hier_heatmap_sampled20x20_{TAG}")
fig.savefig(out_base_sampled + ".png", dpi=300, bbox_inches="tight")
fig.savefig(out_base_sampled + ".svg", bbox_inches="tight")
fig.savefig(out_base_sampled + ".pdf", bbox_inches="tight")
plt.show()
print("[Saved sampled figure]", out_base_sampled + ".png")

# ============== 5) Full heatmap ==============
def downsample_idx(n, max_labels):
    step = max(1, int(np.ceil(n / max_labels)))
    return np.arange(0, n, step)

x_idx_full = downsample_idx(nx, MAX_X_LABELS_FULL)
y_idx_full = downsample_idx(ny, MAX_Y_LABELS_FULL)

fig_h_full = 0.36 * len(y_idx_full) + 2.5
fig_w_full = 0.28 * len(x_idx_full) + 4.5
fig, ax = plt.subplots(figsize=(fig_w_full, fig_h_full), dpi=200, constrained_layout=True)

im = ax.imshow(X_plot_full, aspect="auto")

ax.set_yticks(y_idx_full)
ax.set_yticklabels([env_reordered[i] for i in y_idx_full], fontsize=LABEL_FONTSIZE_FULL)
ax.set_xticks(x_idx_full)
ax.set_xticklabels([species_reordered[i] for i in x_idx_full], rotation=ROTATE_X_FULL,
                   ha="right", fontsize=LABEL_FONTSIZE_FULL)

ax.tick_params(axis='both', length=0)
ax.set_xlabel("Species (genome IDs)")
ax.set_ylabel("Environment")
ax.set_title(
    "Full heatmap after 2D hierarchical ordering\n"
    f"{TAG}, rows: {ROW_METHOD}/{ROW_METRIC}, cols: {COL_METHOD}/{COL_METRIC}"
)

cbar = fig.colorbar(im, ax=ax, shrink=0.9)
cbar.set_label("log10(median per-run PHB mass proxy + 1)" if LOG_TRANSFORM
               else "median per-run PHB mass proxy")

out_base_full = os.path.join(HEATMAP_DIR, f"hier_heatmap_full_{TAG}")
fig.savefig(out_base_full + ".png", dpi=300, bbox_inches="tight")
fig.savefig(out_base_full + ".svg", bbox_inches="tight")
fig.savefig(out_base_full + ".pdf", bbox_inches="tight")
plt.show()
print("[Saved full figure]", out_base_full + ".png")