# -*- coding: utf-8 -*-
"""
Created on Thu Dec 25 20:46:07 2025

@author: Arkin Lab
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


f_ph = Path(r"F:\branchwater\phac_scan\Betaproteobacteria\betaproteobacteria_files_mapped_to_bestlist_rank.csv")
ph = pd.read_csv(f_ph)

r = pd.to_numeric(ph["rank_in_bestlist"], errors="coerce").dropna().astype(int)

n_unique_le800 = (r.drop_duplicates() <= 800).sum()
n_all_le800    = (r <= 800).sum()

print("PhaC files ranks <= 800:")
print("  unique ranks:", n_unique_le800)
print("  all rows    :", n_all_le800)
print("  total matched rows:", len(r), " total unique ranks:", r.nunique())


base = Path(r"F:\branchwater\phac_scan\Betaproteobacteria")

f_no = base / "04_NO_PhaC_TOP50_with_rank_in_bestlist.csv"
f_ph = base / "betaproteobacteria_files_mapped_to_bestlist_rank.csv"

out_png = base / "rank_positions_noPhaC_TOP50_vs_PhaC_files.png"
out_pdf = base / "rank_positions_noPhaC_TOP50_vs_PhaC_files.pdf"

no = pd.read_csv(f_no) 
ph = pd.read_csv(f_ph)

x_no = pd.to_numeric(no["rank_in_bestlist"], errors="coerce").dropna().astype(int)
x_ph = pd.to_numeric(ph["rank_in_bestlist"], errors="coerce").dropna().astype(int)

x_no_u = sorted(set(x_no.tolist()))
x_ph_u = sorted(set(x_ph.tolist()))

y_no = [1] * len(x_no_u)
y_ph = [0] * len(x_ph_u)


fig, ax = plt.subplots(figsize=(12, 3.5))

ax.scatter(x_ph_u, y_ph, s=18, alpha=0.8, label=f"PhaC TOP50 (unique ranks={len(x_ph_u)})")
ax.scatter(x_no_u, y_no, s=18, alpha=0.8, label=f"NO PhaC TOP50 (unique ranks={len(x_no_u)})")

# cutoff = 第38个橙点（unique rank）的位置
if len(x_ph_u) < 38:
    raise ValueError(f"橙点数量不足 38 个：只有 {len(x_ph_u)} 个 unique ranks")
cutoff = x_ph_u[37]

ax.axvline(cutoff, linestyle="--", linewidth=1.5, label=f"cutoff = {cutoff} (38th)")

ax.set_yticks([0, 1])
ax.set_yticklabels(["PhaC files", "NO-PhaC TOP50"])
ax.set_xlabel("Rank in bestlist (1 = highest quality/priority)")
ax.set_ylim(-0.6, 1.6)
ax.grid(True, axis="x", alpha=0.25)
ax.set_title("Rank positions in bestlist")

ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)

fig.tight_layout()
fig.savefig(out_png, dpi=200, bbox_inches="tight")
fig.savefig(out_pdf, bbox_inches="tight")
plt.show()

