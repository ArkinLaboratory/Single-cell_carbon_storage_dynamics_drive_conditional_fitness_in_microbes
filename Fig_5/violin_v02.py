# -*- coding: utf-8 -*-
"""
Violin (blue density) + boxplot overlay for depth_proxy = -ln(1-containment)
Drop containment==0 records. Compare two datasets in two environments.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =======================
# INPUT CSVs (your paths)
# =======================
CSV_1 = Path(r"F:\branchwater\ncbi_cupria\no phac compare\phac_betaproteo_random50_merged_thres_0.csv")
CSV_2 = Path(r"F:\branchwater\ncbi_cupria\no phac compare\TOP50_GENOME_SINGLE_RECORD_thres_0.csv")

LABEL_1 = "PhaC_random50"
LABEL_2 = "TOP50_single_record"

ENVS = ["human gut metagenome", "soil metagenome", "activated sludge metagenome"]

# =======================
# Helpers
# =======================
def load_and_prepare(path: Path, dataset_label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    df = pd.read_csv(path)

    required = {"acc", "organism", "containment"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} missing required columns: {sorted(missing)}")

    df["containment"] = pd.to_numeric(df["containment"], errors="coerce")
    df = df.dropna(subset=["containment"])

    # drop containment==0
    df = df[df["containment"] > 0]

    # avoid inf if containment==1
    df["containment"] = df["containment"].clip(lower=0, upper=1 - 1e-12)

    # depth proxy
    df["depth_proxy"] = -np.log1p(-df["containment"])
    df["dataset"] = dataset_label
    return df

def per_sample_depth(df: pd.DataFrame, envs: list[str]) -> pd.DataFrame:
    sub = df[df["organism"].isin(envs)].copy()
    out = (
        sub.groupby(["dataset", "organism", "acc"], as_index=False)
           .agg(depth_sum=("depth_proxy", "sum"),
                depth_mean=("depth_proxy", "mean"),
                n_hits=("depth_proxy", "size"))
    )
    return out

# =======================
# Main
# =======================
df1 = load_and_prepare(CSV_1, LABEL_1)
df2 = load_and_prepare(CSV_2, LABEL_2)

agg = pd.concat(
    [per_sample_depth(df1, ENVS), per_sample_depth(df2, ENVS)],
    ignore_index=True
)

# optional summary
summary = (agg.groupby(["organism", "dataset"])
             .agg(n_samples=("acc", "nunique"),
                  median_depth=("depth_sum", "median"),
                  mean_depth=("depth_sum", "mean"))
             .reset_index())
print("\n=== Summary (per environment x dataset) ===")
print(summary.to_string(index=False))

# =======================
# Output paths
# =======================
out_dir = CSV_1.parent / "depth_compare_out"
out_dir.mkdir(parents=True, exist_ok=True)

csv_out = out_dir / "per_sample_depth_proxy_humangut_soil.csv"
png_out = out_dir / "violin_box_depth_proxy_humangut_soil.png"
agg.to_csv(csv_out, index=False)
print(f"\nWrote table: {csv_out}")

# =======================
# Plot helper: grouped envs in ONE axis (scheme A style)
# =======================
import matplotlib.patches as mpatches

datasets_order = [LABEL_1, LABEL_2]
dataset_style = {
    LABEL_1: dict(color="tab:blue", alpha=0.25),
    LABEL_2: dict(color="tab:orange", alpha=0.25),
}

def plot_env_group(env_list, ylim, out_prefix):
    fig, ax = plt.subplots(figsize=(8.0, 4.2), constrained_layout=True)

    centers = np.arange(1, len(env_list) + 1)

    offset = 0.20
    width_violin = 0.32
    width_box = 0.10

    all_positions, all_data, all_colors, all_alphas = [], [], [], []

    for i, env in enumerate(env_list):
        sub_env = agg[agg["organism"] == env]
        for j, ds in enumerate(datasets_order):
            vals = sub_env[sub_env["dataset"] == ds]["depth_sum"].dropna().to_numpy()
            if len(vals) == 0:
                continue

            pos = centers[i] + (-offset if j == 0 else offset)
            all_positions.append(pos)
            all_data.append(vals)

            st = dataset_style.get(ds, dict(color="tab:blue", alpha=0.25))
            all_colors.append(st["color"])
            all_alphas.append(st["alpha"])

    if len(all_data) == 0:
        print(f"[WARN] No data for env_list={env_list}, skip plotting.")
        plt.close(fig)
        return

    # violin
    vp = ax.violinplot(
        all_data,
        positions=all_positions,
        widths=width_violin,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )
    for body, c, a in zip(vp["bodies"], all_colors, all_alphas):
        body.set_facecolor(c)
        body.set_edgecolor("none")
        body.set_alpha(a)

    # box overlay
    bp = ax.boxplot(
        all_data,
        positions=all_positions,
        widths=width_box,
        showfliers=False,
        patch_artist=True
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("none")
    for element in ["whiskers", "caps"]:
        for line in bp[element]:
            line.set_color("black")
    
    for line in bp["medians"]:
        line.set_color("tab:red")     # 你想要的颜色
        line.set_linewidth(1.5)   

    # x ticks at env centers
    ax.set_xticks(centers)
    ax.set_xticklabels(env_list, rotation=0)

    ax.set_ylabel("depth proxy = -ln(1 - containment)\n(per sample; summed over hits)")
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(*ylim)

    # legend
    handles = [
        mpatches.Patch(facecolor=dataset_style[LABEL_1]["color"], alpha=dataset_style[LABEL_1]["alpha"], label=LABEL_1),
        mpatches.Patch(facecolor=dataset_style[LABEL_2]["color"], alpha=dataset_style[LABEL_2]["alpha"], label=LABEL_2),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right")

    pdf_out = out_dir / f"{out_prefix}.pdf"
    png_out = out_dir / f"{out_prefix}.png"
    fig.savefig(pdf_out)
    fig.savefig(png_out, dpi=300)
    print(f"Wrote PDF: {pdf_out}")
    print(f"Wrote PNG: {png_out}")

    plt.show()


# =======================
# Make TWO figures:
# 1) human gut only, zoomed
# 2) soil + activated sludge together, broader scale
# =======================
plot_env_group(
    env_list=["human gut metagenome"],
    ylim=(0, 0.06),
    out_prefix="violin_box_depth_proxy_humangut_zoom_0_0p06"
)

plot_env_group(
    env_list=["soil metagenome", "activated sludge metagenome"],
    ylim=(0, 0.4),
    out_prefix="violin_box_depth_proxy_soil_sludge_0_0p4"
)

