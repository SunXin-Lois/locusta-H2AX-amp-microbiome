"""
Table S4 科水平相对丰度堆叠条形图
- 使用表中全部 7 类：5 个优势科 + Unclassified + Others
- 按组均值堆叠：PTW / 1 dpi / 3 dpi / 5 dpi
- 显著性：Mann-Whitney U vs PTW，Bonferroni 校正后标星
"""

from pathlib import Path
import warnings

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

HERE = Path(__file__).resolve().parent
INPUT_FILE = HERE / "Table_S4_family_relative_abundance.tsv"
OUT_PREFIX = HERE / "Family_abundance_S4_Unclassified_Others"

GROUPS = {
    "PTW": ["PTW-1", "PTW-2", "PTW-3", "PTW-4", "PTW-5"],
    "1 dpi": ["1 dpi-1", "1 dpi-2", "1 dpi-3", "1 dpi-4", "1 dpi-5"],
    "3 dpi": ["3 dpi-1", "3 dpi-2", "3 dpi-3", "3 dpi-4", "3 dpi-5"],
    "5 dpi": ["5 dpi-1", "5 dpi-2", "5 dpi-3", "5 dpi-4", "5 dpi-5"],
}
X_TICKS = ["PTW", "1", "3", "5"]
GROUP_LABELS = list(GROUPS.keys())
CONTROL_GROUP = "PTW"
TREAT_GROUPS = [g for g in GROUP_LABELS if g != CONTROL_GROUP]

# 与参考图一致：从下到上
TAXON_ORDER = [
    "Enterobacteriaceae",
    "Streptococcaceae",
    "Leuconostocaceae",
    "Enterococcaceae",
    "Yersiniaceae",
    "Unclassified",
    "Others",
]
TEST_TAXA = [t for t in TAXON_ORDER if t != "Others"]

COLORS = {
    "Enterobacteriaceae": "#E8735A",
    "Streptococcaceae": "#39B5B2",
    "Leuconostocaceae": "#1A7A78",
    "Enterococcaceae": "#7B6FAB",
    "Yersiniaceae": "#5A8FC4",
    "Unclassified": "#8E8E8E",
    "Others": "#D4CEBC",
}

FIG_SIZE = (5.2, 4.8)
BAR_WIDTH = 0.55
STAR_FONTSIZE = 9
MIN_HEIGHT_FOR_LABEL = 0.04


def sig_label(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def text_color_on_bg(hex_color):
    hex_color = hex_color.lstrip("#")
    r, g, b = [int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if luminance < 0.45 else "black"


df = pd.read_csv(INPUT_FILE, sep="\t")
df = df.set_index("Sample")
missing = [s for cols in GROUPS.values() for s in cols if s not in df.index]
if missing:
    raise ValueError(f"表中缺少样本: {missing}")

rel_abund = df[TAXON_ORDER].T  # taxon × sample

group_mean = pd.DataFrame(
    {g: rel_abund[cols].mean(axis=1) for g, cols in GROUPS.items()}
)
group_sd = pd.DataFrame(
    {g: rel_abund[cols].std(axis=1, ddof=1) for g, cols in GROUPS.items()}
)

# 组内再归一化，避免四舍五入导致总和略偏离 1
group_mean = group_mean.div(group_mean.sum(axis=0), axis=1)

n_comp = len(TREAT_GROUPS)
test_results = {}
for fam in TEST_TAXA:
    ck_vals = rel_abund.loc[fam, GROUPS[CONTROL_GROUP]].to_numpy(dtype=float)
    test_results[fam] = {}
    for g in TREAT_GROUPS:
        tr_vals = rel_abund.loc[fam, GROUPS[g]].to_numpy(dtype=float)
        if len(np.unique(np.concatenate([ck_vals, tr_vals]))) < 2:
            test_results[fam][g] = {"p_raw": np.nan, "p_bonf": np.nan, "sig": "na"}
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            _, p_raw = stats.mannwhitneyu(ck_vals, tr_vals, alternative="two-sided")
        if np.isnan(p_raw):
            test_results[fam][g] = {"p_raw": np.nan, "p_bonf": np.nan, "sig": "na"}
        else:
            p_bonf = min(p_raw * n_comp, 1.0)
            test_results[fam][g] = {
                "p_raw": p_raw,
                "p_bonf": p_bonf,
                "sig": sig_label(p_bonf),
            }

stat_rows = []
for fam in TAXON_ORDER:
    row = {"Family": fam}
    for g in GROUP_LABELS:
        m = group_mean.loc[fam, g]
        s = group_sd.loc[fam, g]
        row[f"{g} Mean±SD (%)"] = f"{m * 100:.2f}±{s * 100:.2f}"
        row[f"{g} Mean"] = m
        row[f"{g} SD"] = s
    kw_vals = [rel_abund.loc[fam, GROUPS[g]].to_numpy(dtype=float) for g in GROUP_LABELS]
    all_vals = np.concatenate(kw_vals)
    if len(np.unique(all_vals)) < 2:
        row["Kruskal-Wallis H"] = "na"
        row["Kruskal-Wallis p"] = "na"
        row["KW sig"] = "na"
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            kw_h, kw_p = stats.kruskal(*kw_vals)
        if np.isnan(kw_p):
            row["Kruskal-Wallis H"] = "na"
            row["Kruskal-Wallis p"] = "na"
            row["KW sig"] = "na"
        else:
            row["Kruskal-Wallis H"] = round(kw_h, 3)
            row["Kruskal-Wallis p"] = round(kw_p, 4)
            row["KW sig"] = sig_label(kw_p)
    if fam in test_results:
        for g in TREAT_GROUPS:
            r = test_results[fam][g]
            row[f"{CONTROL_GROUP} vs {g} p_raw"] = (
                round(r["p_raw"], 4) if not np.isnan(r["p_raw"]) else "na"
            )
            row[f"{CONTROL_GROUP} vs {g} p_bonf"] = (
                round(r["p_bonf"], 4) if not np.isnan(r["p_bonf"]) else "na"
            )
            row[f"{CONTROL_GROUP} vs {g} sig"] = r["sig"]
    stat_rows.append(row)

stat_df = pd.DataFrame(stat_rows)
mean_sd_long = []
for fam in TAXON_ORDER:
    for g in GROUP_LABELS:
        mean_sd_long.append({
            "Family": fam,
            "Group": g,
            "N": len(GROUPS[g]),
            "RA_percent_mean": group_mean.loc[fam, g] * 100,
            "RA_percent_sd": group_sd.loc[fam, g] * 100,
        })
mean_sd_df = pd.DataFrame(mean_sd_long)

mean_sd_df.to_csv(
    f"{OUT_PREFIX}_mean_SD_by_group.tsv", sep="\t", index=False
)
stat_df.to_csv(f"{OUT_PREFIX}_significance.tsv", sep="\t", index=False)
print(f"[OK] 统计表 -> {OUT_PREFIX}_mean_SD_by_group.tsv")
print(f"[OK] 显著性表 -> {OUT_PREFIX}_significance.tsv")

x = np.arange(len(GROUP_LABELS))
fig, ax = plt.subplots(figsize=FIG_SIZE)

seg_info = {fam: {} for fam in TAXON_ORDER}
bottoms = np.zeros(len(GROUP_LABELS))
for fam in TAXON_ORDER:
    heights = group_mean.loc[fam, GROUP_LABELS].to_numpy(dtype=float)
    ax.bar(
        x,
        heights,
        BAR_WIDTH,
        bottom=bottoms,
        color=COLORS[fam],
        edgecolor="white",
        linewidth=0.5,
        label=fam,
    )
    for gi, g in enumerate(GROUP_LABELS):
        seg_info[fam][g] = {
            "bottom": bottoms[gi],
            "top": bottoms[gi] + heights[gi],
            "center": bottoms[gi] + heights[gi] / 2,
            "height": heights[gi],
        }
    bottoms += heights

for fam in TEST_TAXA:
    fam_color = COLORS[fam]
    txt_color = text_color_on_bg(fam_color)
    for gi, g in enumerate(GROUP_LABELS):
        if g == CONTROL_GROUP:
            continue
        mark = test_results[fam][g]["sig"]
        if mark in ("ns", "na"):
            continue
        seg = seg_info[fam][g]
        if seg["height"] >= MIN_HEIGHT_FOR_LABEL:
            ax.text(
                gi,
                seg["center"],
                mark,
                ha="center",
                va="center",
                fontsize=STAR_FONTSIZE,
                fontweight="bold",
                color=txt_color,
            )
        else:
            box_pad = 0.012
            box_h = max(seg["height"], 0.018)
            box = mpatches.FancyBboxPatch(
                (gi - BAR_WIDTH / 2 - 0.06, max(seg["bottom"] - box_pad, 0)),
                BAR_WIDTH + 0.12,
                box_h + box_pad,
                boxstyle="square,pad=0",
                linewidth=1.1,
                edgecolor=fam_color,
                facecolor="none",
                zorder=4,
            )
            ax.add_patch(box)
            ax.text(
                gi + BAR_WIDTH / 2 + 0.16,
                max(seg["center"], 0.03),
                mark,
                ha="left",
                va="center",
                fontsize=STAR_FONTSIZE,
                fontweight="bold",
                color="black",
                zorder=5,
            )

ax.set_xlim(-0.55, len(GROUP_LABELS) - 0.05)
ax.set_ylim(0, 1.0)
ax.set_xticks(x)
ax.set_xticklabels(X_TICKS, fontsize=11)
ax.text(
    3.38,
    -0.075,
    "dpi",
    ha="left",
    va="top",
    fontsize=11,
    transform=ax.get_xaxis_transform(),
    clip_on=False,
)

ax.set_ylabel("Relative Abundance (%)", fontsize=11)
ax.set_xlabel("")
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
ax.tick_params(axis="both", labelsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

handles = [mpatches.Patch(color=COLORS[f], label=f) for f in TAXON_ORDER]
ax.legend(
    handles=handles,
    title="Family",
    title_fontsize=9,
    fontsize=8,
    loc="center left",
    bbox_to_anchor=(1.02, 0.55),
    frameon=False,
    handlelength=1.2,
    handleheight=1.0,
)

fig.tight_layout()
fig.savefig(f"{OUT_PREFIX}.png", dpi=300, bbox_inches="tight")
pdf_path = Path(f"{OUT_PREFIX}.pdf")
try:
    fig.savefig(pdf_path, bbox_inches="tight")
except PermissionError:
    pdf_path = Path(f"{OUT_PREFIX}_new.pdf")
    fig.savefig(pdf_path, bbox_inches="tight")
print(f"[OK] 图片已保存 -> {OUT_PREFIX}.png / {pdf_path.name}")
