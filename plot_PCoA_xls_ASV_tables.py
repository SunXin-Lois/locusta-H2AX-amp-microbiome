"""Draw Bray-Curtis PCoA + ANOSIM plots from ASV Excel tables.

Statistical methods are taken from plot_PCoA_Bray_ANOSIM_PTW-DPI.py
(Lingoes-corrected PCoA and rank-based ANOSIM). The figure layout follows
the requested ASV-level PCoA style.
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "output"
PERMUTATIONS = 999
SEED = 123

JOBS = [
    {
        "abundance": SCRIPT_DIR / "ASV_PTW135.xls",
        "prefix": "PCoA_ASV_PTW135",
        "groups": ["PTW", "1dpi", "3dpi", "5dpi"],
        "reference": "PTW",
        "contrast": "5dpi",
        "vertical": "3dpi",
    },
    {
        "abundance": SCRIPT_DIR / "ASV_PTW_RNAi_n1284.xls",
        "prefix": "PCoA_ASV_PTW_RNAi_n1284",
        "groups": ["PTW", "RNAi", "n1284"],
        "reference": "PTW",
        "contrast": "n1284",
        "vertical": "RNAi",
    },
]

GROUP_LABELS = {
    "PTW": "PTW",
    "1dpi": "1dpi",
    "3dpi": "3dpi",
    "5dpi": "5dpi",
    "RNAi": "RNAi",
    "n1284": "n1284",
}
GROUP_COLORS = {
    "PTW": "#B0B0B0",
    "1dpi": "#E6B422",
    "3dpi": "#5B9BD5",
    "5dpi": "#E07B39",
    "RNAi": "#E6B422",
    "n1284": "#E07B39",
}
GROUP_MARKERS = {
    "PTW": "D",
    "1dpi": "s",
    "3dpi": "^",
    "5dpi": "o",
    "RNAi": "s",
    "n1284": "o",
}


def load_pcoa_module():
    script_path = SCRIPT_DIR / "plot_PCoA_Bray_ANOSIM_PTW-DPI.py"
    spec = importlib.util.spec_from_file_location("pcoa_bray_anosim", script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sample_to_group(sample: str) -> str:
    return re.sub(r"[-_]\d+$", "", str(sample))


def read_abundance_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    engine = "xlrd" if suffix == ".xls" else "openpyxl"
    return pd.read_excel(path, engine=engine)


def load_asv_table(path: Path, expected_groups: list[str]) -> tuple[pd.DataFrame, pd.Series, pd.Series, int]:
    table = read_abundance_table(path)
    id_column = table.columns[0]
    sample_columns = [column for column in table.columns if column not in {id_column, "taxonomy"}]
    if not sample_columns:
        raise ValueError(f"No sample columns found in {path.name}")

    abundance = table.set_index(id_column)[sample_columns].apply(pd.to_numeric, errors="coerce")
    if abundance.isna().any().any():
        bad_columns = abundance.columns[abundance.isna().any()].tolist()
        raise ValueError(f"Non-numeric or missing abundance values in columns: {bad_columns}")
    if (abundance < 0).any().any():
        raise ValueError(f"{path.name} contains negative values.")
    if abundance.index.has_duplicates or abundance.columns.has_duplicates:
        raise ValueError(f"Duplicate ASV IDs or sample names in {path.name}.")

    groups = pd.Series([sample_to_group(sample) for sample in abundance.columns], index=abundance.columns, name="Group")
    unexpected = sorted(set(groups) - set(expected_groups))
    missing = [group for group in expected_groups if group not in set(groups)]
    if unexpected:
        raise ValueError(f"Unexpected sample groups in {path.name}: {unexpected}")
    if missing:
        raise ValueError(f"Missing groups in {path.name}: {missing}")

    input_asvs = abundance.shape[0]
    abundance_nonzero = abundance.loc[abundance.sum(axis=1) > 0].copy()
    if abundance_nonzero.shape[0] < 2:
        raise ValueError(f"Too few non-zero ASVs in {path.name} for PCoA.")

    sample_depths = abundance_nonzero.sum(axis=0).rename("Reads")
    relative = abundance_nonzero.T.astype(float)
    relative = relative.div(relative.sum(axis=1), axis=0)
    return relative, groups.loc[relative.index], sample_depths, input_asvs


def orient_coordinates(
    coordinates: np.ndarray,
    groups: pd.Series,
    reference: str,
    contrast: str,
    vertical: str,
) -> np.ndarray:
    oriented = coordinates.copy()
    group_array = groups.to_numpy()
    if reference in group_array:
        if oriented[group_array == reference, 0].mean() > oriented[:, 0].mean():
            oriented[:, 0] *= -1
        if oriented[group_array == reference, 1].mean() > 0:
            oriented[:, 1] *= -1
    elif contrast in group_array and vertical in group_array:
        if oriented[group_array == contrast, 0].mean() < oriented[group_array == vertical, 0].mean():
            oriented[:, 0] *= -1
    return oriented


def padded_limits(values: np.ndarray, padding: float = 0.14) -> tuple[float, float]:
    lower, upper = float(np.min(values)), float(np.max(values))
    span = upper - lower
    if span == 0:
        span = max(abs(lower), 1.0) * 0.1
    return lower - padding * span, upper + padding * span


def format_p_value(p_value: float) -> str:
    if p_value <= 0.001:
        return "0.001"
    return f"{p_value:.4f}".rstrip("0").rstrip(".")


def choose_text_position(coordinates: pd.DataFrame) -> tuple[float, float]:
    xs = coordinates["PCoA1"].to_numpy()
    ys = coordinates["PCoA2"].to_numpy()
    x_span = max(float(xs.max() - xs.min()), 1e-12)
    y_span = max(float(ys.max() - ys.min()), 1e-12)
    nx = (xs - xs.min()) / x_span
    ny = (ys - ys.min()) / y_span
    candidates = {
        (0.03, 0.96): (0.16, 0.88),
        (0.58, 0.96): (0.78, 0.88),
        (0.03, 0.28): (0.16, 0.12),
        (0.58, 0.28): (0.78, 0.12),
    }
    best_anchor = (0.58, 0.96)
    best_distance = -1.0
    for anchor, (cx, cy) in candidates.items():
        distance = float(np.min((nx - cx) ** 2 + (ny - cy) ** 2))
        if distance > best_distance:
            best_anchor = anchor
            best_distance = distance
    return best_anchor


def make_plot(
    coordinates: pd.DataFrame,
    group_order: list[str],
    explained: np.ndarray,
    anosim_r: float,
    anosim_p: float,
    png_path: Path,
    pdf_path: Path,
) -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 1.0,
        }
    )
    fig, ax = plt.subplots(figsize=(6.15, 4.85))
    ax.axhline(0, color="#C8C8C8", linestyle="--", linewidth=0.9, zorder=1)
    ax.axvline(0, color="#C8C8C8", linestyle="--", linewidth=0.9, zorder=1)

    for group in group_order:
        subset = coordinates.loc[coordinates["Group"] == group, ["PCoA1", "PCoA2"]].to_numpy()
        ax.scatter(
            subset[:, 0],
            subset[:, 1],
            s=58,
            marker=GROUP_MARKERS[group],
            c=GROUP_COLORS[group],
            edgecolors="none",
            label=GROUP_LABELS[group],
            zorder=3,
        )

    text_x, text_y = choose_text_position(coordinates)
    ax.text(
        text_x,
        text_y,
        f"ANOSIM\n$R$ = {anosim_r:.4f}, $P$ = {format_p_value(anosim_p)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
    )
    ax.set_title("PCoA on ASV level", fontsize=14, pad=10)
    ax.set_xlabel(f"PCoA1({explained[0]:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PCoA2({explained[1]:.1f}%)", fontsize=11)
    ax.set_xlim(padded_limits(coordinates["PCoA1"].to_numpy()))
    ax.set_ylim(padded_limits(coordinates["PCoA2"].to_numpy()))
    ax.tick_params(axis="both", colors="black", labelsize=9, width=1.0, length=3.5, direction="out")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("black")
        ax.spines[spine].set_linewidth(1.0)

    legend = ax.legend(
        title="Group",
        loc="center left",
        bbox_to_anchor=(1.02, 0.55),
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        handletextpad=0.5,
        borderaxespad=0.0,
        markerscale=1.05,
    )
    legend._legend_title_box.align = "left"
    fig.subplots_adjust(left=0.14, right=0.80, top=0.88, bottom=0.13)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, facecolor="white")
    fig.savefig(pdf_path, facecolor="white")
    plt.close(fig)


def run_job(job: dict, pcoa_module) -> str:
    abundance_path = job["abundance"]
    if not abundance_path.exists():
        raise FileNotFoundError(abundance_path)

    relative, groups, sample_depths, input_asvs = load_asv_table(abundance_path, job["groups"])
    bray_distance = squareform(pdist(relative.to_numpy(), metric="braycurtis"))
    pcoa = pcoa_module.pcoa_lingoes(bray_distance)
    coordinate_array = orient_coordinates(
        np.asarray(pcoa["coordinates"]),
        groups,
        reference=job["reference"],
        contrast=job["contrast"],
        vertical=job["vertical"],
    )
    anosim_r, anosim_p = pcoa_module.anosim(bray_distance, groups, PERMUTATIONS, SEED)
    explained = np.asarray(pcoa["explained"])

    coordinates = pd.DataFrame(
        {
            "Sample": relative.index,
            "Group": groups.loc[relative.index].to_numpy(),
            "Group_label": [GROUP_LABELS[group] for group in groups.loc[relative.index]],
            "PCoA1": coordinate_array[:, 0],
            "PCoA2": coordinate_array[:, 1],
            "Reads": sample_depths.loc[relative.index].to_numpy(),
        }
    )

    prefix = job["prefix"]
    png_path = OUTPUT_DIR / f"{prefix}.png"
    pdf_path = OUTPUT_DIR / f"{prefix}.pdf"
    coordinates_path = OUTPUT_DIR / f"{prefix}_coordinates.tsv"
    summary_path = OUTPUT_DIR / f"{prefix}_summary.txt"

    make_plot(coordinates, job["groups"], explained, anosim_r, anosim_p, png_path, pdf_path)
    coordinates.to_csv(coordinates_path, sep="\t", index=False, float_format="%.10g")

    group_counts = groups.value_counts().reindex(job["groups"])
    count_text = ", ".join(f"{GROUP_LABELS[group]}={int(group_counts[group])}" for group in job["groups"])
    p_text = f"= {format_p_value(anosim_p)}"
    original_eigenvalues = np.asarray(pcoa["original_eigenvalues"])
    positive_original = original_eigenvalues[original_eigenvalues > 1e-10]
    negative_original = original_eigenvalues[original_eigenvalues < -1e-10]
    summary_lines = [
        "Analysis: Bray-Curtis PCoA with Lingoes correction + Bray-Curtis ANOSIM",
        f"Method source: {SCRIPT_DIR.as_posix()}/plot_PCoA_Bray_ANOSIM_PTW-DPI.py",
        f"Abundance input: {abundance_path.as_posix()}",
        f"Samples: {relative.shape[0]} ({count_text})",
        f"ASVs in abundance input: {input_asvs}",
        f"All-zero ASVs removed: {input_asvs - relative.shape[1]}",
        f"Non-zero ASVs used: {relative.shape[1]}",
        f"Sample depth range: {int(sample_depths.min())}-{int(sample_depths.max())} reads",
        "Distance input: within-sample relative abundance",
        "Distance: Bray-Curtis",
        f"Lingoes constant: {float(pcoa['lingoes_constant']):.10g}",
        f"Original positive eigenvalues: {len(positive_original)}",
        f"Original negative eigenvalues: {len(negative_original)}",
        f"PCoA1 corrected variance: {explained[0]:.3f}%",
        f"PCoA2 corrected variance: {explained[1]:.3f}%",
        f"Cumulative PCoA1-PCoA2 corrected variance: {explained[:2].sum():.3f}%",
        f"ANOSIM (Bray-Curtis, {PERMUTATIONS} permutations, seed={SEED}): R = {anosim_r:.4f}, P {p_text}",
        f"PDF: {pdf_path.as_posix()}",
        f"PNG: {png_path.as_posix()}",
    ]
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    return "\n".join(summary_lines)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pcoa_module = load_pcoa_module()
    reports = [run_job(job, pcoa_module) for job in JOBS]
    print("\n\n".join(reports))


if __name__ == "__main__":
    main()
