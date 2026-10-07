"""Bray-Curtis PCoA for PTW, 1 dpi, 3 dpi and 5 dpi ASV communities.

Abundances come from the even-depth ASV table. Taxonomy is joined from the
corrected taxonomy table and is used only for biological interpretation, not
as input to the distance calculation.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata


SCRIPT_DIR = Path(__file__).resolve().parent
EXPECTED_GROUPS = ["PTW", "1dpi", "3dpi", "5dpi"]
GROUP_LABELS = {"PTW": "PTW", "1dpi": "1 dpi", "3dpi": "3 dpi", "5dpi": "5 dpi"}
GROUP_COLORS = {"PTW": "#999999", "1dpi": "#E69F00", "3dpi": "#56B4E9", "5dpi": "#D55E00"}
GROUP_MARKERS = {"PTW": "D", "1dpi": "s", "3dpi": "^", "5dpi": "o"}
TAXONOMY_COLUMNS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def resolve_from_script(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else SCRIPT_DIR / path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--abundance",
        default="table_f_even_depth_w_tax_LmH_PTW-DPI.tsv",
        help="Even-depth ASV abundance table (ASVs in rows, samples in columns).",
    )
    parser.add_argument(
        "--taxonomy",
        default="table_filtered_w_tax_订正肠杆菌科分类_改分组_改样本名.tsv",
        help="ASV table containing corrected Domain-Species taxonomy columns.",
    )
    parser.add_argument(
        "--result-dir",
        default="output/PCoA_PTW-DPI",
        help="Directory for PNG, coordinates, taxonomy tables, log and summary.",
    )
    parser.add_argument(
        "--pdf",
        default="output/pdf/Fig1d_PCoA_Bray_Lingoes.pdf",
        help="Final vector PDF path.",
    )
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=123)
    return parser.parse_args()


def load_even_depth_abundance(path: Path) -> tuple[pd.DataFrame, pd.Series, pd.Series, int]:
    abundance = pd.read_csv(path, sep="\t", index_col=0)
    if abundance.empty:
        raise ValueError("The abundance table is empty.")
    numeric = abundance.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        bad_columns = numeric.columns[numeric.isna().any()].tolist()
        raise ValueError(f"Non-numeric or missing abundance values in columns: {bad_columns}")
    if (numeric < 0).any().any():
        raise ValueError("The abundance table contains negative values.")
    if numeric.index.has_duplicates or numeric.columns.has_duplicates:
        raise ValueError("Duplicate ASV IDs or sample names were found.")

    sample_names = pd.Index(numeric.columns.astype(str))
    group_values = [re.sub(r"-[0-9]+$", "", sample) for sample in sample_names]
    group = pd.Series(group_values, index=sample_names, name="Group")
    unexpected = sorted(set(group) - set(EXPECTED_GROUPS))
    if unexpected:
        raise ValueError(f"Unexpected sample groups: {unexpected}")
    group_counts = group.value_counts().reindex(EXPECTED_GROUPS, fill_value=0)
    if not (group_counts == 5).all():
        raise ValueError(f"Expected five samples per group, observed: {group_counts.to_dict()}")

    input_asvs = numeric.shape[0]
    abundance_nonzero = numeric.loc[numeric.sum(axis=1) > 0].copy()
    if abundance_nonzero.shape[0] < 2:
        raise ValueError("Too few non-zero ASVs for PCoA.")
    sample_depths = abundance_nonzero.sum(axis=0).rename("Reads")
    sample_by_asv = abundance_nonzero.T.astype(float)
    relative = sample_by_asv.div(sample_by_asv.sum(axis=1), axis=0)
    return relative, group, sample_depths, input_asvs


def load_taxonomy(path: Path, asv_ids: pd.Index) -> pd.DataFrame:
    taxonomy_table = pd.read_csv(path, sep="\t", index_col=0, dtype=str, keep_default_na=False)
    missing_columns = [column for column in TAXONOMY_COLUMNS if column not in taxonomy_table.columns]
    if missing_columns:
        raise ValueError(f"Missing taxonomy columns: {missing_columns}")
    if taxonomy_table.index.has_duplicates:
        raise ValueError("The taxonomy table contains duplicate ASV IDs.")
    missing_asvs = asv_ids.difference(taxonomy_table.index)
    if len(missing_asvs):
        raise ValueError(f"Taxonomy is missing for {len(missing_asvs)} ASVs, including {missing_asvs[:5].tolist()}")
    return taxonomy_table.loc[asv_ids, TAXONOMY_COLUMNS].copy()


def centered_gram(distance_matrix: np.ndarray) -> np.ndarray:
    n_samples = distance_matrix.shape[0]
    centering = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    return -0.5 * centering @ (distance_matrix**2) @ centering


def eigendecompose_gram(gram: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    return eigenvalues[order], eigenvectors[:, order]


def pcoa_lingoes(distance_matrix: np.ndarray) -> dict[str, np.ndarray | float]:
    original_eigenvalues, _ = eigendecompose_gram(centered_gram(distance_matrix))
    lingoes_constant = max(0.0, -float(original_eigenvalues.min()))
    if lingoes_constant > 0:
        corrected_squared = distance_matrix**2 + 2.0 * lingoes_constant
        np.fill_diagonal(corrected_squared, 0.0)
        corrected_distance = np.sqrt(np.maximum(corrected_squared, 0.0))
    else:
        corrected_distance = distance_matrix.copy()

    corrected_eigenvalues, corrected_eigenvectors = eigendecompose_gram(centered_gram(corrected_distance))
    tolerance = max(1.0, abs(corrected_eigenvalues[0])) * 1e-10
    positive = corrected_eigenvalues > tolerance
    if positive.sum() < 2:
        raise ValueError("PCoA produced fewer than two positive axes.")
    coordinates_all = corrected_eigenvectors[:, positive] * np.sqrt(corrected_eigenvalues[positive])
    positive_eigenvalues = corrected_eigenvalues[positive]
    explained = 100.0 * positive_eigenvalues / positive_eigenvalues.sum()
    return {
        "coordinates": coordinates_all[:, :2],
        "explained": explained,
        "original_eigenvalues": original_eigenvalues,
        "corrected_eigenvalues": corrected_eigenvalues,
        "lingoes_constant": lingoes_constant,
        "corrected_distance": corrected_distance,
    }


def orient_coordinates(coordinates: np.ndarray, groups: pd.Series) -> np.ndarray:
    oriented = coordinates.copy()
    group_array = groups.to_numpy()
    if oriented[group_array == "5dpi", 0].mean() < oriented[group_array == "PTW", 0].mean():
        oriented[:, 0] *= -1
    if oriented[group_array == "3dpi", 1].mean() < oriented[group_array == "PTW", 1].mean():
        oriented[:, 1] *= -1
    return oriented


def anosim(distance_matrix: np.ndarray, groups: pd.Series, permutations: int, seed: int) -> tuple[float, float]:
    n_samples = distance_matrix.shape[0]
    upper = np.triu_indices(n_samples, 1)
    ranks = rankdata(distance_matrix[upper], method="average")
    group_array = groups.to_numpy()

    def statistic(labels: np.ndarray) -> float:
        within = labels[upper[0]] == labels[upper[1]]
        if not within.any() or within.all():
            raise ValueError("ANOSIM requires both within-group and between-group distances.")
        denominator = n_samples * (n_samples - 1) / 4.0
        return float((ranks[~within].mean() - ranks[within].mean()) / denominator)

    observed = statistic(group_array)
    rng = np.random.default_rng(seed)
    exceedances = 0
    for _ in range(permutations):
        if statistic(rng.permutation(group_array)) >= observed - 1e-12:
            exceedances += 1
    p_value = (exceedances + 1) / (permutations + 1)
    return observed, p_value


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted_ranked = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def asv_axis_associations(
    relative: pd.DataFrame,
    coordinates: np.ndarray,
    taxonomy: pd.DataFrame,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    x = coordinates[:, :2] - coordinates[:, :2].mean(axis=0)
    x_norm = np.linalg.norm(x, axis=0)
    x_standardized = x / x_norm

    y = relative.to_numpy() - relative.to_numpy().mean(axis=0)
    y_norm = np.linalg.norm(y, axis=0)
    valid = y_norm > 0
    y_standardized = np.zeros_like(y)
    y_standardized[:, valid] = y[:, valid] / y_norm[valid]

    correlations = x_standardized.T @ y_standardized
    observed_r2 = np.sum(correlations**2, axis=0)
    rng = np.random.default_rng(seed + 1)
    exceedances = np.zeros(relative.shape[1], dtype=int)
    for _ in range(permutations):
        permuted_correlations = x_standardized[rng.permutation(relative.shape[0]), :].T @ y_standardized
        permuted_r2 = np.sum(permuted_correlations**2, axis=0)
        exceedances += permuted_r2 >= observed_r2 - 1e-12
    p_values = (exceedances + 1) / (permutations + 1)

    result = pd.DataFrame(
        {
            "ASV": relative.columns,
            "PCoA1_correlation": correlations[0],
            "PCoA2_correlation": correlations[1],
            "R2": observed_r2,
            "Permutation_P": p_values,
            "FDR_BH": benjamini_hochberg(p_values),
            "Mean_relative_abundance": relative.mean(axis=0).to_numpy(),
            "Maximum_relative_abundance": relative.max(axis=0).to_numpy(),
        }
    ).set_index("ASV")
    result = result.join(taxonomy, how="left")
    return result.sort_values(["R2", "Permutation_P"], ascending=[False, True])


def dominant_asvs_by_sample(relative: pd.DataFrame, taxonomy: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for sample in relative.index:
        asv = relative.loc[sample].idxmax()
        record: dict[str, object] = {
            "Sample": sample,
            "Group": GROUP_LABELS[groups.loc[sample]],
            "Dominant_ASV": asv,
            "Relative_abundance": float(relative.loc[sample, asv]),
        }
        record.update(taxonomy.loc[asv].to_dict())
        records.append(record)
    return pd.DataFrame(records)


def top_group_taxa(relative: pd.DataFrame, taxonomy: pd.DataFrame, groups: pd.Series, n_top: int = 50) -> pd.DataFrame:
    means = pd.DataFrame(index=relative.columns)
    for group in EXPECTED_GROUPS:
        means[f"Mean_{GROUP_LABELS[group]}"] = relative.loc[groups == group].mean(axis=0)
    means["Overall_mean"] = relative.mean(axis=0)
    means["Maximum_group_mean"] = means[[f"Mean_{GROUP_LABELS[g]}" for g in EXPECTED_GROUPS]].max(axis=1)
    means["Highest_group"] = means[[f"Mean_{GROUP_LABELS[g]}" for g in EXPECTED_GROUPS]].idxmax(axis=1).str.replace("Mean_", "", regex=False)
    means.index.name = "ASV"
    result = means.join(taxonomy, how="left")
    return result.sort_values("Overall_mean", ascending=False).head(n_top)


def padded_limits(values: np.ndarray, padding: float = 0.12) -> tuple[float, float]:
    lower, upper = float(np.min(values)), float(np.max(values))
    span = upper - lower
    if span == 0:
        span = max(abs(lower), 1.0) * 0.1
    return lower - padding * span, upper + padding * span


def make_plot(
    coordinates: pd.DataFrame,
    explained: np.ndarray,
    anosim_r: float,
    anosim_p: float,
    png_path: Path,
    pdf_path: Path,
) -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
        }
    )
    fig, ax = plt.subplots(figsize=(4.25, 4.15))
    for group in EXPECTED_GROUPS:
        subset = coordinates.loc[coordinates["Group"] == group, ["PCoA1", "PCoA2"]].to_numpy()
        ax.scatter(
            subset[:, 0],
            subset[:, 1],
            s=52,
            marker=GROUP_MARKERS[group],
            c=GROUP_COLORS[group],
            edgecolors="white",
            linewidths=0.45,
            label=GROUP_LABELS[group],
            zorder=3,
        )

    p_text = "< 0.001" if anosim_p < 0.001 else f"= {anosim_p:.3f}"
    ax.text(
        0.03,
        0.96,
        rf"$R = {anosim_r:.3f},\ P {p_text}$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
    )
    ax.set_xlabel(f"PCoA1 ({explained[0]:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PCoA2 ({explained[1]:.1f}%)", fontsize=10)
    ax.set_xlim(padded_limits(coordinates["PCoA1"].to_numpy()))
    ax.set_ylim(padded_limits(coordinates["PCoA2"].to_numpy()))
    ax.tick_params(axis="both", colors="black", labelsize=8, width=0.8, length=3)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.8)
    ax.set_box_aspect(1)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=4,
        frameon=False,
        fontsize=8,
        handletextpad=0.35,
        columnspacing=0.9,
        borderaxespad=0.0,
    )
    fig.text(0.028, 0.965, "D", ha="left", va="top", fontsize=12, fontweight="bold")
    fig.subplots_adjust(left=0.20, right=0.97, top=0.96, bottom=0.23)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, facecolor="white")
    fig.savefig(pdf_path, facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    abundance_path = resolve_from_script(args.abundance)
    taxonomy_path = resolve_from_script(args.taxonomy)
    result_dir = resolve_from_script(args.result_dir)
    pdf_path = resolve_from_script(args.pdf)
    if not abundance_path.exists():
        raise FileNotFoundError(abundance_path)
    if not taxonomy_path.exists():
        raise FileNotFoundError(taxonomy_path)
    if args.permutations < 1:
        raise ValueError("At least one permutation is required.")

    result_dir.mkdir(parents=True, exist_ok=True)
    relative, groups, sample_depths, input_asvs = load_even_depth_abundance(abundance_path)
    taxonomy = load_taxonomy(taxonomy_path, relative.columns)

    bray_distance = squareform(pdist(relative.to_numpy(), metric="braycurtis"))
    pcoa = pcoa_lingoes(bray_distance)
    coordinate_array = orient_coordinates(np.asarray(pcoa["coordinates"]), groups)
    anosim_r, anosim_p = anosim(bray_distance, groups, args.permutations, args.seed)

    coordinates = pd.DataFrame(
        {
            "Sample": relative.index,
            "Group": groups.loc[relative.index].to_numpy(),
            "Group_label": [GROUP_LABELS[group] for group in groups.loc[relative.index]],
            "PCoA1": coordinate_array[:, 0],
            "PCoA2": coordinate_array[:, 1],
        }
    ).set_index("Sample", drop=False)

    associations = asv_axis_associations(
        relative,
        coordinate_array,
        taxonomy,
        permutations=args.permutations,
        seed=args.seed,
    )
    dominant = dominant_asvs_by_sample(relative, taxonomy, groups)
    top_taxa = top_group_taxa(relative, taxonomy, groups, n_top=50)

    prefix = "Fig1d_PCoA_Bray_Lingoes"
    png_path = result_dir / f"{prefix}.png"
    coordinates_path = result_dir / f"{prefix}_coordinates.tsv"
    associations_path = result_dir / f"{prefix}_top_ASV_axis_associations.tsv"
    dominant_path = result_dir / f"{prefix}_dominant_ASVs_by_sample.tsv"
    top_taxa_path = result_dir / f"{prefix}_top_taxa_by_group.tsv"
    summary_path = result_dir / f"{prefix}_summary.txt"

    make_plot(coordinates, np.asarray(pcoa["explained"]), anosim_r, anosim_p, png_path, pdf_path)
    coordinates.to_csv(coordinates_path, sep="\t", index=False, float_format="%.10g")
    associations.head(50).to_csv(associations_path, sep="\t", float_format="%.10g")
    dominant.to_csv(dominant_path, sep="\t", index=False, float_format="%.10g")
    top_taxa.to_csv(top_taxa_path, sep="\t", float_format="%.10g")

    original_eigenvalues = np.asarray(pcoa["original_eigenvalues"])
    positive_original = original_eigenvalues[original_eigenvalues > 1e-10]
    negative_original = original_eigenvalues[original_eigenvalues < -1e-10]
    explained = np.asarray(pcoa["explained"])
    top_association = associations.iloc[0]
    max_dominance = dominant.loc[dominant["Relative_abundance"].idxmax()]
    p_summary = "< 0.001" if anosim_p < 0.001 else f"= {anosim_p:.3f}"
    summary_lines = [
        "Analysis: Bray-Curtis PCoA with Lingoes correction + Bray-Curtis ANOSIM",
        f"Abundance input (even depth): {abundance_path.as_posix()}",
        f"Taxonomy input (corrected Domain-Species): {taxonomy_path.as_posix()}",
        f"Samples: {relative.shape[0]} (PTW=5, 1 dpi=5, 3 dpi=5, 5 dpi=5)",
        f"ASVs in abundance input: {input_asvs}",
        f"All-zero ASVs removed: {input_asvs - relative.shape[1]}",
        f"Non-zero ASVs used: {relative.shape[1]}",
        f"Sample depth range: {int(sample_depths.min())}-{int(sample_depths.max())} reads",
        "Distance input: within-sample relative abundance (equivalent to counts here because all depths are equal)",
        "Distance: Bray-Curtis",
        f"Lingoes constant: {float(pcoa['lingoes_constant']):.10g}",
        f"Original positive eigenvalues: {len(positive_original)}",
        f"Original negative eigenvalues: {len(negative_original)}",
        f"Absolute negative/positive eigenvalue sum: {abs(negative_original).sum() / positive_original.sum():.6%}",
        f"PCoA1 corrected variance: {explained[0]:.3f}%",
        f"PCoA2 corrected variance: {explained[1]:.3f}%",
        f"Cumulative PCoA1-PCoA2 corrected variance: {explained[:2].sum():.3f}%",
        f"ANOSIM (Bray-Curtis, {args.permutations} permutations, seed={args.seed}): R = {anosim_r:.3f}, P {p_summary}",
        "Group ellipses omitted because PTW, 1 dpi and 3 dpi had near-singular covariance on the first two axes.",
        "Sample coordinates were not jittered or otherwise altered.",
        "PCoA axis signs were oriented only for consistent display; distances and statistics are unchanged.",
        f"Taxonomy coverage: {taxonomy.shape[0]}/{relative.shape[1]} non-zero ASVs matched.",
        f"Strongest two-axis ASV association: {associations.index[0]} (R2={top_association['R2']:.3f}, P={top_association['Permutation_P']:.3f}, {top_association['Species']})",
        f"Highest single-sample dominance: {max_dominance['Sample']} / {max_dominance['Dominant_ASV']} / {max_dominance['Relative_abundance']:.2%} / {max_dominance['Species']}",
        f"PDF: {pdf_path.as_posix()}",
        f"PNG: {png_path.as_posix()}",
    ]
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
