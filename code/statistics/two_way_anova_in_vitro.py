"""In vitro endpoints (Fig. 5): two-way ANOVA with planned comparisons.

Model per endpoint: value ~ 4-HC concentration (categorical) * iPTH, with
sum-to-zero contrasts. Brown-Forsythe testing across factorial cells selects
classical (p >= 0.05) or HC3 (p < 0.05) covariance. Planned comparisons:
iPTH alone at 0 uM (unadjusted), 4-HC alone vs 0 uM (Dunnett), iPTH vs no iPTH
at each nonzero 4-HC concentration (multivariate-t), and combined treatment vs
baseline control (exploratory, multivariate-t).

Usage (from the repository root):
    python code/statistics/two_way_anova_in_vitro.py
    python code/statistics/build_figure_significance.py

Input:  data/statistics/in_vitro/*.csv
Output: results/statistics/in_vitro_two_way_anova/
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import build_design_matrices
from scipy import stats


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = REPO_ROOT / "data" / "statistics" / "in_vitro"
DEFAULT_OUTPUT = REPO_ROOT / "results" / "statistics" / "in_vitro_two_way_anova"
SEED = 20260615
QMC_POWER = 19


def safe_name(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)


def display_path(path: Path) -> str:
    """Path relative to the repository root when possible, with forward slashes."""
    try:
        return path.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def p_symbol(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "NA"
    if p_value < 0.0001:
        return "****"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def prepare_data(path: Path) -> tuple[pd.DataFrame, str]:
    raw = pd.read_csv(path)
    required = {"CY_conc", "P_trea"}
    if not required.issubset(raw.columns):
        raise ValueError(f"{path.name}: missing CY_conc or P_trea")
    value_columns = [column for column in raw.columns if column not in required]
    if len(value_columns) != 1:
        raise ValueError(
            f"{path.name}: expected one outcome column, found {value_columns}"
        )
    source_value_column = value_columns[0]
    data = raw.rename(
        columns={
            "CY_conc": "concentration_um",
            "P_trea": "teriparatide",
            source_value_column: "value",
        }
    ).copy()
    data["value"] = pd.to_numeric(data["value"], errors="raise")
    data["concentration_um"] = pd.to_numeric(
        data["concentration_um"], errors="raise"
    ).astype(int)
    data["teriparatide"] = data["teriparatide"].astype(str)
    if set(data["teriparatide"]) != {"No", "Yes"}:
        raise ValueError(
            f"{path.name}: expected teriparatide levels No/Yes, "
            f"found {sorted(data['teriparatide'].unique())}"
        )
    concentrations = sorted(data["concentration_um"].unique())
    if 0 not in concentrations:
        raise ValueError(f"{path.name}: concentration 0 is required")
    data["concentration"] = pd.Categorical(
        data["concentration_um"].astype(str),
        categories=[str(value) for value in concentrations],
        ordered=True,
    )
    data["teriparatide"] = pd.Categorical(
        data["teriparatide"], categories=["No", "Yes"], ordered=True
    )
    data.insert(
        0,
        "experimental_unit_id",
        [f"{path.stem}_EU_{index:03d}" for index in range(1, len(data) + 1)],
    )
    data["technical_replicates_averaged"] = 1
    data["analysis_unit"] = "independent experimental observation"
    return data, source_value_column


def simulate_max_abs_t(
    correlation: np.ndarray,
    df_residual: float,
    seed: int,
) -> np.ndarray:
    correlation = np.asarray(correlation, dtype=float)
    correlation = (correlation + correlation.T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    keep = eigenvalues > 1e-10
    loading = eigenvectors[:, keep] @ np.diag(np.sqrt(eigenvalues[keep]))
    uniforms = stats.qmc.Sobol(
        d=loading.shape[1] + 1,
        scramble=True,
        seed=seed,
    ).random_base2(QMC_POWER)
    eps = np.finfo(float).eps
    uniforms = np.clip(uniforms, eps, 1 - eps)
    normals = stats.norm.ppf(uniforms[:, :-1]) @ loading.T
    chi_square = stats.chi2.ppf(uniforms[:, -1], df_residual)
    simulated_t = normals / np.sqrt(chi_square[:, None] / df_residual)
    return np.max(np.abs(simulated_t), axis=1)


def contrast_family(
    endpoint: str,
    family: str,
    labels: list[str],
    contrast_matrix: np.ndarray,
    cell_means: np.ndarray,
    cell_covariance: np.ndarray,
    df_residual: float,
    adjustment: str,
    seed: int,
) -> pd.DataFrame:
    estimates = contrast_matrix @ cell_means
    covariance = contrast_matrix @ cell_covariance @ contrast_matrix.T
    standard_errors = np.sqrt(np.diag(covariance))
    t_values = estimates / standard_errors
    raw_p = 2 * stats.t.sf(np.abs(t_values), df_residual)

    if adjustment == "none" or len(labels) == 1:
        adjusted_p = raw_p.copy()
        critical = stats.t.ppf(0.975, df_residual)
        method = "none; single prespecified comparison"
    else:
        correlation = covariance / np.outer(standard_errors, standard_errors)
        max_abs_t = simulate_max_abs_t(correlation, df_residual, seed)
        adjusted_p = np.asarray(
            [
                (np.sum(max_abs_t >= abs(value)) + 1)
                / (len(max_abs_t) + 1)
                for value in t_values
            ]
        )
        adjusted_p = np.maximum(adjusted_p, raw_p)
        critical = float(np.quantile(max_abs_t, 0.95))
        method = adjustment

    return pd.DataFrame(
        {
            "endpoint": endpoint,
            "family": family,
            "comparison": labels,
            "estimate": estimates,
            "std_error": standard_errors,
            "df": df_residual,
            "t_value": t_values,
            "p_raw": raw_p,
            "p_adjusted": adjusted_p,
            "adjustment": method,
            "ci95_low": estimates - critical * standard_errors,
            "ci95_high": estimates + critical * standard_errors,
            "ci_type": (
                "pointwise 95%" if adjustment == "none" or len(labels) == 1
                else "simultaneous 95%"
            ),
            "significance": [p_symbol(value) for value in adjusted_p],
        }
    )


def build_cell_grid(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for concentration in data["concentration"].cat.categories:
        for treatment in ["No", "Yes"]:
            rows.append(
                {
                    "concentration": concentration,
                    "teriparatide": treatment,
                }
            )
    grid = pd.DataFrame(rows)
    grid["concentration"] = pd.Categorical(
        grid["concentration"],
        categories=data["concentration"].cat.categories,
        ordered=True,
    )
    grid["teriparatide"] = pd.Categorical(
        grid["teriparatide"], categories=["No", "Yes"], ordered=True
    )
    return grid


def vector_for(
    cell_lookup: dict[tuple[str, str], int],
    positive: tuple[str, str],
    negative: tuple[str, str],
) -> np.ndarray:
    vector = np.zeros(len(cell_lookup))
    vector[cell_lookup[positive]] = 1.0
    vector[cell_lookup[negative]] = -1.0
    return vector


def analyse_endpoint(
    path: Path,
    output_dir: Path,
    endpoint_index: int,
) -> dict[str, pd.DataFrame | dict]:
    endpoint = path.stem
    data, source_value_column = prepare_data(path)
    standardized_dir = output_dir / "analysis_level_data"
    standardized_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(
        standardized_dir / f"{safe_name(endpoint)}_analysis_level.csv",
        index=False,
    )

    formula = (
        "value ~ C(concentration, Sum) * C(teriparatide, Sum)"
    )
    classical_model = smf.ols(formula, data=data).fit()

    cell_groups = [
        group["value"].to_numpy()
        for _, group in data.groupby(
            ["concentration", "teriparatide"],
            observed=True,
            sort=False,
        )
    ]
    brown_forsythe = stats.levene(*cell_groups, center="median")
    use_hc3 = bool(brown_forsythe.pvalue < 0.05)
    covariance_method = "HC3" if use_hc3 else "classical"
    inferential_model = (
        classical_model.get_robustcov_results(cov_type="HC3", use_t=True)
        if use_hc3
        else classical_model
    )

    anova = sm.stats.anova_lm(
        classical_model,
        typ=3,
        robust="hc3" if use_hc3 else None,
    ).reset_index()
    anova = anova.rename(
        columns={
            "index": "term",
            "sum_sq": "sum_squares",
            "F": "f_value",
            "PR(>F)": "p_value",
        }
    )
    anova["term"] = anova["term"].replace(
        {
            "C(concentration, Sum)": "4-HC concentration",
            "C(teriparatide, Sum)": "teriparatide",
            "C(concentration, Sum):C(teriparatide, Sum)": (
                "4-HC concentration x teriparatide"
            ),
        }
    )
    anova.insert(0, "endpoint", endpoint)
    anova["covariance_method"] = covariance_method

    grid = build_cell_grid(data)
    design = np.asarray(
        build_design_matrices(
            [classical_model.model.data.design_info], grid
        )[0]
    )
    beta = np.asarray(inferential_model.params)
    beta_covariance = np.asarray(inferential_model.cov_params())
    cell_means = design @ beta
    cell_covariance = design @ beta_covariance @ design.T
    cell_lookup = {
        (str(row.concentration), str(row.teriparatide)): index
        for index, row in grid.iterrows()
    }
    concentrations = [str(value) for value in data["concentration"].cat.categories]
    nonzero = [value for value in concentrations if value != "0"]

    contrasts = []
    teriparatide_alone_matrix = np.vstack(
        [vector_for(cell_lookup, ("0", "Yes"), ("0", "No"))]
    )
    contrasts.append(
        contrast_family(
            endpoint,
            "teriparatide_alone_at_0uM",
            ["0 uM: teriparatide - no teriparatide"],
            teriparatide_alone_matrix,
            cell_means,
            cell_covariance,
            inferential_model.df_resid,
            "none",
            SEED + endpoint_index * 10,
        )
    )

    dunnett_matrix = np.vstack(
        [
            vector_for(cell_lookup, (concentration, "No"), ("0", "No"))
            for concentration in nonzero
        ]
    )
    contrasts.append(
        contrast_family(
            endpoint,
            "4HC_alone_Dunnett",
            [
                f"{concentration} uM - 0 uM, no teriparatide"
                for concentration in nonzero
            ],
            dunnett_matrix,
            cell_means,
            cell_covariance,
            inferential_model.df_resid,
            "Dunnett single-step multivariate-t",
            SEED + endpoint_index * 10 + 1,
        )
    )

    rescue_matrix = np.vstack(
        [
            vector_for(
                cell_lookup,
                (concentration, "Yes"),
                (concentration, "No"),
            )
            for concentration in nonzero
        ]
    )
    contrasts.append(
        contrast_family(
            endpoint,
            "rescue_MVT",
            [
                f"{concentration} uM: teriparatide - no teriparatide"
                for concentration in nonzero
            ],
            rescue_matrix,
            cell_means,
            cell_covariance,
            inferential_model.df_resid,
            "single-step multivariate-t",
            SEED + endpoint_index * 10 + 2,
        )
    )

    diagonal_matrix = np.vstack(
        [
            vector_for(
                cell_lookup,
                (concentration, "Yes"),
                ("0", "No"),
            )
            for concentration in nonzero
        ]
    )
    contrasts.append(
        contrast_family(
            endpoint,
            "exploratory_diagonal_MVT",
            [
                f"{concentration} uM + teriparatide - baseline control"
                for concentration in nonzero
            ],
            diagonal_matrix,
            cell_means,
            cell_covariance,
            inferential_model.df_resid,
            "single-step multivariate-t",
            SEED + endpoint_index * 10 + 3,
        )
    )
    contrasts_df = pd.concat(contrasts, ignore_index=True)
    contrasts_df["covariance_method"] = covariance_method

    descriptions = (
        data.groupby(
            ["concentration_um", "teriparatide"],
            observed=True,
        )["value"]
        .agg(n="count", mean="mean", sd="std", median="median", minimum="min", maximum="max")
        .reset_index()
    )
    descriptions.insert(0, "endpoint", endpoint)

    influence = classical_model.get_influence()
    studentized = influence.resid_studentized_external
    cooks = influence.cooks_distance[0]
    shapiro = stats.shapiro(classical_model.resid)
    diagnostics = pd.DataFrame(
        [
            {
                "endpoint": endpoint,
                "diagnostic": "Shapiro-Wilk residual normality",
                "statistic": shapiro.statistic,
                "p_value": shapiro.pvalue,
                "flag": shapiro.pvalue < 0.05,
                "decision_note": (
                    "Normality flag; inspect residual and Q-Q plots."
                ),
            },
            {
                "endpoint": endpoint,
                "diagnostic": "Brown-Forsythe equal variance across cells",
                "statistic": brown_forsythe.statistic,
                "p_value": brown_forsythe.pvalue,
                "flag": use_hc3,
                "decision_note": (
                    "HC3 used for all endpoint inference."
                    if use_hc3
                    else "Classical covariance retained."
                ),
            },
            {
                "endpoint": endpoint,
                "diagnostic": "Maximum absolute externally studentized residual",
                "statistic": np.nanmax(np.abs(studentized)),
                "p_value": np.nan,
                "flag": np.nanmax(np.abs(studentized)) > 3,
                "decision_note": "Flag threshold: absolute value > 3.",
            },
            {
                "endpoint": endpoint,
                "diagnostic": "Maximum Cook's distance",
                "statistic": np.nanmax(cooks),
                "p_value": np.nan,
                "flag": np.nanmax(cooks) > 4 / len(data),
                "decision_note": f"Flag threshold: 4/n = {4 / len(data):.6g}.",
            },
            {
                "endpoint": endpoint,
                "diagnostic": "Batch effect",
                "statistic": np.nan,
                "p_value": np.nan,
                "flag": False,
                "decision_note": (
                    "Not estimable/not applicable: source data contain no "
                    "batch variable and each row was confirmed independent."
                ),
            },
        ]
    )

    observation_diagnostics = data[
        [
            "experimental_unit_id",
            "concentration_um",
            "teriparatide",
            "value",
        ]
    ].copy()
    observation_diagnostics.insert(0, "endpoint", endpoint)
    observation_diagnostics["fitted"] = classical_model.fittedvalues
    observation_diagnostics["residual"] = classical_model.resid
    observation_diagnostics["studentized_residual"] = studentized
    observation_diagnostics["cooks_distance"] = cooks
    observation_diagnostics["outlier_abs_studentized_gt_3"] = (
        np.abs(studentized) > 3
    )
    observation_diagnostics["influential_cooks_gt_4_over_n"] = (
        cooks > 4 / len(data)
    )

    plot_dir = output_dir / "diagnostic_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    axes[0].scatter(
        classical_model.fittedvalues,
        classical_model.resid,
        color="#1f77b4",
        edgecolor="black",
        linewidth=0.4,
    )
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_xlabel("Fitted value")
    axes[0].set_ylabel("Residual")
    axes[0].set_title("Residuals vs fitted")
    sm.qqplot(
        classical_model.resid,
        line="45",
        ax=axes[1],
        markerfacecolor="#1f77b4",
        markeredgecolor="black",
    )
    axes[1].set_title("Normal Q-Q")
    figure.suptitle(endpoint)
    figure.tight_layout()
    plot_buffer = io.BytesIO()
    figure.savefig(
        plot_buffer,
        format="png",
        dpi=200,
        bbox_inches="tight",
    )
    (
        plot_dir / f"{safe_name(endpoint)}_diagnostics.png"
    ).write_bytes(plot_buffer.getvalue())
    plt.close(figure)

    raw_cell_means = (
        data.groupby(
            ["concentration", "teriparatide"],
            observed=True,
        )["value"]
        .mean()
    )
    rescue_rows = contrasts_df.loc[
        contrasts_df["family"] == "rescue_MVT"
    ].reset_index(drop=True)
    diagonal_rows = contrasts_df.loc[
        contrasts_df["family"] == "exploratory_diagonal_MVT"
    ].reset_index(drop=True)
    classifications = []
    for index, concentration in enumerate(nonzero):
        baseline = float(raw_cell_means.loc[("0", "No")])
        damaged = float(raw_cell_means.loc[(concentration, "No")])
        treated = float(raw_cell_means.loc[(concentration, "Yes")])
        impairment = damaged - baseline
        rescue_change = treated - damaged
        toward_baseline = impairment * rescue_change < 0
        rescue_significant = rescue_rows.loc[index, "p_adjusted"] < 0.05
        diagonal_significant = diagonal_rows.loc[index, "p_adjusted"] < 0.05
        remaining_difference = treated - baseline

        if rescue_significant and toward_baseline:
            if diagonal_significant and impairment * remaining_difference > 0:
                category = "partial rescue"
            elif diagonal_significant and impairment * remaining_difference < 0:
                category = "improvement with overshoot"
            else:
                category = (
                    "improvement; complete rescue not established"
                )
        else:
            category = "no adjusted evidence of improvement"

        classifications.append(
            {
                "endpoint": endpoint,
                "concentration_um": int(concentration),
                "baseline_mean_0uM_no_teriparatide": baseline,
                "4HC_alone_mean": damaged,
                "combined_treatment_mean": treated,
                "4HC_impairment_direction": np.sign(impairment),
                "rescue_change_toward_baseline": toward_baseline,
                "rescue_MVT_p": rescue_rows.loc[index, "p_adjusted"],
                "diagonal_MVT_p": diagonal_rows.loc[index, "p_adjusted"],
                "classification": category,
                "complete_rescue_statement": (
                    "Not established: no prespecified equivalence margin or "
                    "equivalence test; p>0.05 versus baseline is insufficient."
                ),
            }
        )
    classification_df = pd.DataFrame(classifications)

    interaction_row = anova.loc[
        anova["term"] == "4-HC concentration x teriparatide"
    ].iloc[0]
    summary = {
        "endpoint": endpoint,
        "source_file": path.name,
        "source_value_column": source_value_column,
        "n_independent_observations": len(data),
        "n_concentrations": len(concentrations),
        "n_per_cell_min": int(
            data.groupby(
                ["concentration", "teriparatide"], observed=True
            ).size().min()
        ),
        "n_per_cell_max": int(
            data.groupby(
                ["concentration", "teriparatide"], observed=True
            ).size().max()
        ),
        "technical_replicate_handling": (
            "None; each source row is an independent experimental observation."
        ),
        "covariance_method": covariance_method,
        "brown_forsythe_p": brown_forsythe.pvalue,
        "interaction_p": interaction_row["p_value"],
        "residual_shapiro_p": shapiro.pvalue,
        "outlier_count_abs_studentized_gt_3": int(
            np.sum(np.abs(studentized) > 3)
        ),
        "influential_count_cooks_gt_4_over_n": int(
            np.sum(cooks > 4 / len(data))
        ),
    }
    return {
        "anova": anova,
        "contrasts": contrasts_df,
        "descriptions": descriptions,
        "diagnostics": diagnostics,
        "observation_diagnostics": observation_diagnostics,
        "classifications": classification_df,
        "summary": summary,
    }


def write_readme(
    output_dir: Path,
    input_dir: Path,
    summaries: pd.DataFrame,
) -> None:
    hc3_count = int((summaries["covariance_method"] == "HC3").sum())
    text = f"""# Two-way ANOVA of in vitro endpoints

## Scope

- Input: `{display_path(input_dir)}`
- Endpoints: {len(summaries)}
- Analysis unit: each CSV row is one independent experimental observation.
- No technical-replicate averaging or batch pairing was applied.
- The source CSV files were not modified.

## Model

For every endpoint:

`value ~ categorical 4-HC concentration * teriparatide`

Type III tests use sum-to-zero contrasts. Brown-Forsythe testing across
factorial cells was used as a prespecified covariance rule: classical
covariance when p >= 0.05 and HC3 heteroskedasticity-robust covariance when
p < 0.05. HC3 was selected for {hc3_count} of {len(summaries)} endpoints.

## Planned comparisons

1. Teriparatide versus no teriparatide at 0 uM: no multiplicity correction.
2. Each nonzero 4-HC concentration versus 0 uM without teriparatide:
   Dunnett single-step multivariate-t correction.
3. Teriparatide versus no teriparatide at each nonzero concentration:
   one rescue family with single-step multivariate-t correction.
4. Each combined-treatment cell versus 0 uM/no-teriparatide control:
   a separate exploratory family with single-step multivariate-t correction.

The multivariate-t calculations use reproducible scrambled Sobol
quasi-Monte Carlo integration (`seed={SEED}`, `2^{QMC_POWER}` draws per
family). Simultaneous 95% confidence intervals are reported for adjusted
families.

## Rescue interpretation

- `partial rescue`: adjusted rescue comparison is significant in the
  direction of baseline and the adjusted diagonal comparison still differs
  from baseline in the original impairment direction.
- `improvement; complete rescue not established`: rescue is significant
  toward baseline, while the diagonal comparison is not significant.
- `complete rescue` is never concluded here because no equivalence margin
  was prespecified. A nonsignificant difference from baseline does not
  establish equivalence.

## Source column labels

Some source files use the outcome column label `OD450` despite filenames
indicating scratch or tube-network metrics. Results retain both the filename
endpoint and original source column so these labels can be checked before
manuscript use.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def write_endpoint_comparison_summary(
    output_dir: Path,
    summary_df: pd.DataFrame,
    contrasts_df: pd.DataFrame,
    classifications_df: pd.DataFrame,
) -> None:
    rows = []
    for _, summary in summary_df.iterrows():
        endpoint = summary["endpoint"]
        endpoint_contrasts = contrasts_df.loc[
            contrasts_df["endpoint"] == endpoint
        ]
        endpoint_classes = classifications_df.loc[
            classifications_df["endpoint"] == endpoint
        ]

        def family_text(family: str) -> str:
            subset = endpoint_contrasts.loc[
                endpoint_contrasts["family"] == family
            ]
            return " | ".join(
                f"{row.comparison}: estimate={row.estimate:.6g}, "
                f"adjusted p={row.p_adjusted:.6g}"
                for row in subset.itertuples()
            )

        rows.append(
            {
                "endpoint": endpoint,
                "interaction_p": summary["interaction_p"],
                "covariance_method": summary["covariance_method"],
                "teriparatide_alone_at_0uM": family_text(
                    "teriparatide_alone_at_0uM"
                ),
                "4HC_alone_Dunnett": family_text("4HC_alone_Dunnett"),
                "rescue_MVT": family_text("rescue_MVT"),
                "exploratory_diagonal_MVT": family_text(
                    "exploratory_diagonal_MVT"
                ),
                "rescue_interpretation": " | ".join(
                    f"{int(row.concentration_um)} uM: {row.classification}"
                    for row in endpoint_classes.itertuples()
                ),
            }
        )
    pd.DataFrame(rows).to_csv(
        output_dir / "endpoint_planned_comparison_summary.csv",
        index=False,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(input_dir.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {input_dir}")

    all_results: dict[str, list[pd.DataFrame]] = {
        "anova": [],
        "contrasts": [],
        "descriptions": [],
        "diagnostics": [],
        "observation_diagnostics": [],
        "classifications": [],
    }
    summaries = []
    manifest = []

    for endpoint_index, path in enumerate(files):
        result = analyse_endpoint(path, output_dir, endpoint_index)
        endpoint_dir = output_dir / "per_endpoint" / safe_name(path.stem)
        endpoint_dir.mkdir(parents=True, exist_ok=True)
        for key in all_results:
            frame = result[key]
            assert isinstance(frame, pd.DataFrame)
            frame.to_csv(endpoint_dir / f"{key}.csv", index=False)
            all_results[key].append(frame)
        summary = result["summary"]
        assert isinstance(summary, dict)
        summaries.append(summary)
        manifest.append(
            {
                "source_file": path.name,
                "source_path": display_path(path),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
        )
        print(f"Analysed {path.name}")

    combined_dir = output_dir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    for key, frames in all_results.items():
        pd.concat(frames, ignore_index=True).to_csv(
            combined_dir / f"all_endpoints_{key}.csv",
            index=False,
        )

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(output_dir / "endpoint_audit_summary.csv", index=False)
    combined_contrasts = pd.concat(
        all_results["contrasts"], ignore_index=True
    )
    combined_classifications = pd.concat(
        all_results["classifications"], ignore_index=True
    )
    write_endpoint_comparison_summary(
        output_dir,
        summary_df,
        combined_contrasts,
        combined_classifications,
    )
    pd.DataFrame(manifest).to_csv(
        output_dir / "source_sha256_manifest.csv", index=False
    )
    write_readme(output_dir, input_dir, summary_df)

    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
        "seed": SEED,
        "qmc_power": QMC_POWER,
        "input_dir": display_path(input_dir),
        "output_dir": display_path(output_dir),
    }
    (output_dir / "environment.json").write_text(
        json.dumps(environment, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"Wrote results to {output_dir}")


if __name__ == "__main__":
    main()
