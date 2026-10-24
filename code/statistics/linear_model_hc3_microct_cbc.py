"""Micro-CT and complete blood count endpoints (Fig. 1B, 1D, Fig. S1).

2 x 2 factorial linear model (Cy x iPTH) with HC3 heteroskedasticity-robust
covariance. Four prespecified contrasts are adjusted jointly by the
multivariate-t method; Cy+iPTH vs Control is reported as exploratory.

Usage (from the repository root):
    python code/statistics/linear_model_hc3_microct_cbc.py

Input:  data/statistics/microct_cbc/*.csv (one column per group)
Output: results/statistics/microct_cbc_hc3/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import build_design_matrices
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]

GROUP_MAP = {
    "Control": ("No", "No"),
    "Cyclophosphamide": ("Yes", "No"),
    "Cyclophosphamide + Teriparatide": ("Yes", "Yes"),
    "Teriparatide": ("No", "Yes"),
}

# Cell order: Control, Cy, iPTH, Cy+iPTH.
PRIMARY_CONTRASTS = {
    "Cy - Control": np.array([-1.0, 1.0, 0.0, 0.0]),
    "iPTH - Control": np.array([-1.0, 0.0, 1.0, 0.0]),
    "Cy+iPTH - Cy": np.array([0.0, -1.0, 0.0, 1.0]),
    "Cy+iPTH - iPTH": np.array([0.0, 0.0, -1.0, 1.0]),
}
EXPLORATORY_CONTRAST = {
    "Cy+iPTH - Control": np.array([-1.0, 0.0, 0.0, 1.0])
}


def significance_symbol(p_value: float) -> str:
    if p_value < 0.0001:
        return "#"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def prepare_data(path: Path) -> pd.DataFrame:
    wide = pd.read_csv(path)
    missing = set(GROUP_MAP).difference(wide.columns)
    unexpected = set(wide.columns).difference(GROUP_MAP)
    if missing or unexpected:
        raise ValueError(
            f"{path.name}: expected exactly {list(GROUP_MAP)}; "
            f"missing={sorted(missing)}, unexpected={sorted(unexpected)}"
        )

    data = wide.melt(var_name="group", value_name="value").dropna()
    data["value"] = pd.to_numeric(data["value"], errors="raise")
    factors = data["group"].map(GROUP_MAP)
    data["cy"] = [item[0] for item in factors]
    data["ipth"] = [item[1] for item in factors]
    data["cy"] = pd.Categorical(data["cy"], ["No", "Yes"])
    data["ipth"] = pd.Categorical(data["ipth"], ["No", "Yes"])
    return data


def mvt_adjusted_p_values(
    t_values: np.ndarray,
    correlation: np.ndarray,
    df_residual: float,
    seed: int,
) -> np.ndarray:
    correlation = (correlation + correlation.T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    keep = eigenvalues > 1e-10
    loading = eigenvectors[:, keep] @ np.diag(np.sqrt(eigenvalues[keep]))

    # The four planned contrasts span a rank-three contrast space. A
    # scrambled Sobol sequence integrates the resulting singular
    # multivariate-t distribution reproducibly without forcing it into a
    # nonsingular four-dimensional rectangle.
    power = 20
    uniforms = stats.qmc.Sobol(
        d=loading.shape[1] + 1,
        scramble=True,
        seed=seed,
    ).random_base2(power)
    eps = np.finfo(float).eps
    uniforms = np.clip(uniforms, eps, 1 - eps)
    independent_normals = stats.norm.ppf(uniforms[:, :-1])
    chi_square = stats.chi2.ppf(uniforms[:, -1], df_residual)
    correlated_normals = independent_normals @ loading.T
    simulated_t = correlated_normals / np.sqrt(
        chi_square[:, None] / df_residual
    )
    max_distribution = np.max(np.abs(simulated_t), axis=1)
    adjusted = np.asarray(
        [np.mean(max_distribution >= abs(t_value)) for t_value in t_values]
    )
    raw = 2 * stats.t.sf(np.abs(t_values), df_residual)
    bonferroni_upper_bound = np.minimum(1.0, len(t_values) * raw)
    tolerance = 5 / len(max_distribution)
    if np.any(adjusted > bonferroni_upper_bound + tolerance):
        raise RuntimeError("MVT integration failed the Bonferroni upper-bound check.")
    return adjusted


def analyse_endpoint(
    path: Path,
    output_dir: Path,
    seed: int,
    covariance: str,
) -> pd.DataFrame:
    data = prepare_data(path)
    classical_model = smf.ols(
        "value ~ C(cy, Sum) * C(ipth, Sum)", data=data
    ).fit()
    if covariance == "hc3":
        model = classical_model.get_robustcov_results(cov_type="HC3", use_t=True)
        anova = sm.stats.anova_lm(
            classical_model, typ=3, robust="hc3"
        ).reset_index()
    else:
        model = classical_model
        anova = sm.stats.anova_lm(classical_model, typ=3).reset_index()
    anova = anova.rename(
        columns={
            "index": "term",
            "sum_sq": "sum_squares",
            "df": "df",
            "F": "f_value",
            "PR(>F)": "p_value",
        }
    )
    term_names = {
        "C(cy, Sum)": "Cy main effect",
        "C(ipth, Sum)": "iPTH main effect",
        "C(cy, Sum):C(ipth, Sum)": "Cy x iPTH interaction",
        "Residual": "Residual",
        "Intercept": "Intercept",
    }
    anova["term"] = anova["term"].replace(term_names)
    anova["endpoint"] = path.stem
    anova["covariance"] = covariance.upper()
    anova.to_csv(
        output_dir / f"{path.stem}_two_way_ANOVA_{covariance.upper()}.csv",
        index=False,
    )

    cell_order = [
        ("No", "No"),
        ("Yes", "No"),
        ("No", "Yes"),
        ("Yes", "Yes"),
    ]
    grid = pd.DataFrame(cell_order, columns=["cy", "ipth"])
    grid["cy"] = pd.Categorical(grid["cy"], ["No", "Yes"])
    grid["ipth"] = pd.Categorical(grid["ipth"], ["No", "Yes"])
    design = np.asarray(
        build_design_matrices([classical_model.model.data.design_info], grid)[0]
    )
    beta = np.asarray(model.params)
    beta_cov = np.asarray(model.cov_params())
    cell_means = design @ beta
    cell_cov = design @ beta_cov @ design.T

    primary_matrix = np.vstack(list(PRIMARY_CONTRASTS.values()))
    estimates = primary_matrix @ cell_means
    contrast_cov = primary_matrix @ cell_cov @ primary_matrix.T
    standard_errors = np.sqrt(np.diag(contrast_cov))
    t_values = estimates / standard_errors
    raw_p_values = 2 * stats.t.sf(np.abs(t_values), model.df_resid)
    correlation = contrast_cov / np.outer(standard_errors, standard_errors)
    adjusted_p_values = mvt_adjusted_p_values(
        t_values,
        correlation,
        model.df_resid,
        seed,
    )

    rows = []
    for index, comparison in enumerate(PRIMARY_CONTRASTS):
        rows.append(
            {
                "endpoint": path.stem,
                "comparison_family": "four_primary",
                "comparison": comparison,
                "estimate": estimates[index],
                "std_error": standard_errors[index],
                "df": model.df_resid,
                "t_value": t_values[index],
                "p_value_raw": raw_p_values[index],
                "p_value_reported": adjusted_p_values[index],
                "multiplicity_method": (
                    "multivariate-t across four prespecified comparisons; "
                    "Sobol quasi-Monte Carlo integration"
                ),
                "significance": significance_symbol(adjusted_p_values[index]),
            }
        )

    exploratory_name, exploratory_vector = next(iter(EXPLORATORY_CONTRAST.items()))
    exploratory_estimate = float(exploratory_vector @ cell_means)
    exploratory_variance = float(
        exploratory_vector @ cell_cov @ exploratory_vector
    )
    exploratory_se = np.sqrt(exploratory_variance)
    exploratory_t = exploratory_estimate / exploratory_se
    exploratory_p = float(2 * stats.t.sf(abs(exploratory_t), model.df_resid))
    rows.append(
        {
            "endpoint": path.stem,
            "comparison_family": "exploratory",
            "comparison": exploratory_name,
            "estimate": exploratory_estimate,
            "std_error": exploratory_se,
            "df": model.df_resid,
            "t_value": exploratory_t,
            "p_value_raw": exploratory_p,
            "p_value_reported": exploratory_p,
            "multiplicity_method": "none; exploratory nominal p value",
            "significance": significance_symbol(exploratory_p),
        }
    )

    contrasts = pd.DataFrame(rows)
    contrasts["ci_95_low"] = contrasts["estimate"] - stats.t.ppf(
        0.975, model.df_resid
    ) * contrasts["std_error"]
    contrasts["ci_95_high"] = contrasts["estimate"] + stats.t.ppf(
        0.975, model.df_resid
    ) * contrasts["std_error"]
    contrasts["model"] = "value ~ Cy * iPTH"
    contrasts["covariance"] = covariance.upper()
    contrasts["symbol_definition"] = (
        "*=p<0.05; **=p<0.01; ***=p<0.001; #=p<0.0001; ns=p>=0.05"
    )
    contrasts.to_csv(
        output_dir
        / (
            f"{path.stem}_{covariance.upper()}_"
            "primary4_MVT_exploratory_contrasts.csv"
        ),
        index=False,
    )

    descriptions = (
        data.groupby("group", observed=True)["value"]
        .agg(n="count", mean="mean", sd="std", median="median", minimum="min", maximum="max")
        .reset_index()
    )
    descriptions["endpoint"] = path.stem
    descriptions.to_csv(
        output_dir / f"{path.stem}_descriptive_statistics.csv", index=False
    )

    groups = [
        data.loc[data["group"] == group, "value"].to_numpy()
        for group in GROUP_MAP
    ]
    shapiro = stats.shapiro(model.resid)
    levene = stats.levene(*groups, center="median")
    diagnostics = pd.DataFrame(
        [
            {
                "endpoint": path.stem,
                "test": "Shapiro-Wilk test of model residuals",
                "statistic": shapiro.statistic,
                "p_value": shapiro.pvalue,
            },
            {
                "endpoint": path.stem,
                "test": "Levene test across four groups (median centered)",
                "statistic": levene.statistic,
                "p_value": levene.pvalue,
            },
        ]
    )
    diagnostics.to_csv(
        output_dir / f"{path.stem}_model_diagnostics.csv", index=False
    )
    return contrasts


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=REPO_ROOT / "data" / "statistics" / "microct_cbc",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="default: results/statistics/microct_cbc_<covariance>",
    )
    parser.add_argument("--seed", type=int, default=20260615)
    parser.add_argument(
        "--covariance",
        choices=("classical", "hc3"),
        default="hc3",
    )
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = (
            REPO_ROOT / "results" / "statistics" / f"microct_cbc_{args.covariance}"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.input_dir.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {args.input_dir}")

    results = []
    for path in files:
        results.append(
            analyse_endpoint(path, args.output_dir, args.seed, args.covariance)
        )
        print(f"Analysed: {path.name}")

    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(
        args.output_dir
        / (
            f"all_endpoints_{args.covariance.upper()}_"
            "primary4_MVT_exploratory_contrasts.csv"
        ),
        index=False,
    )
    print(f"Wrote results to: {args.output_dir}")


if __name__ == "__main__":
    main()
