"""Assemble the tables used for the significance labels in the figures.

Run after linear_model_hc3_microct_cbc.py and two_way_anova_in_vitro.py.
In addition to the planned comparisons of two_way_anova_in_vitro.py, this adds
one family for the in vitro endpoints: each nonzero 4-HC concentration vs
0 uM in the presence of iPTH (Dunnett single-step multivariate-t).

Usage (from the repository root):
    python code/statistics/build_figure_significance.py

Output: results/statistics/figure_significance/
    microCT_CBC_nonLMM_HC3_primary4_MVT_contrasts.csv   (Fig. 1B/1D, Fig. S1)
    vitro_all_endpoints_anova.csv                       (Fig. 5)
    vitro_all_endpoints_contrasts.csv                   (Fig. 5)
    vitro_endpoint_planned_comparison_summary.csv       (Fig. 5)
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import build_design_matrices
from scipy import stats

from two_way_anova_in_vitro import (
    DEFAULT_INPUT,
    REPO_ROOT,
    SEED,
    build_cell_grid,
    contrast_family,
    prepare_data,
    vector_for,
)

RESULTS = REPO_ROOT / "results" / "statistics"
FAMILY = "4HC_under_teriparatide_Dunnett"


def dunnett_with_ipth(path: Path, endpoint_index: int) -> pd.DataFrame:
    """Each nonzero 4-HC concentration vs 0 uM, all with iPTH."""
    data, _ = prepare_data(path)
    classical_model = smf.ols(
        "value ~ C(concentration, Sum) * C(teriparatide, Sum)", data=data
    ).fit()
    cell_groups = [
        group["value"].to_numpy()
        for _, group in data.groupby(
            ["concentration", "teriparatide"], observed=True, sort=False
        )
    ]
    # Same covariance rule as the main analysis (Brown-Forsythe across cells).
    use_hc3 = bool(stats.levene(*cell_groups, center="median").pvalue < 0.05)
    model = (
        classical_model.get_robustcov_results(cov_type="HC3", use_t=True)
        if use_hc3
        else classical_model
    )

    grid = build_cell_grid(data)
    design = np.asarray(
        build_design_matrices([classical_model.model.data.design_info], grid)[0]
    )
    beta = np.asarray(model.params)
    cell_means = design @ beta
    cell_covariance = design @ np.asarray(model.cov_params()) @ design.T
    cell_lookup = {
        (str(row.concentration), str(row.teriparatide)): index
        for index, row in grid.iterrows()
    }

    nonzero = [str(c) for c in data["concentration"].cat.categories if str(c) != "0"]
    rows = contrast_family(
        endpoint=path.stem,
        family=FAMILY,
        labels=[f"{c} uM - 0 uM, with teriparatide" for c in nonzero],
        contrast_matrix=np.vstack(
            [vector_for(cell_lookup, (c, "Yes"), ("0", "Yes")) for c in nonzero]
        ),
        cell_means=cell_means,
        cell_covariance=cell_covariance,
        df_residual=model.df_resid,
        adjustment="Dunnett single-step multivariate-t",
        seed=SEED + endpoint_index * 10 + 4,
    )
    rows["covariance_method"] = "HC3" if use_hc3 else "classical"
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--in-vitro-input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--microct-results", type=Path, default=RESULTS / "microct_cbc_hc3")
    parser.add_argument("--in-vitro-results", type=Path, default=RESULTS / "in_vitro_two_way_anova")
    parser.add_argument("--output-dir", type=Path, default=RESULTS / "figure_significance")
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    shutil.copyfile(
        args.microct_results / "all_endpoints_HC3_primary4_MVT_exploratory_contrasts.csv",
        out / "microCT_CBC_nonLMM_HC3_primary4_MVT_contrasts.csv",
    )
    shutil.copyfile(
        args.in_vitro_results / "combined" / "all_endpoints_anova.csv",
        out / "vitro_all_endpoints_anova.csv",
    )
    shutil.copyfile(
        args.in_vitro_results / "endpoint_planned_comparison_summary.csv",
        out / "vitro_endpoint_planned_comparison_summary.csv",
    )

    # Planned comparisons from the main analysis, followed by the extra family.
    contrasts_path = out / "vitro_all_endpoints_contrasts.csv"
    shutil.copyfile(
        args.in_vitro_results / "combined" / "all_endpoints_contrasts.csv", contrasts_path
    )
    columns = pd.read_csv(contrasts_path, nrows=0).columns
    extra = pd.concat(
        [
            dunnett_with_ipth(path, index)
            for index, path in enumerate(sorted(args.in_vitro_input.glob("*.csv")))
        ],
        ignore_index=True,
    )
    extra.loc[:, columns].to_csv(contrasts_path, mode="a", index=False, header=False)
    print(f"Wrote {out} ({len(extra)} additional contrast rows)")


if __name__ == "__main__":
    main()
