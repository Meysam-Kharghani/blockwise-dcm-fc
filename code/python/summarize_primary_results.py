#!/usr/bin/env python3
"""Reproduce the primary run-averaged Late contrasts from distributed CSV exports."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

BLOCKS = ["Block-01", "Block-02", "Block-03"]


def paired_late_contrasts(
    df: pd.DataFrame,
    subject_col: str,
    run_col: str,
    condition_col: str,
    value_cols: list[str],
) -> pd.DataFrame:
    rows = []
    for value_col in value_cols:
        tmp = df[[subject_col, run_col, condition_col, value_col]].copy()
        tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce")

        valid_rows = tmp.dropna(subset=[value_col])
        run_counts = (
            valid_rows.groupby([subject_col, condition_col])[run_col]
            .nunique()
            .unstack(fill_value=0)
        )
        for block in BLOCKS:
            if block not in run_counts.columns:
                run_counts[block] = 0
        complete_subjects = run_counts.index[(run_counts[BLOCKS] == 2).all(axis=1)]
        tmp = tmp[tmp[subject_col].isin(complete_subjects)]
        tmp = (
            tmp.groupby([subject_col, condition_col], as_index=False)[value_col]
            .mean()
            .pivot(index=subject_col, columns=condition_col, values=value_col)
        )
        for block in BLOCKS:
            if block not in tmp.columns:
                tmp[block] = np.nan
        late = tmp["Block-03"] - tmp[["Block-01", "Block-02"]].mean(axis=1)
        late = late.dropna()
        n = int(late.size)
        mean_diff = float(late.mean()) if n else np.nan
        sd_diff = float(late.std(ddof=1)) if n > 1 else np.nan
        dz = mean_diff / sd_diff if n > 1 and sd_diff > 0 else np.nan
        t_stat, p_value = stats.ttest_1samp(late.to_numpy(), 0.0, nan_policy="omit") if n > 1 else (np.nan, np.nan)
        rows.append(
            {
                "parameter": value_col,
                "N": n,
                "mean_late_contrast": mean_diff,
                "sd_late_contrast": sd_diff,
                "cohen_dz": dz,
                "t": float(t_stat),
                "p": float(p_value),
            }
        )

    out = pd.DataFrame(rows)
    valid = np.isfinite(out["p"].to_numpy(dtype=float))
    q = np.full(len(out), np.nan)
    if valid.any():
        q[valid] = fdrcorrection(out.loc[valid, "p"].to_numpy(dtype=float))[1]
    out["q_FDR"] = q
    out["FDR_sig_05"] = out["q_FDR"] < 0.05
    return out


def dcm_offdiagonal_columns(df: pd.DataFrame) -> list[str]:
    cols = []
    for col in df.columns:
        if not col.startswith("A_") or "_to_" not in col:
            continue
        source, target = col[2:].split("_to_", 1)
        if source != target:
            cols.append(col)
    return cols


def fc_edge_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if "__" in col]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data" / "derived",
        help="Path to the repository data/derived directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "results" / "reproduced",
        help="Directory for regenerated summary tables.",
    )
    args = parser.parse_args()

    data_root = args.data_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    dcm = pd.read_csv(data_root / "dcm" / "DCM_A_allsubjects_wide.csv")
    dcm = dcm[dcm["Condition"].isin(BLOCKS)].copy()
    dcm_results = paired_late_contrasts(
        dcm,
        subject_col="Subject",
        run_col="Run",
        condition_col="Condition",
        value_cols=dcm_offdiagonal_columns(dcm),
    )
    dcm_results.to_csv(output_dir / "dcm_late_contrasts_runaveraged.csv", index=False)

    for label, relative_path in {
        "raw_fc": Path("fc_raw") / "FC_allsubjects_wide.csv",
        "residual_fc": Path("fc_residual") / "ResidualFC_allsubjects_wide.csv",
    }.items():
        fc = pd.read_csv(data_root / relative_path)
        fc = fc[fc["tag"].isin(BLOCKS)].copy()
        fc_results = paired_late_contrasts(
            fc,
            subject_col="subject",
            run_col="run",
            condition_col="tag",
            value_cols=fc_edge_columns(fc),
        )
        fc_results.to_csv(output_dir / f"{label}_late_contrasts_runaveraged.csv", index=False)

    print(f"Wrote regenerated summaries to {output_dir}")


if __name__ == "__main__":
    main()
