# Derived analysis data

This directory contains anonymized derived numerical outputs used to summarize and visualize the reported analyses. Raw BOLD and structural imaging files are not redistributed here; they should be obtained from the source dataset described in the manuscript and repository README.

Final ROI NIfTI masks are distributed separately under `resources/rois/`.

Key participant-level files include:

- `derived/dcm/DCM_A_allsubjects_wide.csv` — primary full-run A-matrix export (45 participants × 2 runs × 4 analytical levels).
- `derived/fc_raw/FC_allsubjects_wide.csv` — raw functional-connectivity export.
- `derived/fc_residual/ResidualFC_allsubjects_wide.csv` — task-regressed residual functional-connectivity export.
- `derived/peb/` — PEB/BMA summary exports.
- `derived/dcm_model_fit/` and `derived/dcm_qc/` — model-fit and posterior-QC outputs.
- `derived/block_timing/` — three-block timing table.
