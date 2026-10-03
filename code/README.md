# Analysis code

## MATLAB pipeline

The MATLAB scripts implement ROI reconstruction, block timing, first-level GLM specification, DCM estimation, DCM export and QC, PEB/BMR/BMA group analysis, raw FC, and task-regressed residual FC.

Recommended execution order:

1. `roi/make_analysis_roi_masks.m`
2. `block_timing/build_three_block_timing_table.m`
3. `glm_dcm/run_fullrun_block_specific_glm.m`
4. `glm_dcm/build_fullrun_block_specific_dcm.m`
5. `dcm_export/export_dcm_a_matrices.m`
6. `dcm_export/summarize_dcm_a_matrices.m`
7. `qc/run_dcm_parameter_qc.m`
8. `qc/run_dcm_model_fit_qc.m`
9. `peb/run_peb_group_analysis.m`
10. `fc/run_raw_fc_analysis.m`
11. `residual_fc/run_residual_fc_analysis.m`

External paths are configured through `project_paths.m`; no machine-specific absolute paths are embedded in the scripts.

## Python summary script

`python/summarize_primary_results.py` reproduces the participant-level run-averaged Late contrasts for the 12 directed between-region DCM parameters and the six raw/residual FC edges from the distributed CSV exports. It applies FDR correction separately within each analysis family.
