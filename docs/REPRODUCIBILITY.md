# Reproducibility notes

The repository distributes analysis code, final ROI masks, anonymized derived numerical outputs, figure source data, and sensitivity-analysis materials. Raw participant imaging files are not redistributed.

## External inputs for complete first-level reprocessing

Re-running the distributed first-level analysis requires:

1. the source task-fMRI dataset reported by Duncan et al. (2009), prepared with the preprocessing workflow described in the manuscript;
2. the corresponding motion, outlier, and tissue-mask files expected by the analysis scripts;
3. MATLAB with SPM12 available on the path;
4. a writable processing directory.

CONN was used in the upstream preprocessing workflow described in the manuscript. The code distributed here begins from the preprocessed imaging and nuisance-regressor files used by the reported GLM/DCM and FC analyses.

The four ROI NIfTI masks used in the analyses are included under `resources/rois/`. They can also be regenerated using `code/matlab/roi/make_analysis_roi_masks.m`.

## Path configuration

Machine-specific absolute paths are not stored in the repository. External locations are configured through `code/matlab/project_paths.m` and environment variables described in the top-level README.

## Primary and sensitivity analyses

The primary DCM analysis uses block-specific GLM/input definitions while retaining the complete run time series. The physically cropped analysis under `sensitivity/cropped_dcm/` is distributed as an auxiliary sensitivity analysis and includes its participant-level DCM export and QA tables.

## Included verification material

- ROI voxelwise regeneration QC: `resources/rois/qc/`
- DCM model-fit and posterior QC: `data/derived/dcm_model_fit/` and `data/derived/dcm_qc/`
- Figure source data: `figures/source_data/`
- Repository file hashes: `FILE_MANIFEST.csv`
