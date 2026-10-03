# Block-Specific Full-Run DCM and Functional Connectivity Analysis

This repository contains analysis code, final ROI masks, anonymized derived outputs, figure source data, quality-control summaries, and a cropped-segment sensitivity analysis accompanying the study **“Block-Specific Full-Run Effective Connectivity Estimates in a Visual Occipitotemporal-Frontal Model During a Visual One-Back Task.”**

## Authors

1. **Amir Hakimjavadi** — Department of Psychology, Universitat Rovira i Virgili, Tarragona, Spain
2. **Meysam Kharghani** — Department of Cognitive Neuroscience, Faculty of Education and Psychology, University of Tabriz, Tabriz, Iran
3. **Arash Zare-Sadeghi** — Medical Physics Department, Iran University of Medical Sciences, Tehran, Iran — corresponding author: zare.a@iums.ac.ir

Full contact information is provided in `AUTHORS.md`.

## Study overview

The project is a secondary analysis of the open-access visual one-back task-fMRI dataset reported by Duncan et al. (2009). The modeled left-hemisphere network contains V1_L, lateral occipitotemporal cortex (lOTC), ventral occipitotemporal cortex (vOTC), and left inferior frontal gyrus (IFG_L).

The primary effective-connectivity analysis uses separate block-specific GLMs and DCMs estimated from the complete run time series. For Block-01, Block-02, and Block-03, the target analytical block is represented explicitly while scans outside that block remain in the acquired run. The primary inferential target is the DCM A-matrix. Group-level inference uses Parametric Empirical Bayes (PEB), Bayesian model reduction (BMR), and Bayesian model averaging (BMA).

Complementary analyses include raw functional connectivity (FC), task-regressed residual FC, DCM model-fit and posterior summaries, and a physically cropped DCM sensitivity analysis.

## Repository structure

```text
blockwise-dcm-fc/
├── code/
│   ├── matlab/                 # GLM/DCM, PEB, FC, QC, ROI reconstruction
│   └── python/                 # verification summaries from distributed CSV exports
├── data/derived/               # anonymized numerical exports
├── results/tables_csv/         # analysis summary tables
├── figures/
│   ├── main/                   # manuscript main figures
│   ├── supplementary/          # supplementary and diagnostic figures
│   └── source_data/            # tabular figure source data
├── resources/rois/             # final NIfTI masks and ROI QC
├── sensitivity/cropped_dcm/    # physically cropped DCM sensitivity analysis
├── docs/                       # reproducibility and methods notes
├── AUTHORS.md
├── CITATION.cff
├── FILE_MANIFEST.csv
├── LICENSE
└── VERSION
```

## Final ROI masks

The exact masks used for VOI extraction are distributed in `resources/rois/`.

| ROI | Definition | Nonzero voxels | Approx. volume |
|---|---|---:|---:|
| IFG_L | 6-mm sphere centered at MNI `[-48, 28, 0]` | 123 | 984 mm³ |
| V1_L | 6-mm sphere centered at MNI `[-16, -98, -10]` | 123 | 984 mm³ |
| lOTC | X `-58:-36`, Y `-86:-66`, Z `-18:6` | 1716 | 13,728 mm³ |
| vOTC | X `-56:-32`, Y `-68:-44`, Z `-28:-2` | 2366 | 18,928 mm³ |

Voxelwise reconstruction QC produced exact matches for all four masks. lOTC and vOTC overlap by 198 voxels; all other ROI pairs have zero overlap. Additional mask metadata and QC are provided in `resources/rois/README.md` and `resources/rois/qc/`.

## Primary DCM configuration

- 45 participants × 2 runs × 4 analytical levels = 360 primary DCMs.
- Block-wise PEB uses 270 DCMs (45 participants × 2 runs × 3 blocks).
- Echo time used for DCM specification: **0.05 s (50 ms)**.
- All 12 directed off-diagonal A-matrix connections are enabled.
- Standard SPM self-inhibition parameters are estimated.
- No condition-specific B-matrix modulation is included.
- Task-driving input enters all four modeled regions through the C-matrix.

The participant-level primary A-matrix export is `data/derived/dcm/DCM_A_allsubjects_wide.csv`.

## Software and path configuration

The distributed MATLAB analysis scripts require MATLAB and SPM12. CONN was used in the upstream preprocessing workflow described in the manuscript; the distributed analysis scripts operate on the resulting preprocessed images and associated nuisance files. The Python verification script uses NumPy, pandas, SciPy, and statsmodels.

```bash
pip install -r requirements.txt
```

MATLAB path configuration is handled by `code/matlab/project_paths.m`. Optional environment variables are:

- `DCM_FC_DATA_ROOT` — root of the external source dataset.
- `DCM_FC_WORK_ROOT` — writable processing directory; defaults to `<repository>/work`.
- `DCM_FC_ROI_ROOT` — ROI directory; defaults to `<repository>/resources/rois`.
- `SPM12_DIR` — SPM12 installation directory.
- `CONN_DIR` — CONN installation directory.

No machine-specific absolute paths are stored in this repository.

From MATLAB, add `code/matlab/` to the MATLAB path and run:

```matlab
cfg = setup_repository();
```

## MATLAB execution order

1. `roi/make_analysis_roi_masks.m` — reconstruct the distributed ROI masks if required.
2. `block_timing/build_three_block_timing_table.m` — construct the three-block timing table.
3. `glm_dcm/run_fullrun_block_specific_glm.m` — estimate the full-run block-specific GLMs.
4. `glm_dcm/build_fullrun_block_specific_dcm.m` — extract VOIs and estimate four-node DCMs.
5. `dcm_export/export_dcm_a_matrices.m` — export participant-level A-matrix estimates.
6. `dcm_export/summarize_dcm_a_matrices.m` — generate group A-matrix summaries.
7. `qc/run_dcm_parameter_qc.m` — summarize posterior parameter evidence.
8. `qc/run_dcm_model_fit_qc.m` — summarize model fit and free energy.
9. `peb/run_peb_group_analysis.m` — run PEB/BMR/BMA group inference.
10. `fc/run_raw_fc_analysis.m` — calculate raw functional connectivity.
11. `residual_fc/run_residual_fc_analysis.m` — calculate task-regressed residual functional connectivity.

The cropped sensitivity pipeline is documented under `sensitivity/cropped_dcm/`.

## Derived results and figures

The repository includes the numerical exports used for the reported neuroimaging analyses. Main and supplementary figures are stored in `figures/`, and their principal tabular inputs are stored in `figures/source_data/`. Detailed numerical summaries are available in `results/tables_csv/`.

Raw FC values are near ceiling and are therefore interpreted cautiously as strongly task-locked shared BOLD synchrony. Task-regressed residual FC is supplied separately. DCM, raw FC, and residual FC estimate different quantities and should not be treated as interchangeable measures.

The physically cropped analysis in `sensitivity/cropped_dcm/` is provided as a sensitivity analysis of temporal model specification rather than as the primary DCM pipeline.

## Source data

Raw BOLD and structural imaging files are not redistributed. The source imaging dataset should be obtained from the study reported by:

Duncan, K. J., Pattamadilok, C., Knierim, I., & Devlin, J. T. (2009). *Consistency and variability in functional localisers*. NeuroImage, 46(4), 1018–1026. https://doi.org/10.1016/j.neuroimage.2009.03.014

This repository focuses on the neuroimaging analyses and derived outputs used in the accompanying manuscript.

## Reproducibility

- File-level SHA-256 hashes are listed in `FILE_MANIFEST.csv`.
- ROI reconstruction and voxelwise QC are provided under `resources/rois/` and `code/matlab/roi/`.
- Primary participant-level DCM values are in `data/derived/dcm/DCM_A_allsubjects_wide.csv`.
- DCM model-fit and posterior QC outputs are in `data/derived/dcm_model_fit/` and `data/derived/dcm_qc/`.
- Figure source data are in `figures/source_data/`.
- Cropped-segment sensitivity outputs and QA tables are in `sensitivity/cropped_dcm/`.

See `docs/REPRODUCIBILITY.md` for additional details.

## Citation

Citation metadata are provided in `CITATION.cff`.

Repository URL: https://github.com/Meysam-Kharghani/blockwise-dcm-fc

## License

Code and repository materials are released under the MIT License. See `LICENSE`.
