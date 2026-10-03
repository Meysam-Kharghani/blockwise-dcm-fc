# Methods and output notes

## Primary DCM analysis

- DCM sample: 45 participants.
- Runs: two per participant.
- Estimated analytical levels per run: ALL, Block-01, Block-02, Block-03.
- Total primary DCMs: 360.
- Network: V1_L, lOTC, vOTC, IFG_L.
- Echo time used for DCM specification: 0.05 s (50 ms).
- Endogenous architecture: all 12 directed off-diagonal connections enabled, with standard SPM self-inhibition parameters.
- Primary temporal contrast: Block-03 minus mean(Block-01, Block-02).
- Primary PEB: 270 block-specific DCMs (45 participants × 2 runs × 3 blocks); ALL models are excluded from block-wise PEB inference.
- Frequentist FDR correction is applied across the 12 directed between-region DCM edges for each contrast family; self-connections are summarized separately.

## Full-run block-specific specification

Block labels refer to separate block-specific GLM/input definitions estimated from the complete run time series. The primary block-specific DCMs are not physically cropped. In each block-specific GLM, the target analytical block is represented explicitly while the complete run remains in the estimation data; activity not explained by that design remains in residual variance.

The A-matrix is the primary inferential target. Condition-specific B-matrix modulation is not part of this analysis.

## ROI definitions

The final masks used in the analyses are distributed under `resources/rois/`. IFG_L and V1_L are 6-mm spheres centered at MNI `[-48, 28, 0]` and `[-16, -98, -10]`, respectively. lOTC and vOTC are Duncan-derived box ROIs with the exact bounds documented in `resources/rois/README.md`. Voxelwise regeneration QC and overlap information are included in `resources/rois/qc/`.

## Functional connectivity

- Raw FC: Run 1 has 44 available participants and Run 2 has 45 because sub-37/run_1 is absent from the supplied FC exports.
- Residual FC is calculated after task-effect regression.
- Frequentist FDR correction is applied separately across the six undirected edges for raw FC and residual FC.
- Raw FC, residual FC, and DCM quantify different properties and are not treated as interchangeable estimands.

## Sensitivity analysis

The `sensitivity/cropped_dcm/` directory contains the physically segmented DCM analysis, participant-level DCM exports, and QA summaries. It is distinct from the primary full-run pipeline.
