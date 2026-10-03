# ROI masks used in the reported analyses

This directory contains the four binary NIfTI masks used for VOI extraction in the reported analyses. The masks are distributed so that the ROI definitions can be inspected and reproduced directly.

| ROI | Construction | Final MNI definition | Nonzero voxels | Approx. volume |
|---|---|---|---:|---:|
| IFG_L | 6-mm sphere | center `[-48, 28, 0]` | 123 | 984 mm³ |
| V1_L | 6-mm sphere | center `[-16, -98, -10]` | 123 | 984 mm³ |
| lOTC | box ROI | X `-58:-36`, Y `-86:-66`, Z `-18:6` | 1716 | 13,728 mm³ |
| vOTC | box ROI | X `-56:-32`, Y `-68:-44`, Z `-28:-2` | 2366 | 18,928 mm³ |

The sphere masks use a 2-mm MNI grid with dimensions 91×109×91 and affine origin `[90, -126, -72]`. The lOTC/vOTC masks use a 2-mm MNI grid with dimensions 91×109×91 and affine origin `[88, -124, -70]`.

Voxelwise regeneration QC found zero differing voxels for all four masks. The lOTC and vOTC masks overlap by 198 voxels (1,584 mm³); all other ROI pairs have zero overlap. Full QC is provided in `qc/`.

Reproducibility scripts are in `code/matlab/roi/`. The combined script `make_analysis_roi_masks.m` regenerates all four masks.
