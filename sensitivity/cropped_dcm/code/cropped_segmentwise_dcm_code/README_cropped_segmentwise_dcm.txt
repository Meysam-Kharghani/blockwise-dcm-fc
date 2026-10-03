Cropped segment-wise DCM sensitivity analysis
==============================================

Purpose
-------
This auxiliary analysis estimates Block-01, Block-02, and Block-03 DCMs from
physically cropped time series. It is separate from the primary full-run,
block-specific analysis.

Scripts
-------
1. run_cropped_segmentwise_glm_3blocks.m
   Builds block-specific GLMs containing only scans from the selected segment.

2. build_dcm_spm12_cropped_segments.m
   Estimates four-node DCMs from the cropped GLMs. The default echo time is
   0.05 s, matching the acquisition TE of 50 ms.

3. check_cropped_dcm_timepoints.m
   Verifies that block-specific DCMs contain fewer time points than ALL models.

Configuration
-------------
The functions accept name-value path arguments. Repository-level defaults are
derived from code/matlab/project_paths.m when that helper is on the MATLAB path.
No machine-specific absolute paths are required.

Interpretation
--------------
This directory is provided as a sensitivity analysis and does not replace the
primary full-run DCM pipeline.

Derived outputs
---------------
The `sensitivity/cropped_dcm/derived/dcm_exports/` directory contains the participant-level and group DCM export set associated with this cropped analysis.
