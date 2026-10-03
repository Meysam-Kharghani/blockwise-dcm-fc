function check_ROI_masks_against_reference(generated_dir, reference_dir)
% Compare regenerated ROI masks against the distributed reference masks voxel-by-voxel.
%
% Usage:
%   check_ROI_masks_against_reference(fullfile(pwd,'generated_roi_masks'), ...
%                                     fullfile(pwd,'resources','rois'))

spm('defaults','FMRI');
files = {'IFG_L_mask.nii','V1_L_mask.nii','lOTC_mask.nii','vOTC_mask.nii'};
for i = 1:numel(files)
    gf = fullfile(generated_dir, files{i});
    rf = fullfile(reference_dir, files{i});
    if ~exist(gf,'file')
        fprintf('%s: generated file missing.\n', files{i}); continue;
    end
    if ~exist(rf,'file')
        fprintf('%s: reference file missing.\n', files{i}); continue;
    end
    Vg = spm_vol(gf); Yg = spm_read_vols(Vg) > 0.5;
    Vr = spm_vol(rf); Yr = spm_read_vols(Vr) > 0.5;
    same_dim = isequal(Vg.dim, Vr.dim);
    same_voxels = same_dim && isequal(Yg, Yr);
    mat_diff = max(abs(Vg.mat(:)-Vr.mat(:)));
    if same_voxels
        fprintf('%s: IDENTICAL voxel set; generated voxels=%d; reference voxels=%d; max affine diff=%.6g\n', ...
            files{i}, nnz(Yg), nnz(Yr), mat_diff);
    else
        fprintf('%s: DIFFERENT; differing voxels=%d; generated voxels=%d; reference voxels=%d; max affine diff=%.6g\n', ...
            files{i}, nnz(xor(Yg,Yr)), nnz(Yg), nnz(Yr), mat_diff);
    end
end
end
