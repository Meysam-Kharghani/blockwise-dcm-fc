function make_analysis_roi_masks(out_dir, reference_dir)
%MAKE_ANALYSIS_ROI_MASKS Generate the four ROI masks used in the reported analyses.
%
% This function reconstructs the voxel sets corresponding to the ROI masks
% used in the analyses:
%   1) IFG_L_mask.nii  : 6-mm sphere, center [-48  28   0]
%   2) V1_L_mask.nii   : 6-mm sphere, center [-16 -98 -10]
%   3) lOTC_mask.nii   : Duncan-derived box ROI
%   4) vOTC_mask.nii   : Duncan-derived box ROI
%
% Usage:
%   make_analysis_roi_masks(fullfile(pwd, 'generated_roi_masks'))
%   make_analysis_roi_masks(fullfile(pwd, 'generated_roi_masks'), fullfile(pwd, 'resources', 'rois'))
%
% If reference_dir is provided and contains the original masks, the script will
% compare the newly generated files voxel-by-voxel against the originals.
%
% Requirements:
%   SPM12 must be on the MATLAB path.

if nargin < 1 || isempty(out_dir)
    out_dir = fullfile(pwd, 'generated_roi_masks');
end
if nargin < 2
    reference_dir = '';
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

spm('defaults','FMRI');

% -------------------------------------------------------------------------
% Geometry 1: used by IFG_L and V1_L masks
% This matches the distributed IFG_L_mask.nii and V1_L_mask.nii geometry:
% dim = [91 109 91], voxel size = 2 mm, affine rows:
% [-2 0 0 90; 0 2 0 -126; 0 0 2 -72; 0 0 0 1]
% -------------------------------------------------------------------------
dim_sphere = [91 109 91];
mat_sphere = [-2 0 0 90; 0 2 0 -126; 0 0 2 -72; 0 0 0 1];

make_sphere_mask_exact([-48  28   0], 6, dim_sphere, mat_sphere, ...
    fullfile(out_dir, 'IFG_L_mask.nii'), 'sphere ROI');

make_sphere_mask_exact([-16 -98 -10], 6, dim_sphere, mat_sphere, ...
    fullfile(out_dir, 'V1_L_mask.nii'), 'Sphere ROI V1/2_L (6.0mm) at [-16 -98 -10]');

% -------------------------------------------------------------------------
% Geometry 2: used by Duncan-derived lOTC/vOTC box masks
% This matches the distributed lOTC_mask.nii and vOTC_mask.nii geometry:
% dim = [91 109 91], voxel size = 2 mm, affine rows:
% [-2 0 0 88; 0 2 0 -124; 0 0 2 -70; 0 0 0 1]
% -------------------------------------------------------------------------
dim_box = [91 109 91];
mat_box = [-2 0 0 88; 0 2 0 -124; 0 0 2 -70; 0 0 0 1];

% lOTC box: X = -58:-36, Y = -86:-66, Z = -18:6 in 2-mm MNI grid.
make_box_mask_exact([-58 -36], [-86 -66], [-18  6], dim_box, mat_box, ...
    fullfile(out_dir, 'lOTC_mask.nii'), 'lOTC Box ROI (Duncan 2009)');

% vOTC box: X = -56:-32, Y = -68:-44, Z = -28:-2 in 2-mm MNI grid.
make_box_mask_exact([-56 -32], [-68 -44], [-28 -2], dim_box, mat_box, ...
    fullfile(out_dir, 'vOTC_mask.nii'), 'vOTC Box ROI (Duncan 2009)');

% Report QC for all generated masks.
expected.IFG_L_mask = 123;
expected.V1_L_mask  = 123;
expected.lOTC_mask  = 1716;
expected.vOTC_mask  = 2366;

files = {'IFG_L_mask.nii','V1_L_mask.nii','lOTC_mask.nii','vOTC_mask.nii'};
fprintf('\n===== ROI QC summary =====\n');
for i = 1:numel(files)
    f = fullfile(out_dir, files{i});
    V = spm_vol(f);
    Y = spm_read_vols(V) > 0.5;
    nvox = nnz(Y);
    key = erase(files{i}, '.nii');
    fprintf('%s: %d voxels', files{i}, nvox);
    if isfield(expected, key)
        if nvox == expected.(key)
            fprintf('  [OK]\n');
        else
            fprintf('  [WARNING: expected %d]\n', expected.(key));
        end
    else
        fprintf('\n');
    end
end

% Optional voxelwise comparison against reference masks.
if ~isempty(reference_dir)
    fprintf('\n===== Voxelwise comparison with reference_dir =====\n');
    for i = 1:numel(files)
        ref_file = fullfile(reference_dir, files{i});
        gen_file = fullfile(out_dir, files{i});
        if exist(ref_file, 'file')
            compare_masks_exact(gen_file, ref_file);
        else
            fprintf('%s: reference file not found, skipped.\n', files{i});
        end
    end
end

fprintf('\nDone. ROI masks written to:\n%s\n', out_dir);
end

% =========================================================================
function make_sphere_mask_exact(center_mm, radius_mm, dim, mat, out_nii, descrip)
[X,Y,Z] = make_mni_grid_zero_based(dim, mat);
D = sqrt((X-center_mm(1)).^2 + (Y-center_mm(2)).^2 + (Z-center_mm(3)).^2);
mask = uint8(D <= radius_mm);
write_uint8_nifti(mask, dim, mat, out_nii, descrip);
end

% =========================================================================
function make_box_mask_exact(xlim_mm, ylim_mm, zlim_mm, dim, mat, out_nii, descrip)
[X,Y,Z] = make_mni_grid_zero_based(dim, mat);
mask = uint8(X >= xlim_mm(1) & X <= xlim_mm(2) & ...
             Y >= ylim_mm(1) & Y <= ylim_mm(2) & ...
             Z >= zlim_mm(1) & Z <= zlim_mm(2));
write_uint8_nifti(mask, dim, mat, out_nii, descrip);
end

% =========================================================================
function [X,Y,Z] = make_mni_grid_zero_based(dim, mat)
% Use the same voxel-coordinate convention used in the distributed masks:
% voxel index starts at 0 for conversion to the NIfTI sform coordinates.
[ix,iy,iz] = ndgrid(0:dim(1)-1, 0:dim(2)-1, 0:dim(3)-1);
X = mat(1,1).*ix + mat(1,2).*iy + mat(1,3).*iz + mat(1,4);
Y = mat(2,1).*ix + mat(2,2).*iy + mat(2,3).*iz + mat(2,4);
Z = mat(3,1).*ix + mat(3,2).*iy + mat(3,3).*iz + mat(3,4);
end

% =========================================================================
function write_uint8_nifti(mask, dim, mat, out_nii, descrip)
V = struct();
V.fname = out_nii;
V.dim = dim;
V.dt = [spm_type('uint8') 0];
V.mat = mat;
V.pinfo = [1; 0; 0];
V.descrip = descrip;
spm_write_vol(V, mask);
end

% =========================================================================
function compare_masks_exact(gen_file, ref_file)
Vg = spm_vol(gen_file); Yr = [];
Yg = spm_read_vols(Vg) > 0.5;
Vr = spm_vol(ref_file);
Yr = spm_read_vols(Vr) > 0.5;

same_dim = isequal(Vg.dim, Vr.dim);
same_voxels = same_dim && isequal(Yg, Yr);
mat_diff = max(abs(Vg.mat(:) - Vr.mat(:)));

[~,name,ext] = fileparts(gen_file);
if same_voxels
    fprintf('%s%s: voxel set IDENTICAL; max affine difference = %.6g\n', name, ext, mat_diff);
else
    fprintf('%s%s: NOT identical; differing voxels = %d; max affine difference = %.6g\n', ...
        name, ext, nnz(xor(Yg,Yr)), mat_diff);
end
end
