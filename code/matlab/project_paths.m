function cfg = project_paths()
%PROJECT_PATHS Repository and external-data path configuration.
%
% Environment variables:
%   DCM_FC_DATA_ROOT  Root of the external neuroimaging dataset (required for
%                     first-level reprocessing).
%   DCM_FC_WORK_ROOT  Writable processing directory. Defaults to <repo>/work.
%   DCM_FC_ROI_ROOT   Directory containing ROI NIfTI masks. Defaults to
%                     <repo>/resources/rois.
%   SPM12_DIR         Optional SPM12 installation directory.
%   CONN_DIR          Optional CONN installation directory.
%
% No machine-specific absolute paths are stored in the repository.

this_file = mfilename('fullpath');
matlab_dir = fileparts(this_file);
code_dir = fileparts(matlab_dir);
cfg.repo_root = fileparts(code_dir);

cfg.data_root = getenv('DCM_FC_DATA_ROOT');
if isempty(cfg.data_root)
    cfg.data_root = fullfile(cfg.repo_root, 'external_data');
end

cfg.work_root = getenv('DCM_FC_WORK_ROOT');
if isempty(cfg.work_root)
    cfg.work_root = fullfile(cfg.repo_root, 'work');
end

cfg.roi_root = getenv('DCM_FC_ROI_ROOT');
if isempty(cfg.roi_root)
    cfg.roi_root = fullfile(cfg.repo_root, 'resources', 'rois');
end

cfg.spm12_dir = getenv('SPM12_DIR');
cfg.conn_dir = getenv('CONN_DIR');

cfg.block_timing_csv = fullfile(cfg.repo_root, 'data', 'derived', ...
    'block_timing', 'FC_blocks_COND.csv');
cfg.glm_root = fullfile(cfg.work_root, 'glm');
cfg.cropped_glm_root = fullfile(cfg.work_root, 'glm_cropped');
cfg.dcm_export_root = fullfile(cfg.work_root, 'dcm_exports');
cfg.raw_fc_root = fullfile(cfg.work_root, 'fc_raw');
cfg.residual_fc_root = fullfile(cfg.work_root, 'fc_residual');
cfg.peb_root = fullfile(cfg.work_root, 'peb');
cfg.figure_root = fullfile(cfg.work_root, 'figures');

cfg.task_json = fullfile(cfg.data_root, 'task-onebacktask_bold.json');
end
