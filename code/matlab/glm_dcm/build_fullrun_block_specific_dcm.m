function build_fullrun_block_specific_dcm
% =========================================================================
% build_fullrun_block_specific_dcm
%
% Purpose:
%   Build and estimate one separate fully-connected DCM for:
%     1) each run-wise ALL GLM
%     2) each block-wise GLM (block-01, block-02, block-03)
%
% For every subject:
%   GLM_run_1_ALL
%   GLM_run_2_ALL
%   GLM_run_1_block-01
%   GLM_run_1_block-02
%   GLM_run_1_block-03
%   GLM_run_2_block-01
%   GLM_run_2_block-02
%   GLM_run_2_block-03
%
% Output:
%   DCM_*.mat saved inside each GLM folder
%
% Notes:
%   - Fully connected DCM: all directed between-region connections ON
%   - All off-diagonal endogenous connections are enabled; standard SPM self-inhibition parameters are estimated internally.
%   - For later histogram analysis use:
%         DCM.Ep.A
% =========================================================================

clear;
clc;

% --------------------------- PATH CONFIGURATION ---------------------------
cfg = project_paths();
glm_root = cfg.glm_root;
roi_root = cfg.roi_root;

% choose what to process
process_all_models   = true;   % GLM_run_*_ALL
process_block_models = true;   % GLM_run_*_block-01/02/03

% fMRI / DCM settings
TE_value      = 0.05;   % seconds (50 ms)
voi_adjust    = 0;      % 0 = no adjustment
voi_session   = 1;      % each GLM folder assumed single-session
overwrite_dcm = false;  % false = skip already existing DCM

% ROI definitions
roi_defs = { ...
    'vOTC', fullfile(roi_root, 'vOTC_mask.nii'); ...
    'lOTC', fullfile(roi_root, 'lOTC_mask.nii'); ...
    'V1',   fullfile(roi_root, 'V1_L_mask.nii'); ...
    'IFG',  fullfile(roi_root, 'IFG_L_mask.nii') ...
};

block_names = {'block-01','block-02','block-03'};

% --------------------------- STARTUP CHECKS ------------------------------
if ~exist(glm_root, 'dir')
    error('GLM root folder not found:\n%s', glm_root);
end

if ~exist(roi_root, 'dir')
    error('ROI root folder not found:\n%s', roi_root);
end

for k = 1:size(roi_defs,1)
    if ~exist(roi_defs{k,2}, 'file')
        error('ROI mask not found: %s', roi_defs{k,2});
    end
end

% Initialize SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

% --------------------------- FIND SUBJECTS -------------------------------
sub_dirs = dir(fullfile(glm_root, 'sub-*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);
sub_dirs = sub_dirs(~ismember({sub_dirs.name}, {'.','..'}));

if isempty(sub_dirs)
    error('No subject folders found in:\n%s', glm_root);
end

fprintf('\n============================================================\n');
fprintf('Starting DCM estimation for ALL + BLOCK models\n');
fprintf('GLM root : %s\n', glm_root);
fprintf('ROI root : %s\n', roi_root);
fprintf('Subjects : %d\n', numel(sub_dirs));
fprintf('============================================================\n');

% logs
n_total_glm   = 0;
n_success_dcm = 0;
n_skipped     = 0;
n_failed      = 0;

% --------------------------- SUBJECT LOOP --------------------------------
for s = 1:numel(sub_dirs)

    sub_id  = sub_dirs(s).name;
    sub_dir = fullfile(sub_dirs(s).folder, sub_id);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Subject: %s\n', sub_id);
    fprintf('Subject folder: %s\n', sub_dir);
    fprintf('------------------------------------------------------------\n');

    % find all folders inside subject and then filter by name
    all_items = dir(sub_dir);
    all_items = all_items([all_items.isdir]);
    all_items = all_items(~ismember({all_items.name}, {'.','..'}));

    fprintf('Folders found inside subject directory:\n');
    for ii = 1:numel(all_items)
        fprintf('   - %s\n', all_items(ii).name);
    end

    glm_folders = [];
    for ii = 1:numel(all_items)
        this_name = all_items(ii).name;

        if is_target_glm_name(this_name, process_all_models, process_block_models, block_names)
            glm_folders = [glm_folders; all_items(ii)]; %#ok<AGROW>
        else
            fprintf('   [SKIP-NONGLM] %s\n', this_name);
        end
    end

    if isempty(glm_folders)
        fprintf('  [WARN] No target GLM folders found for %s that match requested patterns.\n', sub_id);
        continue;
    end

    % sort for stable processing order
    [~, idx_sort] = sort(lower({glm_folders.name}));
    glm_folders = glm_folders(idx_sort);

    fprintf('\nMatched GLM folders:\n');
    for ii = 1:numel(glm_folders)
        fprintf('   + %s\n', glm_folders(ii).name);
    end

    % ----------------------- GLM LOOP ------------------------------------
    for g = 1:numel(glm_folders)

        n_total_glm = n_total_glm + 1;

        glm_name = glm_folders(g).name;
        glm_dir  = fullfile(glm_folders(g).folder, glm_name);
        spm_file = fullfile(glm_dir, 'SPM.mat');

        fprintf('\n  >>> Processing GLM: %s\n', glm_name);
        fprintf('      Path: %s\n', glm_dir);

        % Check SPM.mat existence
        if ~exist(spm_file, 'file')
            fprintf('      [SKIP] Missing SPM.mat\n');
            n_skipped = n_skipped + 1;
            continue;
        end

        % Define DCM output name and path
        dcm_name = sprintf('DCM_%s_%s_%dROI.mat', sub_id, glm_name, size(roi_defs,1));
        dcm_path = fullfile(glm_dir, dcm_name);

        % Skip if DCM already exists and overwrite is false
        if exist(dcm_path, 'file') && ~overwrite_dcm
            fprintf('      [SKIP] DCM already exists: %s\n', dcm_name);
            n_skipped = n_skipped + 1;
            continue;
        end

        try
            % ------------------- Load SPM.mat ----------------------------
            tmp = load(spm_file, 'SPM');

            if ~isfield(tmp, 'SPM') || isempty(tmp.SPM)
                fprintf('      [FAIL] SPM structure missing in SPM.mat\n');
                n_failed = n_failed + 1;
                continue;
            end
            SPM = tmp.SPM;

            % working directory
            if ~isfield(SPM, 'swd') || isempty(SPM.swd)
                SPM.swd = glm_dir;
            end

            % session check
            if ~isfield(SPM, 'Sess') || isempty(SPM.Sess)
                fprintf('      [FAIL] SPM.Sess missing or empty\n');
                n_failed = n_failed + 1;
                continue;
            end

            % TR check
            if ~isfield(SPM, 'xY') || ~isfield(SPM.xY, 'RT') || isempty(SPM.xY.RT)
                fprintf('      [FAIL] TR not found in SPM.xY.RT\n');
                n_failed = n_failed + 1;
                continue;
            end
            TR = SPM.xY.RT;

            % experimental inputs
            if ~isfield(SPM.Sess(1), 'U') || isempty(SPM.Sess(1).U)
                fprintf('      [FAIL] No experimental input found in SPM.Sess(1).U\n');
                n_failed = n_failed + 1;
                continue;
            end
            nU = numel(SPM.Sess(1).U);

            fprintf('      TR = %.4f sec | Inputs = %d\n', TR, nU);

            % ------------------- Build or load VOIs ----------------------
            voi_files = cell(1, size(roi_defs,1));
            for r = 1:size(roi_defs,1)
                roi_label = roi_defs{r,1};
                roi_mask  = roi_defs{r,2};

                fprintf('      ROI %d/%d: %s\n', r, size(roi_defs,1), roi_label);

                voi_files{r} = ensure_voi_file( ...
                    glm_dir, spm_file, roi_label, roi_mask, voi_adjust, voi_session);

                fprintf('         VOI file: %s\n', voi_files{r});
            end

            % ------------------- Load VOI xY structs ----------------------
            xY_cell = cell(1, numel(voi_files));
            for r = 1:numel(voi_files)
                voi_tmp = load(voi_files{r}, 'xY');
                if ~isfield(voi_tmp, 'xY') || isempty(voi_tmp.xY)
                    error('VOI file missing valid xY: %s', voi_files{r});
                end
                voi_tmp.xY.name = roi_defs{r,1};
                xY_cell{r} = voi_tmp.xY;
            end
            xY = unify_xY_structs_local(xY_cell);

            % ------------------- DCM specification -----------------------
            nROI = numel(xY);

            % fully connected, no self-connections
            a = ones(nROI, nROI) - eye(nROI);

            % no modulatory / nonlinear effects
            b = zeros(nROI, nROI, nU);
            d = zeros(nROI, nROI, 0);
            c = ones(nROI, nU);

            include_inputs = ones(nU, 1);

            s_dcm = struct();
            s_dcm.name      = erase(dcm_name, '.mat');
            s_dcm.u         = include_inputs;
            s_dcm.delays    = repmat(TR, 1, nROI);
            s_dcm.TE        = TE_value;
            s_dcm.nonlinear = 0;
            s_dcm.two_state = 0;
            s_dcm.stochastic = 0;
            s_dcm.centre    = 1;
            s_dcm.induced   = 0;
            s_dcm.a         = a;
            s_dcm.b         = b;
            s_dcm.c         = c;
            s_dcm.d         = d;

            fprintf('      Specifying DCM...\n');
            DCM = spm_dcm_specify(SPM, xY, s_dcm);

            fprintf('      Estimating DCM...\n');
            DCM = spm_dcm_estimate(DCM);

            fprintf('      Saving DCM...\n');
            save(dcm_path, 'DCM', '-v7.3');

            fprintf('      [OK] Saved: %s\n', dcm_name);
            n_success_dcm = n_success_dcm + 1;

        catch ME
            fprintf('      [FAIL] %s\n', ME.message);
            n_failed = n_failed + 1;
        end
    end
end

% --------------------------- SUMMARY -------------------------------------
fprintf('\n============================================================\n');
fprintf('DCM BUILD SUMMARY\n');
fprintf('============================================================\n');
fprintf('Total matched GLM folders : %d\n', n_total_glm);
fprintf('Successfully estimated    : %d\n', n_success_dcm);
fprintf('Skipped                   : %d\n', n_skipped);
fprintf('Failed                    : %d\n', n_failed);
fprintf('============================================================\n');

end

% =========================================================================
% Helper: decide whether folder name is a target GLM
% =========================================================================
function tf = is_target_glm_name(folder_name, process_all_models, process_block_models, block_names)

tf = false;

% Match GLM_run_X_ALL
if process_all_models
    tok = regexp(folder_name, '^GLM_run_(\d+)_ALL$', 'tokens', 'once');
    if ~isempty(tok)
        tf = true;
        return;
    end
end

% Match GLM_run_X_block-YY
if process_block_models
    tok = regexp(folder_name, '^GLM_run_(\d+)_(block-\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        block_part = tok{2};
        if any(strcmp(block_part, block_names))
            tf = true;
            return;
        end
    end
end

end

% =========================================================================
% Helper: create VOI if missing, otherwise reuse
% =========================================================================
function voi_file = ensure_voi_file(glm_dir, spm_file, roi_label, roi_mask, voi_adjust, voi_session)

voi_file = fullfile(glm_dir, sprintf('VOI_%s_1.mat', roi_label));

if exist(voi_file, 'file')
    fprintf('         Reusing existing VOI\n');
    return;
end

if ~exist(roi_mask, 'file')
    error('ROI mask not found: %s', roi_mask);
end

matlabbatch = {};
matlabbatch{1}.spm.util.voi.spmmat = {spm_file};
matlabbatch{1}.spm.util.voi.adjust = voi_adjust;
matlabbatch{1}.spm.util.voi.session = voi_session;
matlabbatch{1}.spm.util.voi.name = roi_label;

% ROI from mask only
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {roi_mask};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;

% expression
matlabbatch{1}.spm.util.voi.expression = 'i1';

spm_jobman('run', matlabbatch);

if ~exist(voi_file, 'file')
    error('VOI file was not created: %s', voi_file);
end

end

% =========================================================================
% Helper: unify xY structs so all have same fields
% =========================================================================
function xY = unify_xY_structs_local(xY_cell)

all_fields = {};
for i = 1:numel(xY_cell)
    all_fields = union(all_fields, fieldnames(xY_cell{i}));
end

for i = 1:numel(xY_cell)
    missing = setdiff(all_fields, fieldnames(xY_cell{i}));
    for j = 1:numel(missing)
        xY_cell{i}.(missing{j}) = [];
    end
    xY_cell{i} = orderfields(xY_cell{i}, all_fields);
end

xY = [xY_cell{:}];

end
