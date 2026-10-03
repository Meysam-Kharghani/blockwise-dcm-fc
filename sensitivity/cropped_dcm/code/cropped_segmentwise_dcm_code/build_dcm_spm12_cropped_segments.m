function build_dcm_spm12_cropped_segments(varargin)
% =========================================================================
% build_dcm_spm12_cropped_segments
% =========================================================================
% Estimate DCMs from GLMs produced by run_cropped_segmentwise_glm_3blocks.
%
% Key assumption:
%   Block GLMs now contain only the scans from the relevant block segment.
%   Therefore each DCM.Y.y should have approximately one-third of the run.
%
% Default paths are obtained from project_paths.m.
%
% Output:
%   DCM_*.mat saved inside each GLM folder.
%
% Optional name-value inputs:
%   'glm_root', 'roi_root', 'TE_value', 'overwrite_dcm',
%   'process_all_models', 'process_block_models'
% =========================================================================

cfg = project_paths();
P = inputParser;
P.addParameter('glm_root', cfg.cropped_glm_root, @ischar);
P.addParameter('roi_root', cfg.roi_root, @ischar);
P.addParameter('TE_value', 0.05, @isnumeric);       % acquisition TE = 50 ms
P.addParameter('overwrite_dcm', true, @islogical);
P.addParameter('process_all_models', true, @islogical);
P.addParameter('process_block_models', true, @islogical);
P.addParameter('voi_adjust', 0, @isnumeric);
P.addParameter('voi_session', 1, @isnumeric);
P.parse(varargin{:});
S = P.Results;

roi_defs = { ...
    'vOTC',  fullfile(S.roi_root, 'vOTC_mask.nii'); ...
    'lOTC',  fullfile(S.roi_root, 'lOTC_mask.nii'); ...
    'V1_L',  fullfile(S.roi_root, 'V1_L_mask.nii'); ...
    'IFG_L', fullfile(S.roi_root, 'IFG_L_mask.nii') ...
};
block_names = {'block-01','block-02','block-03'};

if ~exist(S.glm_root, 'dir'), error('GLM root not found: %s', S.glm_root); end
if ~exist(S.roi_root, 'dir'), error('ROI root not found: %s', S.roi_root); end
for k = 1:size(roi_defs,1)
    if ~exist(roi_defs{k,2}, 'file')
        error('ROI mask not found: %s', roi_defs{k,2});
    end
end

spm('defaults','FMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

fprintf('\n============================================================\n');
fprintf('CROPPED SEGMENT-WISE DCM ESTIMATION\n');
fprintf('GLM root : %s\n', S.glm_root);
fprintf('ROI root : %s\n', S.roi_root);
fprintf('TE       : %.4f sec\n', S.TE_value);
fprintf('============================================================\n');

sub_dirs = dir(fullfile(S.glm_root, 'sub-*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);

n_total = 0; n_ok = 0; n_skip = 0; n_fail = 0;
log_rows = {};

for si = 1:numel(sub_dirs)
    sub_id = sub_dirs(si).name;
    sub_dir = fullfile(sub_dirs(si).folder, sub_id);
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Subject: %s\n', sub_id);

    all_items = dir(sub_dir);
    all_items = all_items([all_items.isdir]);
    all_items = all_items(~ismember({all_items.name}, {'.','..','GLM_blocks_summaries'}));

    glm_folders = [];
    for ii = 1:numel(all_items)
        nm = all_items(ii).name;
        if is_target_glm_name(nm, S.process_all_models, S.process_block_models, block_names)
            glm_folders = [glm_folders; all_items(ii)]; %#ok<AGROW>
        end
    end
    [~, ix] = sort(lower({glm_folders.name}));
    glm_folders = glm_folders(ix);

    for gi = 1:numel(glm_folders)
        n_total = n_total + 1;
        glm_name = glm_folders(gi).name;
        glm_dir  = fullfile(glm_folders(gi).folder, glm_name);
        spm_file = fullfile(glm_dir, 'SPM.mat');
        fprintf('\n  >>> %s\n', glm_name);

        if ~exist(spm_file, 'file')
            fprintf('      [SKIP] missing SPM.mat\n');
            n_skip = n_skip + 1; continue;
        end

        dcm_name = sprintf('DCM_%s_%s_%dROI.mat', sub_id, glm_name, size(roi_defs,1));
        dcm_path = fullfile(glm_dir, dcm_name);
        if exist(dcm_path, 'file') && ~S.overwrite_dcm
            fprintf('      [SKIP] DCM exists\n');
            n_skip = n_skip + 1; continue;
        end

        try
            tmp = load(spm_file, 'SPM');
            SPM = tmp.SPM;
            if ~isfield(SPM, 'swd') || isempty(SPM.swd), SPM.swd = glm_dir; end
            if ~isfield(SPM, 'xY') || ~isfield(SPM.xY, 'RT') || isempty(SPM.xY.RT)
                error('SPM.xY.RT not found');
            end
            TR = SPM.xY.RT;
            if ischar(SPM.xY.P)
                nScans = size(SPM.xY.P, 1);
            else
                nScans = numel(SPM.xY.P);
            end
            if ~isfield(SPM.Sess(1), 'U') || isempty(SPM.Sess(1).U)
                error('No input found in SPM.Sess(1).U');
            end
            nU = numel(SPM.Sess(1).U);
            fprintf('      TR=%.3f | nScans=%d | nU=%d\n', TR, nScans, nU);

            % Create/reuse VOIs from masks
            voi_files = cell(1, size(roi_defs,1));
            for r = 1:size(roi_defs,1)
                roi_label = roi_defs{r,1};
                roi_mask = roi_defs{r,2};
                voi_files{r} = ensure_voi_file(glm_dir, spm_file, roi_label, roi_mask, S.voi_adjust, S.voi_session);
            end

            xY_cell = cell(1, numel(voi_files));
            for r = 1:numel(voi_files)
                V = load(voi_files{r}, 'xY');
                V.xY.name = roi_defs{r,1};
                xY_cell{r} = V.xY;
            end
            xY = unify_xY_structs_local(xY_cell);
            nROI = numel(xY);

            % Fully-connected intrinsic A matrix; no B modulation in this model.
            a = ones(nROI, nROI) - eye(nROI);
            b = zeros(nROI, nROI, nU);
            c = ones(nROI, nU);
            d = zeros(nROI, nROI, 0);

            s_dcm = struct();
            s_dcm.name       = erase(dcm_name, '.mat');
            s_dcm.u          = ones(nU, 1);
            s_dcm.delays     = repmat(TR, 1, nROI);
            s_dcm.TE         = S.TE_value;
            s_dcm.nonlinear  = 0;
            s_dcm.two_state  = 0;
            s_dcm.stochastic = 0;
            s_dcm.centre     = 1;
            s_dcm.induced    = 0;
            s_dcm.a          = a;
            s_dcm.b          = b;
            s_dcm.c          = c;
            s_dcm.d          = d;

            DCM = spm_dcm_specify(SPM, xY, s_dcm);
            DCM = spm_dcm_estimate(DCM);
            save(dcm_path, 'DCM', '-v7.3');

            fprintf('      [OK] saved: %s\n', dcm_name);
            n_ok = n_ok + 1;
            log_rows(end+1,:) = {sub_id, glm_name, nScans, TR, nU, S.TE_value, dcm_path, 'OK', ''}; %#ok<AGROW>
        catch ME
            fprintf('      [FAIL] %s\n', ME.message);
            n_fail = n_fail + 1;
            log_rows(end+1,:) = {sub_id, glm_name, NaN, NaN, NaN, S.TE_value, dcm_path, 'FAIL', ME.message}; %#ok<AGROW>
        end
    end
end

if ~isempty(log_rows)
    L = cell2table(log_rows, 'VariableNames', ...
        {'subject','glm_name','n_scans','TR','n_inputs','TE','dcm_path','status','message'});
    writetable(L, fullfile(S.glm_root, 'cropped_segmentwise_dcm_estimation_log.csv'));
end

fprintf('\n============================================================\n');
fprintf('CROPPED DCM SUMMARY\n');
fprintf('Matched GLMs : %d\n', n_total);
fprintf('Estimated    : %d\n', n_ok);
fprintf('Skipped      : %d\n', n_skip);
fprintf('Failed       : %d\n', n_fail);
fprintf('Log file     : %s\n', fullfile(S.glm_root, 'cropped_segmentwise_dcm_estimation_log.csv'));
fprintf('============================================================\n');
end

%% =========================================================================
function tf = is_target_glm_name(folder_name, process_all_models, process_block_models, block_names)
    tf = false;
    if process_all_models && ~isempty(regexp(folder_name, '^GLM_run_(\d+)_ALL$', 'once'))
        tf = true; return;
    end
    tok = regexp(folder_name, '^GLM_run_(\d+)_(block-\d+)$', 'tokens', 'once');
    if process_block_models && ~isempty(tok) && any(strcmp(tok{2}, block_names))
        tf = true; return;
    end
end

function voi_file = ensure_voi_file(glm_dir, spm_file, roi_label, roi_mask, voi_adjust, voi_session)
    voi_file = fullfile(glm_dir, sprintf('VOI_%s_1.mat', roi_label));
    if exist(voi_file, 'file')
        fprintf('      VOI reuse: %s\n', roi_label);
        return;
    end
    matlabbatch = {};
    matlabbatch{1}.spm.util.voi.spmmat = {spm_file};
    matlabbatch{1}.spm.util.voi.adjust = voi_adjust;
    matlabbatch{1}.spm.util.voi.session = voi_session;
    matlabbatch{1}.spm.util.voi.name = roi_label;
    matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {roi_mask};
    matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
    matlabbatch{1}.spm.util.voi.expression = 'i1';
    spm_jobman('run', matlabbatch);
    if ~exist(voi_file, 'file')
        error('VOI file was not created: %s', voi_file);
    end
end

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
