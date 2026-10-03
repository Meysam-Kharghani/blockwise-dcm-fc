function run_residual_fc_analysis(varargin)
% RUN_RESIDUAL_FC_ANALYSIS
%
% Computes residual functional connectivity after regressing out nuisance
% effects AND task-evoked effects from ROI time series.
%
% Task/event regressors are added to the nuisance model before FC is computed:
%
%   Y_residual = Y_raw - X_full * beta
%
% where X_full includes intercept, motion, ART outliers, WM/CSF PCs,
% high-pass DCT bases, and task regressors convolved with the canonical HRF.
%
% FC is then computed on Y_residual for ALL, Block-01, Block-02, Block-03.
%
% -------------------------------------------------------------------------
% Example:
%   cfg = project_paths();
%   run_residual_fc_analysis('base_dir', cfg.data_root, ...
%       'output_root', cfg.residual_fc_root, ...
%       'cond_csv', cfg.block_timing_csv, ...
%       'json_file', cfg.task_json, ...
%       'roi_dir', cfg.roi_root);
%
% Optional settings:
%   'task_events_dir'       root folder to search for events.tsv files.
%                           default: base_dir
%   'task_regressor_mode'   'events_tsv' | 'block_csv' | 'none' | 'auto'
%                           default: 'auto'
%   'task_categories'       cellstr of category labels to model. If empty,
%                           all detected event categories are used.
%   'include_derivatives'   true/false, include temporal derivatives of task
%                           regressors. default: false
%
% Residual FC is computed after removal of task and nuisance effects.


%% =========================
% DEFAULT SETTINGS
% ==========================
p = inputParser;
p.FunctionName = mfilename;
cfg = project_paths();

addParameter(p, 'base_dir', cfg.data_root, @ischar);
addParameter(p, 'output_root', cfg.residual_fc_root, @ischar);
addParameter(p, 'cond_csv', cfg.block_timing_csv, @ischar);
addParameter(p, 'json_file', cfg.task_json, @ischar);
addParameter(p, 'roi_dir', cfg.roi_root, @ischar);
addParameter(p, 'task_events_dir', '', @ischar);

addParameter(p, 'roi_files', { ...
    'IFG_L_mask.nii', ...
    'lOTC_mask.nii', ...
    'V1_L_mask.nii', ...
    'vOTC_mask.nii'}, @iscell);

addParameter(p, 'roi_names', { ...
    'IFG_L', ...
    'lOTC', ...
    'V1_L', ...
    'vOTC'}, @iscell);

addParameter(p, 'make_figs', true, @islogical);
addParameter(p, 'prefer_smoothed_functional', true, @islogical);
addParameter(p, 'do_zscore_after_denoise', true, @islogical);
addParameter(p, 'use_dct_highpass', true, @islogical);
addParameter(p, 'highpass_cutoff_sec', 128, @isnumeric);
addParameter(p, 'wm_n_components', 5, @isnumeric);
addParameter(p, 'csf_n_components', 5, @isnumeric);
addParameter(p, 'min_valid_timepoints_fc', 10, @isnumeric);
addParameter(p, 'low_variance_threshold', 1e-8, @isnumeric);
addParameter(p, 'clip_r_for_fisher', 0.999999, @isnumeric);

% Task regression settings
addParameter(p, 'task_regressor_mode', 'auto', @ischar); % auto/events_tsv/block_csv/none
addParameter(p, 'task_categories', {}, @iscell);         % if empty, use all detected categories
addParameter(p, 'event_label_column', 'auto', @ischar);  % auto/trial_type/condition/etc
addParameter(p, 'include_derivatives', false, @islogical);
addParameter(p, 'hrf_oversampling', 16, @isnumeric);
addParameter(p, 'exclude_event_labels', {'fixation','rest','baseline','n/a','na','null','none'}, @iscell);

% Optional subset selection. Empty lists process all available subjects and runs.
addParameter(p, 'restrict_subjects', {}, @iscell);
addParameter(p, 'restrict_runs', {}, @iscell);

parse(p, varargin{:});
S = p.Results;

base_dir    = S.base_dir;
output_root = S.output_root;
cond_csv    = S.cond_csv;
json_file   = S.json_file;
roi_dir     = S.roi_dir;
roi_files   = S.roi_files;
roi_names   = S.roi_names;

if isempty(S.task_events_dir)
    task_events_dir = base_dir;
else
    task_events_dir = S.task_events_dir;
end

%% =========================
% INIT
% ==========================
clc;
close all;

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

log_file = fullfile(output_root, 'ResidualFC_processing_log.txt');
append_log(log_file, '=====================================================');
append_log(log_file, sprintf('[START] %s', datestr(now, 31)));
append_log(log_file, sprintf('[INFO] output_root = %s', output_root));
append_log(log_file, sprintf('[INFO] task_regressor_mode = %s', S.task_regressor_mode));

if isempty(which('spm'))
    error('SPM12 is not on the MATLAB path.');
end

spm('Defaults','fMRI');
spm_jobman('initcfg');

%% =========================
% LOAD TR FROM JSON
% ==========================
if ~exist(json_file, 'file')
    error('JSON file not found: %s', json_file);
end

json_txt = fileread(json_file);
json_dat = jsondecode(json_txt);

if isfield(json_dat, 'RepetitionTime')
    TR = double(json_dat.RepetitionTime);
else
    error('RepetitionTime not found in JSON.');
end

if isfield(json_dat, 'StartTime')
    StartTime = double(json_dat.StartTime);
else
    StartTime = 0;
end

msg = sprintf('[INFO] TR = %.6f sec | StartTime = %.6f sec', TR, StartTime);
fprintf('%s\n', msg);
append_log(log_file, msg);

%% =========================
% LOAD CONDITION CSV FOR BLOCK SELECTION
% ==========================
if ~exist(cond_csv, 'file')
    error('Condition CSV not found: %s', cond_csv);
end

Tcond = readtable_robust(cond_csv);

required_cols = {'subject','session','condition','onset','duration'};
for i = 1:numel(required_cols)
    if isempty(find_varname(Tcond, required_cols{i}))
        error('Missing required column in FC_blocks_COND.csv: %s', required_cols{i});
    end
end

%% =========================
% FIND SUBJECTS
% ==========================
d = dir(fullfile(base_dir, 'sub-*'));
d = d([d.isdir]);
subjects = sort({d.name}');

if isempty(subjects)
    error('No subject folders found under: %s', base_dir);
end

fprintf('[INFO] Found %d subjects.\n', numel(subjects));
append_log(log_file, sprintf('[INFO] Found %d subject folders.', numel(subjects)));

if ~isempty(S.restrict_subjects)
    subjects = intersect(subjects, S.restrict_subjects, 'stable');
    append_log(log_file, sprintf('[INFO] Restricted subjects count = %d', numel(subjects)));
end

%% =========================
% MAIN LOOP
% ==========================
for s = 1:numel(subjects)
    subj = subjects{s};
    fprintf('\n=====================================================\n');
    fprintf('[SUBJECT] %s\n', subj);
    fprintf('=====================================================\n');
    append_log(log_file, sprintf('[SUBJECT] %s', subj));

    subj_dir = fullfile(base_dir, subj, 'func');
    if ~exist(subj_dir, 'dir')
        warning('No func directory for %s. Skipping.', subj);
        append_log(log_file, sprintf('[WARN] No func directory for %s. Skipping.', subj));
        continue;
    end

    run_dirs = dir(fullfile(subj_dir, 'run_*'));
    run_dirs = run_dirs([run_dirs.isdir]);
    run_names = sort({run_dirs.name}');

    if isempty(run_names)
        warning('No run_* directories for %s. Skipping.', subj);
        append_log(log_file, sprintf('[WARN] No run_* directories for %s. Skipping.', subj));
        continue;
    end

    if ~isempty(S.restrict_runs)
        run_names = intersect(run_names, S.restrict_runs, 'stable');
    end

    for rr = 1:numel(run_names)
        run_name = run_names{rr};
        fprintf('\n-------------------------------\n');
        fprintf('[RUN] %s | %s\n', subj, run_name);
        fprintf('-------------------------------\n');
        append_log(log_file, sprintf('[RUN] %s | %s', subj, run_name));

        try
            run_dir = fullfile(subj_dir, run_name);
            run_outdir = fullfile(output_root, subj, run_name);
            if ~exist(run_outdir, 'dir')
                mkdir(run_outdir);
            end

            files = get_subject_run_files(run_dir, subj, run_name, S.prefer_smoothed_functional);

            fprintf('[FILES] func:   %s\n', files.func);
            fprintf('[FILES] motion: %s\n', files.motion);
            fprintf('[FILES] ART:    %s\n', files.art);
            fprintf('[FILES] WM:     %s\n', files.wm);
            fprintf('[FILES] CSF:    %s\n', files.csf);

            append_log(log_file, sprintf('[FILES] func=%s', files.func));
            append_log(log_file, sprintf('[FILES] motion=%s', files.motion));
            append_log(log_file, sprintf('[FILES] art=%s', files.art));
            append_log(log_file, sprintf('[FILES] wm=%s', files.wm));
            append_log(log_file, sprintf('[FILES] csf=%s', files.csf));

            if isempty(files.func) || ~exist(files.func, 'file')
                error('Functional file not found.');
            end

            Vfunc = spm_vol(files.func);
            nScans = numel(Vfunc);
            append_log(log_file, sprintf('[INFO] %s %s nScans=%d', subj, run_name, nScans));

            % ===== load motion =====
            motion = [];
            if ~isempty(files.motion) && exist(files.motion, 'file')
                motion = load(files.motion);
                motion = match_rows_to_scans(motion, nScans);
            end

            % ===== load ART =====
            art_regs = [];
            if ~isempty(files.art) && exist(files.art, 'file')
                art_regs = load_art_regressors(files.art, nScans, log_file);
            end

            % ===== extract ROI signals =====
            Yraw = nan(nScans, numel(roi_files));
            roi_used_paths = cell(numel(roi_files),1);

            for r = 1:numel(roi_files)
                roi_path = fullfile(roi_dir, roi_files{r});
                if ~exist(roi_path, 'file')
                    error('ROI file not found: %s', roi_path);
                end

                roi_use = roi_path;
                if ~same_space(roi_path, files.func)
                    roi_use = reslice_roi_to_func_space(roi_path, files.func, run_outdir, ['roi_' roi_names{r} '_resliced'], log_file);
                else
                    append_log(log_file, sprintf('[ROI] same-space used directly: %s', roi_path));
                end

                roi_used_paths{r} = roi_use;
                Yraw(:,r) = extract_mean_timeseries_from_mask(files.func, roi_use);
            end

            raw_tbl = array2table(Yraw, 'VariableNames', matlab.lang.makeValidName(roi_names));
            writetable(raw_tbl, fullfile(run_outdir, 'roi_timeseries_raw.csv'));

            % ===== WM / CSF nuisance =====
            wm_ts = [];
            csf_ts = [];

            if ~isempty(files.wm) && exist(files.wm, 'file')
                wm_use = files.wm;
                if ~same_space(files.wm, files.func)
                    wm_use = reslice_roi_to_func_space(files.wm, files.func, run_outdir, 'wm_resliced', log_file);
                end
                wm_ts = extract_voxelwise_pca(files.func, wm_use, S.wm_n_components);
            else
                append_log(log_file, '[WARN] WM file missing.');
            end

            if ~isempty(files.csf) && exist(files.csf, 'file')
                csf_use = files.csf;
                if ~same_space(files.csf, files.func)
                    csf_use = reslice_roi_to_func_space(files.csf, files.func, run_outdir, 'csf_resliced', log_file);
                end
                csf_ts = extract_voxelwise_pca(files.func, csf_use, S.csf_n_components);
            else
                append_log(log_file, '[WARN] CSF file missing.');
            end

            % ===== build base nuisance design =====
            X = ones(nScans,1);
            X_names = {'Intercept'};

            if ~isempty(motion)
                X = [X, motion]; %#ok<AGROW>
                for k = 1:size(motion,2)
                    X_names{end+1} = sprintf('motion_%02d', k); %#ok<AGROW>
                end
            end

            if ~isempty(art_regs)
                X = [X, art_regs]; %#ok<AGROW>
                for k = 1:size(art_regs,2)
                    X_names{end+1} = sprintf('art_%02d', k); %#ok<AGROW>
                end
            end

            if ~isempty(wm_ts)
                X = [X, wm_ts]; %#ok<AGROW>
                for k = 1:size(wm_ts,2)
                    X_names{end+1} = sprintf('wmPC_%02d', k); %#ok<AGROW>
                end
            end

            if ~isempty(csf_ts)
                X = [X, csf_ts]; %#ok<AGROW>
                for k = 1:size(csf_ts,2)
                    X_names{end+1} = sprintf('csfPC_%02d', k); %#ok<AGROW>
                end
            end

            if S.use_dct_highpass
                Xhp = make_dct_basis(nScans, TR, S.highpass_cutoff_sec);
                X = [X, Xhp]; %#ok<AGROW>
                for k = 1:size(Xhp,2)
                    X_names{end+1} = sprintf('DCT_%02d', k); %#ok<AGROW>
                end
            end

            % ===== build task design regressors =====
            [Xtask, task_names, task_source] = build_task_regressors_for_run( ...
                S.task_regressor_mode, task_events_dir, run_dir, subj, run_name, ...
                Tcond, nScans, TR, StartTime, S.event_label_column, ...
                S.task_categories, S.exclude_event_labels, S.hrf_oversampling, ...
                S.include_derivatives, log_file);

            if ~isempty(Xtask)
                X = [X, Xtask]; %#ok<AGROW>
                X_names = [X_names, task_names]; %#ok<AGROW>
            end

            % ===== clean design matrix =====
            [X, X_names, removed_cols] = clean_design_matrix(X, X_names, log_file);

            design_tbl = array2table(X, 'VariableNames', matlab.lang.makeValidName(X_names));
            writetable(design_tbl, fullfile(run_outdir, 'residual_fc_design_matrix.csv'));

            design_info = table((1:numel(X_names))', X_names(:), 'VariableNames', {'column','name'});
            writetable(design_info, fullfile(run_outdir, 'residual_fc_design_columns.csv'));

            if ~isempty(removed_cols)
                writecell(removed_cols(:), fullfile(run_outdir, 'residual_fc_removed_design_columns.csv'));
            end

            append_log(log_file, sprintf('[TASK] source=%s | n_task_cols=%d', task_source, size(Xtask,2)));

            % ===== residualize ROI time series =====
            Yresid = residualize_timeseries(Yraw, X, S.do_zscore_after_denoise);

            resid_tbl = array2table(Yresid, 'VariableNames', matlab.lang.makeValidName(roi_names));
            writetable(resid_tbl, fullfile(run_outdir, 'roi_timeseries_residual_taskregressed.csv'));

            save(fullfile(run_outdir, 'run_workspace_residual_fc.mat'), ...
                'subj','run_name','files','TR','StartTime','roi_names','roi_files','roi_used_paths', ...
                'Yraw','Yresid','motion','art_regs','wm_ts','csf_ts','X','X_names','Xtask','task_names','task_source');

            % ===== Residual FC ALL =====
            compute_and_save_fc(Yresid, roi_names, run_outdir, 'ALL', ...
                S.make_figs, S.min_valid_timepoints_fc, S.low_variance_threshold, S.clip_r_for_fisher, ...
                subj, run_name, log_file);

            % ===== Residual block-wise FC =====
            session_num = parse_run_number(run_name);
            Tsub = get_subject_condition_rows(Tcond, subj, session_num);

            if isempty(Tsub)
                warning('No condition rows found for %s session %d', subj, session_num);
                append_log(log_file, sprintf('[WARN] No condition rows found for %s session %d', subj, session_num));
            else
                cond_var = find_varname(Tsub, 'condition');
                cond_values = unique(Tsub.(cond_var));
                cond_values = cond_values(~isnan(cond_values));

                for cc = 1:numel(cond_values)
                    cond_num = cond_values(cc);
                    rows = Tsub(Tsub.(cond_var) == cond_num, :);

                    block_idx = make_block_scan_index(rows, nScans, TR, StartTime);
                    append_log(log_file, sprintf('[BLOCK] %s %s Block-%02d selected_scans=%d', subj, run_name, cond_num, sum(block_idx)));

                    if sum(block_idx) < S.min_valid_timepoints_fc
                        warning('Too few scans for %s %s Block-%02d. Skipping.', subj, run_name, cond_num);
                        append_log(log_file, sprintf('[WARN] Too few scans for %s %s Block-%02d', subj, run_name, cond_num));
                        continue;
                    end

                    Yblk = Yresid(block_idx, :);
                    compute_and_save_fc(Yblk, roi_names, run_outdir, sprintf('Block-%02d', cond_num), ...
                        S.make_figs, S.min_valid_timepoints_fc, S.low_variance_threshold, S.clip_r_for_fisher, ...
                        subj, run_name, log_file);
                end
            end

            fprintf('[DONE] %s %s\n', subj, run_name);
            append_log(log_file, sprintf('[DONE] %s %s', subj, run_name));

        catch ME
            warning('[FAILED] %s %s: %s', subj, run_name, ME.message);
            fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            append_log(log_file, sprintf('[FAILED] %s %s: %s', subj, run_name, ME.message));
            append_log(log_file, getReport(ME, 'extended', 'hyperlinks', 'off'));
        end
    end
end

%% =========================
% GROUP-LEVEL EXPORT
% ==========================
try
    append_log(log_file, '[GROUP] Starting group-level export...');
    export_group_level_fc(output_root, roi_names, log_file);
    append_log(log_file, '[GROUP] Group-level export finished.');
catch MEg
    warning('Group-level export failed: %s', MEg.message);
    append_log(log_file, sprintf('[GROUP-FAILED] %s', MEg.message));
    append_log(log_file, getReport(MEg, 'extended', 'hyperlinks', 'off'));
end

append_log(log_file, sprintf('[FINISH] %s', datestr(now, 31)));
fprintf('\n[FINISHED] Residual FC processing complete.\n');

end

%% ========================================================================
function T = readtable_robust(fname)
try
    opts = detectImportOptions(fname, 'FileType', 'text');
    try
        opts.VariableNamingRule = 'preserve';
    catch
    end
    T = readtable(fname, opts);
catch
    T = readtable(fname);
end
end

%% ========================================================================
function v = find_varname(T, candidate)
vars = T.Properties.VariableNames;
canon_candidate = canon_name(candidate);
v = '';
for i = 1:numel(vars)
    if strcmp(canon_name(vars{i}), canon_candidate)
        v = vars{i};
        return;
    end
end
end

%% ========================================================================
function c = canon_name(s)
s = char(string(s));
s = lower(strtrim(s));
s = regexprep(s, '[^a-z0-9]+', '');
c = s;
end

%% ========================================================================
function files = get_subject_run_files(run_dir, subj, run_name, prefer_smoothed_functional)

files.func   = '';
files.motion = '';
files.art    = '';
files.wm     = '';
files.csf    = '';

func_candidates = {};
ext_patterns = {'*.nii','*.NIFTI','*.img'};

for e = 1:numel(ext_patterns)
    dd = dir(fullfile(run_dir, ext_patterns{e}));
    for i = 1:numel(dd)
        nm = dd(i).name;
        low = lower(nm);
        if contains(low, 'bold') && ...
           ~contains(low, 'mean') && ...
           ~contains(low, 'mask') && ...
           ~contains(low, 'wc1') && ...
           ~contains(low, 'wc2') && ...
           ~contains(low, 'wc3') && ...
           ~contains(low, 'art_')
            func_candidates{end+1} = fullfile(run_dir, nm); %#ok<AGROW>
        end
    end
end

if isempty(func_candidates)
    error('No functional candidates found in %s', run_dir);
end

scores = zeros(numel(func_candidates),1);
for i = 1:numel(func_candidates)
    nm = lower(func_candidates{i});
    sc = 0;
    if contains(nm, [lower(subj) '_task-onebacktask'])
        sc = sc + 50;
    end
    if contains(nm, lower(strrep(run_name,'run_','run-')))
        sc = sc + 50;
    end
    if contains(nm, 'swa')
        sc = sc + (prefer_smoothed_functional * 100);
    end
    if contains(nm, 'wa')
        sc = sc + 40;
    end
    if contains(nm, 'ua')
        sc = sc + 20;
    end
    if contains(nm, 'a')
        sc = sc + 5;
    end
    scores(i) = sc;
end

[~,ix] = max(scores);
files.func = func_candidates{ix};

motion_patterns = { ...
    sprintf('rp_%s_task-onebacktask_%s_bold.txt', subj, strrep(run_name,'run_','run-')), ...
    'rp_*bold*.txt', ...
    'rp_*.txt'};
files.motion = find_first_existing(run_dir, motion_patterns);

art_patterns = { ...
    sprintf('art_regression_outliers_and_movement_*%s*%s*bold*.mat', subj, strrep(run_name,'run_','run-')), ...
    'art_regression_outliers_and_movement_*.mat', ...
    'art_*.mat'};
files.art = find_first_existing(run_dir, art_patterns);

wm_patterns = { ...
    sprintf('wc2*%s*%s*bold*.nii', subj, strrep(run_name,'run_','run-')), ...
    sprintf('wc2*%s*%s*bold*.NIFTI', subj, strrep(run_name,'run_','run-')), ...
    'wc2*.nii', ...
    'wc2*.NIFTI', ...
    'wc2*.img'};
files.wm = find_first_existing(run_dir, wm_patterns);

csf_patterns = { ...
    sprintf('wc3*%s*%s*bold*.nii', subj, strrep(run_name,'run_','run-')), ...
    sprintf('wc3*%s*%s*bold*.NIFTI', subj, strrep(run_name,'run_','run-')), ...
    'wc3*.nii', ...
    'wc3*.NIFTI', ...
    'wc3*.img'};
files.csf = find_first_existing(run_dir, csf_patterns);

end

%% ========================================================================
function fpath = find_first_existing(run_dir, patterns)
fpath = '';
for p = 1:numel(patterns)
    dd = dir(fullfile(run_dir, patterns{p}));
    if ~isempty(dd)
        [~,ix] = sort({dd.name});
        dd = dd(ix);
        fpath = fullfile(run_dir, dd(1).name);
        return;
    end
end
end

%% ========================================================================
function X = match_rows_to_scans(X, nScans)
if isempty(X)
    return;
end
if size(X,1) ~= nScans
    nMin = min(size(X,1), nScans);
    X = X(1:nMin,:);
    if nMin < nScans
        X = [X; zeros(nScans-nMin, size(X,2))];
    end
end
end

%% ========================================================================
function art_regs = load_art_regressors(art_file, nScans, log_file)
art_regs = [];
try
    S = load(art_file);
    fns = fieldnames(S);
    for k = 1:numel(fns)
        val = S.(fns{k});
        if isnumeric(val) && size(val,1) == nScans
            art_regs = val;
            return;
        end
    end
    % fallback: any numeric matrix with approximately matching rows
    for k = 1:numel(fns)
        val = S.(fns{k});
        if isnumeric(val) && size(val,1) > 1
            art_regs = match_rows_to_scans(val, nScans);
            append_log(log_file, sprintf('[WARN] ART rows adjusted from file: %s', art_file));
            return;
        end
    end
    append_log(log_file, sprintf('[WARN] Could not identify ART regressors in %s', art_file));
catch ME
    append_log(log_file, sprintf('[WARN] Failed loading ART file %s: %s', art_file, ME.message));
end
end

%% ========================================================================
function tf = same_space(file1, file2)
try
    V1 = spm_vol(file1);
    V2 = spm_vol(file2);
    if numel(V1) > 1, V1 = V1(1); end
    if numel(V2) > 1, V2 = V2(1); end
    tf = isequal(V1.dim, V2.dim) && max(abs(V1.mat(:)-V2.mat(:))) < 1e-6;
catch
    tf = false;
end
end

%% ========================================================================
function out_roi = reslice_roi_to_func_space(roi_file, func_file, outdir, out_prefix, log_file)

Vf = spm_vol(func_file);
if numel(Vf) > 1
    ref = [func_file ',1'];
else
    ref = func_file;
end

[~,~,ext] = fileparts(roi_file);
if isempty(ext), ext = '.nii'; end
out_roi = fullfile(outdir, [out_prefix ext]);

if exist(out_roi, 'file')
    try
        if same_space(out_roi, func_file)
            append_log(log_file, sprintf('[RESLICE-SKIP] Existing resliced file reused: %s', out_roi));
            return;
        else
            delete(out_roi);
        end
    catch
        delete(out_roi);
    end
end

tmp_copy = fullfile(outdir, ['tmp_' num2str(round(now*1e8)) ext]);
copyfile(roi_file, tmp_copy);

P = char(ref, tmp_copy);
flags = struct();
flags.mask   = true;
flags.mean   = false;
flags.interp = 0;
flags.which  = 1;
flags.wrap   = [0 0 0];
flags.prefix = 'r';

append_log(log_file, sprintf('[RESLICE] %s -> %s', roi_file, out_roi));
spm_reslice(P, flags);

[pth,nm,ex] = fileparts(tmp_copy);
resliced_tmp = fullfile(pth, ['r' nm ex]);

if ~exist(resliced_tmp, 'file')
    error('Resliced ROI not created: %s', resliced_tmp);
end

movefile(resliced_tmp, out_roi, 'f');
if exist(tmp_copy, 'file'), delete(tmp_copy); end
end

%% ========================================================================
function ts = extract_mean_timeseries_from_mask(func_file, mask_file)
Vf = spm_vol(func_file);
Vm = spm_vol(mask_file);
Ym = spm_read_vols(Vm);
mask = Ym > 0.5;

if ~any(mask(:))
    warning('Mask is empty: %s', mask_file);
    ts = nan(numel(Vf),1);
    return;
end

idx = find(mask);
ts = nan(numel(Vf),1);
for t = 1:numel(Vf)
    Y = spm_read_vols(Vf(t));
    vals = Y(idx);
    vals = vals(isfinite(vals));
    if isempty(vals)
        ts(t) = NaN;
    else
        ts(t) = mean(vals);
    end
end
end

%% ========================================================================
function pcs = extract_voxelwise_pca(func_file, mask_file, ncomp)
Vf = spm_vol(func_file);
Vm = spm_vol(mask_file);
Ym = spm_read_vols(Vm);
mask = Ym > 0.5;

if ~any(mask(:))
    warning('Mask is empty for PCA: %s', mask_file);
    pcs = [];
    return;
end

idx = find(mask);
nT = numel(Vf);
X = nan(nT, numel(idx));
for t = 1:nT
    Y = spm_read_vols(Vf(t));
    X(t,:) = Y(idx)';
end

goodvox = all(isfinite(X),1) & std(X,0,1) > 0;
X = X(:,goodvox);
if isempty(X)
    warning('No valid voxels in PCA mask: %s', mask_file);
    pcs = [];
    return;
end

X = detrend(X, 'constant');
ncomp = min([ncomp, size(X,1)-1, size(X,2)]);
if ncomp < 1
    pcs = [];
    return;
end

try
    [~,score,~] = pca(X, 'NumComponents', ncomp);
catch
    [U,~,~] = svd(bsxfun(@minus, X, mean(X,1)), 'econ');
    score = U(:,1:ncomp);
end
pcs = score;
end

%% ========================================================================
function Xdct = make_dct_basis(nScans, TR, cutoff_sec)
if cutoff_sec <= 0
    Xdct = [];
    return;
end
N = nScans;
L = N * TR;
K = floor(2 * L / cutoff_sec);
if K < 1
    Xdct = [];
    return;
end
Xdct = zeros(N, K);
n = (0:N-1)';
for k = 1:K
    Xdct(:,k) = cos(pi * (2*n + 1) * k / (2*N));
end
end

%% ========================================================================
function [Xtask, names, source] = build_task_regressors_for_run(mode, task_events_dir, run_dir, subj, run_name, Tcond, nScans, TR, StartTime, label_column, task_categories, exclude_labels, oversampling, include_derivatives, log_file)

Xtask = [];
names = {};
source = 'none';
mode = lower(strtrim(mode));

if strcmp(mode, 'none')
    append_log(log_file, '[TASK] mode=none; no task regressors added.');
    return;
end

if strcmp(mode, 'auto') || strcmp(mode, 'events_tsv')
    ev_file = find_events_file(task_events_dir, run_dir, subj, run_name);
    if ~isempty(ev_file) && exist(ev_file, 'file')
        try
            [Xtask, names] = build_task_design_from_events_tsv(ev_file, nScans, TR, StartTime, label_column, task_categories, exclude_labels, oversampling, include_derivatives, log_file);
            source = ['events_tsv:' ev_file];
            return;
        catch ME
            append_log(log_file, sprintf('[TASK-WARN] Failed events TSV task design: %s', ME.message));
            if strcmp(mode, 'events_tsv')
                rethrow(ME);
            end
        end
    else
        append_log(log_file, sprintf('[TASK-WARN] No events.tsv found for %s %s.', subj, run_name));
        if strcmp(mode, 'events_tsv')
            error('task_regressor_mode=events_tsv but no events.tsv found for %s %s.', subj, run_name);
        end
    end
end

if strcmp(mode, 'auto') || strcmp(mode, 'block_csv')
    try
        session_num = parse_run_number(run_name);
        Tsub = get_subject_condition_rows(Tcond, subj, session_num);
        if isempty(Tsub)
            append_log(log_file, '[TASK-WARN] No block rows for block_csv task design.');
            return;
        end
        [Xtask, names] = build_task_design_from_block_csv(Tsub, nScans, TR, StartTime, oversampling, include_derivatives);
        source = 'block_csv';
        return;
    catch ME
        append_log(log_file, sprintf('[TASK-WARN] Failed block_csv task design: %s', ME.message));
        if strcmp(mode, 'block_csv')
            rethrow(ME);
        end
    end
end

end

%% ========================================================================
function ev_file = find_events_file(task_events_dir, run_dir, subj, run_name)
ev_file = '';
run_dash = strrep(run_name, 'run_', 'run-');
run_num = parse_run_number(run_name);
patterns = { ...
    sprintf('%s*task-onebacktask*%s*events.tsv', subj, run_dash), ...
    sprintf('%s*%s*events.tsv', subj, run_dash), ...
    sprintf('*task-onebacktask*%s*events.tsv', run_dash), ...
    sprintf('*run-%02d*events.tsv', run_num), ...
    sprintf('*run-%d*events.tsv', run_num), ...
    '*events.tsv'};

search_dirs = {run_dir, fileparts(run_dir), fullfile(task_events_dir, subj, 'func'), task_events_dir};

for d = 1:numel(search_dirs)
    if isempty(search_dirs{d}) || ~exist(search_dirs{d}, 'dir'), continue; end
    for p = 1:numel(patterns)
        dd = dir(fullfile(search_dirs{d}, patterns{p}));
        if ~isempty(dd)
            [~,ix] = sort({dd.name});
            dd = dd(ix);
            ev_file = fullfile(search_dirs{d}, dd(1).name);
            return;
        end
    end
end

% recursive fallback, limited to subject folder
subj_root = fullfile(task_events_dir, subj);
if exist(subj_root, 'dir')
    dd = dir(fullfile(subj_root, '**', '*events.tsv'));
    if ~isempty(dd)
        names = {dd.name};
        keep = contains(lower(names), lower(run_dash)) | contains(lower(names), sprintf('run-%02d',run_num));
        dd = dd(keep);
        if ~isempty(dd)
            [~,ix] = sort({dd.name});
            dd = dd(ix);
            ev_file = fullfile(dd(1).folder, dd(1).name);
        end
    end
end
end

%% ========================================================================
function [Xtask, names] = build_task_design_from_events_tsv(ev_file, nScans, TR, StartTime, label_column, task_categories, exclude_labels, oversampling, include_derivatives, log_file)

T = readtable_robust(ev_file);
onset_var = find_varname(T, 'onset');
duration_var = find_varname(T, 'duration');
if isempty(onset_var) || isempty(duration_var)
    error('events.tsv must contain onset and duration columns: %s', ev_file);
end

if strcmpi(label_column, 'auto')
    label_candidates = {'trial_type','condition','stim_type','stimulus_type','stimulus','category','event_type'};
    label_var = '';
    for k = 1:numel(label_candidates)
        label_var = find_varname(T, label_candidates{k});
        if ~isempty(label_var), break; end
    end
    if isempty(label_var)
        % If no category column exists, create one all-task regressor.
        T.all_task_label = repmat({'task'}, height(T), 1);
        label_var = 'all_task_label';
    end
else
    label_var = find_varname(T, label_column);
    if isempty(label_var)
        error('Requested event_label_column not found: %s', label_column);
    end
end

labels = table_column_to_cellstr(T.(label_var));
labels = strtrim(labels);
labels_lower = lower(labels);

% exclude labels such as fixation/rest/baseline
exclude_lower = lower(cellfun(@char, exclude_labels, 'UniformOutput', false));
keep = true(numel(labels),1);
for k = 1:numel(exclude_lower)
    keep = keep & ~strcmp(labels_lower, exclude_lower{k});
end

onsets = double(T.(onset_var));
durs   = double(T.(duration_var));
keep = keep & isfinite(onsets) & isfinite(durs);

labels = labels(keep);
onsets = onsets(keep);
durs   = durs(keep);

if isempty(labels)
    Xtask = [];
    names = {};
    append_log(log_file, sprintf('[TASK-WARN] No usable events in %s', ev_file));
    return;
end

if isempty(task_categories)
    cats = unique(labels, 'stable');
else
    cats = task_categories(:)';
end

nCat = numel(cats);
Xtask = [];
names = {};
for c = 1:nCat
    cat = char(string(cats{c}));
    idx = strcmpi(labels, cat);
    if ~any(idx)
        append_log(log_file, sprintf('[TASK-WARN] Category not found in events file: %s', cat));
        continue;
    end
    x = events_to_hrf_regressor(onsets(idx), durs(idx), nScans, TR, StartTime, oversampling);
    if std(x) > 0
        Xtask = [Xtask, x]; %#ok<AGROW>
        names{end+1} = ['task_' matlab.lang.makeValidName(cat)]; %#ok<AGROW>
        if include_derivatives
            xd = [0; diff(x)];
            if std(xd) > 0
                Xtask = [Xtask, xd]; %#ok<AGROW>
                names{end+1} = ['taskDeriv_' matlab.lang.makeValidName(cat)]; %#ok<AGROW>
            end
        end
    end
end

append_log(log_file, sprintf('[TASK] events file used: %s | n categories modeled=%d', ev_file, numel(names)));
end

%% ========================================================================
function [Xtask, names] = build_task_design_from_block_csv(Tsub, nScans, TR, StartTime, oversampling, include_derivatives)
cond_var = find_varname(Tsub, 'condition');
onset_var = find_varname(Tsub, 'onset');
dur_var = find_varname(Tsub, 'duration');

cond_values = unique(Tsub.(cond_var));
cond_values = cond_values(~isnan(cond_values));
Xtask = [];
names = {};

for cc = 1:numel(cond_values)
    cond_num = cond_values(cc);
    rows = Tsub(Tsub.(cond_var) == cond_num, :);
    x = events_to_hrf_regressor(double(rows.(onset_var)), double(rows.(dur_var)), nScans, TR, StartTime, oversampling);
    if std(x) > 0
        Xtask = [Xtask, x]; %#ok<AGROW>
        names{end+1} = sprintf('task_Block%02d', cond_num); %#ok<AGROW>
        if include_derivatives
            xd = [0; diff(x)];
            if std(xd) > 0
                Xtask = [Xtask, xd]; %#ok<AGROW>
                names{end+1} = sprintf('taskDeriv_Block%02d', cond_num); %#ok<AGROW>
            end
        end
    end
end
end

%% ========================================================================
function x = events_to_hrf_regressor(onsets, durs, nScans, TR, StartTime, oversampling)
if nargin < 6 || isempty(oversampling) || oversampling < 1
    oversampling = 16;
end

dt = TR / oversampling;
total_time = nScans * TR;
nt = ceil(total_time / dt) + 1;
u = zeros(nt,1);

for k = 1:numel(onsets)
    t0 = onsets(k) - StartTime;
    dur = durs(k);
    if ~isfinite(dur) || dur < 0
        dur = 0;
    end
    t0 = max(0, t0);
    t1 = min(total_time, t0 + dur);

    if dur <= 0
        ix = round(t0 / dt) + 1;
        ix = min(max(ix,1),nt);
        u(ix) = u(ix) + 1;
    else
        ix0 = floor(t0 / dt) + 1;
        ix1 = max(ix0, ceil(t1 / dt));
        ix0 = min(max(ix0,1),nt);
        ix1 = min(max(ix1,1),nt);
        u(ix0:ix1) = u(ix0:ix1) + 1;
    end
end

hrf = spm_hrf(dt);
y = conv(u, hrf);
y = y(1:nt);

scan_times = (0:nScans-1)' * TR;
fine_times = (0:nt-1)' * dt;
x = interp1(fine_times, y, scan_times, 'linear', 0);

% demean but do not zscore; scaling does not affect regression residuals when
% using ordinary least squares, but demeaning helps numerical conditioning.
x = x(:);
x = x - mean(x(isfinite(x)));
if any(~isfinite(x))
    x(~isfinite(x)) = 0;
end
end

%% ========================================================================
function c = table_column_to_cellstr(x)
if iscell(x)
    c = cell(size(x));
    for i = 1:numel(x)
        c{i} = char(string(x{i}));
    end
elseif isstring(x)
    c = cellstr(x);
elseif ischar(x)
    c = cellstr(x);
elseif isnumeric(x)
    c = cell(size(x));
    for i = 1:numel(x)
        c{i} = num2str(x(i));
    end
elseif iscategorical(x)
    c = cellstr(x);
else
    c = cell(size(x));
    for i = 1:numel(x)
        c{i} = char(string(x(i)));
    end
end
c = c(:);
end

%% ========================================================================
function [Xclean, names_clean, removed_cols] = clean_design_matrix(X, names, log_file)
if isempty(X)
    Xclean = [];
    names_clean = {};
    removed_cols = {};
    return;
end

% Replace non-finite values with zeros after warning.
if any(~isfinite(X(:)))
    append_log(log_file, '[WARN] Non-finite values in design matrix replaced with zero.');
    X(~isfinite(X)) = 0;
end

bad = false(1,size(X,2));
for c = 1:size(X,2)
    if all(abs(X(:,c)) < eps) || std(X(:,c)) == 0
        % keep intercept column if it is the first column
        if c ~= 1
            bad(c) = true;
        end
    end
end

removed_cols = names(bad);
X(:,bad) = [];
names(bad) = [];

% Remove near-duplicate/rank-deficient columns conservatively using QR.
try
    [~,R,E] = qr(X,0);
    tol = max(size(X)) * eps(norm(R,'fro')) * 100;
    rnk = sum(abs(diag(R)) > tol);
    if rnk < size(X,2)
        keep_idx = sort(E(1:rnk));
        drop_idx = setdiff(1:size(X,2), keep_idx);
        removed_cols = [removed_cols, names(drop_idx)]; %#ok<AGROW>
        X = X(:,keep_idx);
        names = names(keep_idx);
        append_log(log_file, sprintf('[WARN] Removed %d rank-deficient design columns.', numel(drop_idx)));
    end
catch
    append_log(log_file, '[WARN] QR rank check failed; proceeding without rank pruning.');
end

Xclean = X;
names_clean = names;
end

%% ========================================================================
function Yclean = residualize_timeseries(Yraw, X, do_zscore)
Yclean = nan(size(Yraw));
for c = 1:size(Yraw,2)
    y = Yraw(:,c);
    idx = isfinite(y) & all(isfinite(X),2);
    yc = nan(size(y));
    if sum(idx) > size(X,2) + 2
        beta = X(idx,:) \ y(idx);
        yc(idx) = y(idx) - X(idx,:)*beta;
        if do_zscore
            mu = mean(yc(idx));
            sd = std(yc(idx));
            if sd > 0
                yc(idx) = (yc(idx) - mu) ./ sd;
            end
        end
    end
    Yclean(:,c) = yc;
end
end

%% ========================================================================
function session_num = parse_run_number(run_name)
tok = regexp(run_name, 'run_(\d+)', 'tokens', 'once');
if isempty(tok)
    tok = regexp(run_name, 'run-(\d+)', 'tokens', 'once');
end
if isempty(tok)
    error('Could not parse session number from run name: %s', run_name);
end
session_num = str2double(tok{1});
end

%% ========================================================================
function Tsub = get_subject_condition_rows(Tcond, subj_label, session_num)
subject_var = find_varname(Tcond, 'subject');
real_subject_var = find_varname(Tcond, 'real subject');
session_var = find_varname(Tcond, 'session');

subnum = sscanf(subj_label, 'sub-%d');
cand_rows = false(height(Tcond),1);

if ~isempty(subject_var)
    cand_rows = cand_rows | match_subject_column(Tcond.(subject_var), subj_label, subnum);
end
if ~isempty(real_subject_var)
    cand_rows = cand_rows | match_subject_column(Tcond.(real_subject_var), subj_label, subnum);
end

sess_ok = match_session_column(Tcond.(session_var), session_num);
Tsub = Tcond(cand_rows & sess_ok, :);
end

%% ========================================================================
function tf = match_subject_column(col, subj_label, subnum)
tf = false(numel(col),1);
if isnumeric(col)
    tf = (col == subnum);
else
    txt = table_column_to_cellstr(col);
    tf = strcmp(txt, subj_label) | strcmp(txt, sprintf('%d',subnum)) | strcmp(txt, sprintf('sub-%02d',subnum));
end
end

%% ========================================================================
function tf = match_session_column(col, session_num)
tf = false(numel(col),1);
if isnumeric(col)
    tf = (col == session_num);
else
    txt = table_column_to_cellstr(col);
    tf = strcmp(txt, sprintf('%d',session_num)) | strcmp(txt, sprintf('%.0f',session_num));
end
end

%% ========================================================================
function block_idx = make_block_scan_index(rows, nScans, TR, StartTime)
onset_var = find_varname(rows, 'onset');
dur_var = find_varname(rows, 'duration');
block_idx = false(nScans,1);
for k = 1:height(rows)
    onset_sec = double(rows.(onset_var)(k)) - StartTime;
    dur_sec   = double(rows.(dur_var)(k));
    t0 = max(0, onset_sec);
    t1 = max(t0, onset_sec + dur_sec);
    scan_idx = ((0:nScans-1)' * TR >= t0) & ((0:nScans-1)' * TR < t1);
    block_idx = block_idx | scan_idx;
end
end

%% ========================================================================
function compute_and_save_fc(Y, roi_names, outdir, tag, make_figs, min_valid_timepoints_fc, low_variance_threshold, clip_r_for_fisher, subj, run_name, log_file)

if ~exist(outdir,'dir')
    mkdir(outdir);
end

nrois = size(Y,2);
R = nan(nrois, nrois);

for i = 1:nrois
    yi = Y(:,i);
    for j = i:nrois
        yj = Y(:,j);
        idx = isfinite(yi) & isfinite(yj);
        if sum(idx) < min_valid_timepoints_fc
            r = NaN;
        elseif std(yi(idx)) < low_variance_threshold || std(yj(idx)) < low_variance_threshold
            r = NaN;
        else
            C = corrcoef(yi(idx), yj(idx));
            r = C(1,2);
        end
        R(i,j) = r;
        R(j,i) = r;
    end
end

R(1:nrois+1:end) = 0;
Rclip = min(max(R, -clip_r_for_fisher), clip_r_for_fisher);
Z = atanh(Rclip);
Z(1:nrois+1:end) = 0;

r_cell = cell(nrois+1, nrois+1);
z_cell = cell(nrois+1, nrois+1);
r_cell{1,1} = 'ROI';
z_cell{1,1} = 'ROI';
for i = 1:nrois
    r_cell{1,i+1} = roi_names{i};
    z_cell{1,i+1} = roi_names{i};
    r_cell{i+1,1} = roi_names{i};
    z_cell{i+1,1} = roi_names{i};
end
for i = 1:nrois
    for j = 1:nrois
        r_cell{i+1,j+1} = R(i,j);
        z_cell{i+1,j+1} = Z(i,j);
    end
end

writecell(r_cell, fullfile(outdir, ['ResidualFC_' tag '_r.csv']));
writecell(z_cell, fullfile(outdir, ['ResidualFC_' tag '_z.csv']));

vec_names = cell(0,1);
r_vec = [];
z_vec = [];
for i = 1:nrois-1
    for j = i+1:nrois
        vec_names{end+1,1} = sprintf('%s__%s', roi_names{i}, roi_names{j}); %#ok<AGROW>
        r_vec(end+1,1) = R(i,j); %#ok<AGROW>
        z_vec(end+1,1) = Z(i,j); %#ok<AGROW>
    end
end

edge_tbl = table(vec_names, r_vec, z_vec, 'VariableNames', {'edge','r','z'});
writetable(edge_tbl, fullfile(outdir, ['ResidualFC_' tag '_zvec.csv']));

save(fullfile(outdir, ['ResidualFC_' tag '.mat']), 'R','Z','roi_names','tag','subj','run_name');

if make_figs
    try
        f1 = figure('Visible','off','Color','w','Position',[100 100 700 620]);
        imagesc(R);
        axis square;
        colorbar;
        caxis([-1 1]);
        title(['Residual FC r - ' tag], 'Interpreter', 'none');
        set(gca,'XTick',1:nrois,'XTickLabel',roi_names,'XTickLabelRotation',45);
        set(gca,'YTick',1:nrois,'YTickLabel',roi_names);
        saveas(f1, fullfile(outdir, ['ResidualFC_' tag '_r_heatmap.png']));
        close(f1);

        f2 = figure('Visible','off','Color','w','Position',[100 100 700 620]);
        imagesc(Z);
        axis square;
        colorbar;
        title(['Residual FC z - ' tag], 'Interpreter', 'none');
        set(gca,'XTick',1:nrois,'XTickLabel',roi_names,'XTickLabelRotation',45);
        set(gca,'YTick',1:nrois,'YTickLabel',roi_names);
        saveas(f2, fullfile(outdir, ['ResidualFC_' tag '_z_heatmap.png']));
        close(f2);

        validz = z_vec(isfinite(z_vec));
        f3 = figure('Visible','off','Color','w','Position',[100 100 700 500]);
        histogram(validz, 20);
        title(['Residual FC z histogram - ' tag], 'Interpreter', 'none');
        xlabel('Fisher z');
        ylabel('Count');
        saveas(f3, fullfile(outdir, ['ResidualFC_' tag '_z_hist.png']));
        close(f3);
    catch figME
        warning('Figure saving failed for %s: %s', tag, figME.message);
        append_log(log_file, sprintf('[WARN] Figure saving failed for %s: %s', tag, figME.message));
    end
end

append_log(log_file, sprintf('[SAVED] %s | %s | %s', subj, run_name, tag));
fprintf('[SAVED] %s\n', fullfile(outdir, ['ResidualFC_' tag '.mat']));
end

%% ========================================================================
function export_group_level_fc(output_root, roi_names, log_file)

group_dir = fullfile(output_root, 'ResidualFC_group_exports');
if ~exist(group_dir, 'dir')
    mkdir(group_dir);
end

subdirs = dir(fullfile(output_root, 'sub-*'));
subdirs = subdirs([subdirs.isdir]);

if isempty(subdirs)
    append_log(log_file, '[GROUP] No subject output folders found.');
    return;
end

run_labels = {'run_1','run_2'};
tag_labels = {'ALL','Block-01','Block-02','Block-03'};

allsubject_rows = {};
allsubject_header_done = false;
edge_summary_rows = table();
edge_summary_header_done = false;

for rr = 1:numel(run_labels)
    run_name = run_labels{rr};
    run_num = parse_run_number(run_name);

    for tt = 1:numel(tag_labels)
        tag = tag_labels{tt};

        mats = {};
        subj_list = {};
        edge_table_this = table();

        for s = 1:numel(subdirs)
            subj = subdirs(s).name;
            matfile = fullfile(output_root, subj, run_name, ['ResidualFC_' tag '.mat']);
            if exist(matfile, 'file')
                S = load(matfile);
                if isfield(S, 'Z')
                    mats{end+1} = S.Z; %#ok<AGROW>
                    subj_list{end+1} = subj; %#ok<AGROW>
                    zvec_file = fullfile(output_root, subj, run_name, ['ResidualFC_' tag '_zvec.csv']);
                    if exist(zvec_file, 'file')
                        T = readtable(zvec_file);
                        T.subject = repmat({subj}, height(T), 1);
                        T.run = repmat({run_name}, height(T), 1);
                        T.tag = repmat({tag}, height(T), 1);
                        if isempty(edge_table_this)
                            edge_table_this = T;
                        else
                            edge_table_this = [edge_table_this; T]; %#ok<AGROW>
                        end
                    end
                end
            end
        end

        nSub = numel(mats);
        append_log(log_file, sprintf('[GROUP] %s %s subjects=%d', run_name, tag, nSub));
        if nSub == 0, continue; end

        nroi = size(mats{1},1);
        Zstack = nan(nroi, nroi, nSub);
        for k = 1:nSub
            Zstack(:,:,k) = mats{k};
        end

        Zmean = nanmean_compat(Zstack, 3);
        Zmedian = nanmedian_compat(Zstack, 3);
        Zmode = matrix_mode_3d(Zstack);

        write_matrix_csv_with_labels(fullfile(group_dir, sprintf('ResidualFC_group_R%d_%s_mean.csv', run_num, tag)), Zmean, roi_names);
        write_matrix_csv_with_labels(fullfile(group_dir, sprintf('ResidualFC_group_R%d_%s_median.csv', run_num, tag)), Zmedian, roi_names);
        write_matrix_csv_with_labels(fullfile(group_dir, sprintf('ResidualFC_group_R%d_%s_mode.csv', run_num, tag)), Zmode, roi_names);

        save(fullfile(group_dir, sprintf('ResidualFC_group_R%d_%s_stack.mat', run_num, tag)), ...
            'Zstack','Zmean','Zmedian','Zmode','roi_names','subj_list','run_name','tag');

        if ~isempty(edge_table_this)
            edge_names = unique(edge_table_this.edge, 'stable');
            summary_edge = table();
            edge_col = cell(numel(edge_names),1);
            mean_col = nan(numel(edge_names),1);
            median_col = nan(numel(edge_names),1);
            mode_col = nan(numel(edge_names),1);
            sd_col = nan(numel(edge_names),1);
            sem_col = nan(numel(edge_names),1);
            run_col = cell(numel(edge_names),1);
            tag_col = cell(numel(edge_names),1);
            nsub_col = nan(numel(edge_names),1);

            for e = 1:numel(edge_names)
                edge_name = edge_names{e};
                idx = strcmp(edge_table_this.edge, edge_name);
                vals = edge_table_this.z(idx);
                vals = vals(isfinite(vals));
                edge_col{e} = edge_name;
                run_col{e} = run_name;
                tag_col{e} = tag;
                nsub_col(e) = numel(vals);
                if ~isempty(vals)
                    mean_col(e) = mean(vals);
                    median_col(e) = median(vals);
                    mode_col(e) = safe_mode(vals);
                    sd_col(e) = std(vals);
                    sem_col(e) = std(vals) ./ sqrt(numel(vals));
                end
            end

            summary_edge.edge = edge_col;
            summary_edge.run = run_col;
            summary_edge.tag = tag_col;
            summary_edge.n = nsub_col;
            summary_edge.mean_z = mean_col;
            summary_edge.median_z = median_col;
            summary_edge.mode_z = mode_col;
            summary_edge.sd_z = sd_col;
            summary_edge.sem_z = sem_col;
            writetable(summary_edge, fullfile(group_dir, sprintf('ResidualFC_group_R%d_%s_edge_summary.csv', run_num, tag)));

            if ~edge_summary_header_done
                edge_summary_rows = summary_edge;
                edge_summary_header_done = true;
            else
                edge_summary_rows = [edge_summary_rows; summary_edge]; %#ok<AGROW>
            end
        end

        [wideTbl, headerCells] = make_allsubjects_wide_table(subj_list, mats, roi_names, run_name, tag);
        writetable(wideTbl, fullfile(group_dir, sprintf('ResidualFC_group_R%d_%s_allsubjects_wide.csv', run_num, tag)));

        if ~allsubject_header_done
            allsubject_rows = headerCells;
            allsubject_header_done = true;
        else
            allsubject_rows = [allsubject_rows; headerCells(2:end,:)]; %#ok<AGROW>
        end
    end
end

if edge_summary_header_done
    writetable(edge_summary_rows, fullfile(group_dir, 'ResidualFC_group_edge_summary.csv'));
end
if allsubject_header_done
    writecell(allsubject_rows, fullfile(group_dir, 'ResidualFC_allsubjects_wide.csv'));
end

inventory = build_fc_inventory(output_root);
writetable(inventory, fullfile(group_dir, 'ResidualFC_file_inventory.csv'));
append_log(log_file, sprintf('[GROUP] Outputs saved under %s', group_dir));
end

%% ========================================================================
function write_matrix_csv_with_labels(outfile, M, roi_names)
nrois = numel(roi_names);
C = cell(nrois+1, nrois+1);
C{1,1} = 'ROI';
for i = 1:nrois
    C{1,i+1} = roi_names{i};
    C{i+1,1} = roi_names{i};
end
for i = 1:nrois
    for j = 1:nrois
        C{i+1,j+1} = M(i,j);
    end
end
writecell(C, outfile);
end

%% ========================================================================
function Mmode = matrix_mode_3d(Zstack)
[n1,n2,~] = size(Zstack);
Mmode = nan(n1,n2);
for i = 1:n1
    for j = 1:n2
        vals = squeeze(Zstack(i,j,:));
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            Mmode(i,j) = safe_mode(vals);
        end
    end
end
end

%% ========================================================================
function m = safe_mode(vals)
vals = vals(:);
vals = vals(isfinite(vals));
if isempty(vals)
    m = NaN;
    return;
end
vals_round = round(vals, 6);
u = unique(vals_round);
counts = zeros(size(u));
for i = 1:numel(u)
    counts(i) = sum(vals_round == u(i));
end
[~,ix] = max(counts);
m = u(ix);
end

%% ========================================================================
function [wideTbl, cellOut] = make_allsubjects_wide_table(subj_list, mats, roi_names, run_name, tag)
edge_names = {};
for i = 1:numel(roi_names)-1
    for j = i+1:numel(roi_names)
        edge_names{end+1,1} = sprintf('%s__%s', roi_names{i}, roi_names{j}); %#ok<AGROW>
    end
end

nSub = numel(subj_list);
nEdge = numel(edge_names);
subject_col = cell(nSub,1);
run_col = cell(nSub,1);
tag_col = cell(nSub,1);
data = nan(nSub, nEdge);

for s = 1:nSub
    subject_col{s} = subj_list{s};
    run_col{s} = run_name;
    tag_col{s} = tag;
    Z = mats{s};
    c = 0;
    for i = 1:numel(roi_names)-1
        for j = i+1:numel(roi_names)
            c = c + 1;
            data(s,c) = Z(i,j);
        end
    end
end

wideTbl = table(subject_col, run_col, tag_col, 'VariableNames', {'subject','run','tag'});
for c = 1:nEdge
    vname = matlab.lang.makeValidName(edge_names{c});
    wideTbl.(vname) = data(:,c);
end
cellOut = [wideTbl.Properties.VariableNames; table2cell(wideTbl)];
end

%% ========================================================================
function T = build_fc_inventory(output_root)
subdirs = dir(fullfile(output_root, 'sub-*'));
subdirs = subdirs([subdirs.isdir]);
subject_col = {};
run_col = {};
tag_col = {};
exists_col = [];
path_col = {};
tags = {'ALL','Block-01','Block-02','Block-03'};
runs = {'run_1','run_2'};
for s = 1:numel(subdirs)
    subj = subdirs(s).name;
    for r = 1:numel(runs)
        for t = 1:numel(tags)
            f = fullfile(output_root, subj, runs{r}, ['ResidualFC_' tags{t} '.mat']);
            subject_col{end+1,1} = subj; %#ok<AGROW>
            run_col{end+1,1} = runs{r}; %#ok<AGROW>
            tag_col{end+1,1} = tags{t}; %#ok<AGROW>
            exists_col(end+1,1) = exist(f, 'file') ~= 0; %#ok<AGROW>
            path_col{end+1,1} = f; %#ok<AGROW>
        end
    end
end
T = table(subject_col, run_col, tag_col, exists_col, path_col, ...
    'VariableNames', {'subject','run','tag','exists','path'});
end

%% ========================================================================
function M = nanmean_compat(X, dim)
try
    M = mean(X, dim, 'omitnan');
catch
    good = ~isnan(X);
    X0 = X;
    X0(~good) = 0;
    n = sum(good, dim);
    M = sum(X0, dim) ./ n;
    M(n==0) = NaN;
end
end

%% ========================================================================
function M = nanmedian_compat(X, dim)
try
    M = median(X, dim, 'omitnan');
catch
    M = nanmedian(X, dim);
end
end

%% ========================================================================
function append_log(log_file, msg)
fid = fopen(log_file, 'a');
if fid == -1
    fprintf(2, '[LOG-FAIL] %s\n', msg);
    return;
end
fprintf(fid, '%s | %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg);
fclose(fid);
end
