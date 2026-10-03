function run_raw_fc_analysis()
clc;
clear;
close all;

%% =========================
% PATH CONFIGURATION
% ==========================
cfg = project_paths();
base_dir = cfg.data_root;
output_root = cfg.raw_fc_root;
cond_csv = cfg.block_timing_csv;
json_file = cfg.task_json;
roi_dir = cfg.roi_root;

roi_files = { ...
    'IFG_L_mask.nii', ...
    'lOTC_mask.nii', ...
    'V1_L_mask.nii', ...
    'vOTC_mask.nii'};

roi_names = { ...
    'IFG_L', ...
    'lOTC', ...
    'V1_L', ...
    'vOTC'};

make_figs = true;
prefer_smoothed_functional = true;
do_zscore_after_denoise = true;
use_dct_highpass = true;
highpass_cutoff_sec = 128;
wm_n_components = 5;
csf_n_components = 5;
min_valid_timepoints_fc = 10;
low_variance_threshold = 1e-8;
clip_r_for_fisher = 0.999999;

% Optional subset selection. Empty lists process all available subjects and runs.
restrict_subjects = {};
restrict_runs     = {};

%% =========================
% INIT
% ==========================
if ~exist(output_root, 'dir')
    mkdir(output_root);
end

log_file = fullfile(output_root, 'FC_processing_log.txt');
append_log(log_file, '=====================================================');
append_log(log_file, sprintf('[START] %s', datestr(now, 31)));
append_log(log_file, sprintf('[INFO] output_root = %s', output_root));

if isempty(which('spm'))
    error('SPM12 is not on MATLAB path.');
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
% LOAD CONDITION CSV
% ==========================
if ~exist(cond_csv, 'file')
    error('Condition CSV not found: %s', cond_csv);
end

Tcond = readtable(cond_csv);

required_cols = {'subject','session','condition','onset','duration'};
for i = 1:numel(required_cols)
    if ~ismember(required_cols{i}, Tcond.Properties.VariableNames)
        error('Missing required column in FC_blocks_COND.csv: %s', required_cols{i});
    end
end

has_real_subject = ismember('real subject', Tcond.Properties.VariableNames);

%% =========================
% FIND SUBJECTS
% ==========================
d = dir(fullfile(base_dir, 'sub-*'));
d = d([d.isdir]);
subjects = {d.name};

if isempty(subjects)
    error('No subject folders found under: %s', base_dir);
end

subjects = sort(subjects(:));
fprintf('[INFO] Found %d subjects.\n', numel(subjects));
disp(subjects);
append_log(log_file, sprintf('[INFO] Found %d subject folders.', numel(subjects)));

if ~isempty(restrict_subjects)
    subjects = intersect(subjects, restrict_subjects, 'stable');
    fprintf('[INFO] Restricted subjects:\n');
    disp(subjects);
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
    run_names = {run_dirs.name};

    if isempty(run_names)
        warning('No run_* directories for %s. Skipping.', subj);
        append_log(log_file, sprintf('[WARN] No run_* directories for %s. Skipping.', subj));
        continue;
    end

    run_names = sort(run_names(:));

    if ~isempty(restrict_runs)
        run_names = intersect(run_names, restrict_runs, 'stable');
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

            files = get_subject_run_files(run_dir, subj, run_name, prefer_smoothed_functional);

            fprintf('\n[FILE MATCH]\n');
            fprintf('  Chosen functional: %s\n', files.func);
            fprintf('  Chosen motion:     %s\n', files.motion);
            fprintf('  Chosen ART:        %s\n', files.art);
            fprintf('  Chosen WM:         %s\n', files.wm);
            fprintf('  Chosen CSF:        %s\n', files.csf);

            append_log(log_file, sprintf('[FILES] func=%s', files.func));
            append_log(log_file, sprintf('[FILES] motion=%s', files.motion));
            append_log(log_file, sprintf('[FILES] art=%s', files.art));
            append_log(log_file, sprintf('[FILES] wm=%s', files.wm));
            append_log(log_file, sprintf('[FILES] csf=%s', files.csf));

            if isempty(files.func) || ~exist(files.func, 'file')
                error('Functional file not found.');
            end

            fprintf('[INFO] Functional file: %s\n', files.func);

            Vfunc = spm_vol(files.func);
            nScans = numel(Vfunc);
            fprintf('[INFO] Number of scans: %d\n', nScans);
            append_log(log_file, sprintf('[INFO] %s %s nScans=%d', subj, run_name, nScans));

            % ===== load motion =====
            motion = [];
            if ~isempty(files.motion) && exist(files.motion, 'file')
                motion = load(files.motion);
                if size(motion,1) ~= nScans
                    warning('Motion rows (%d) != nScans (%d). Using min rows.', size(motion,1), nScans);
                    append_log(log_file, sprintf('[WARN] Motion rows mismatch: %d vs %d', size(motion,1), nScans));
                    nMin = min(size(motion,1), nScans);
                    motion = motion(1:nMin,:);
                end
            end

            % ===== load ART =====
            art_regs = [];
            if ~isempty(files.art) && exist(files.art, 'file')
                S = load(files.art);
                fns = fieldnames(S);
                found_mat = false;
                for k = 1:numel(fns)
                    val = S.(fns{k});
                    if isnumeric(val) && size(val,1) == nScans
                        art_regs = val;
                        found_mat = true;
                        break;
                    end
                end
                if ~found_mat
                    warning('Could not identify ART regressors in: %s', files.art);
                    append_log(log_file, sprintf('[WARN] Could not identify ART regressors in %s', files.art));
                end
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
                else
                    append_log(log_file, sprintf('[WM] same-space used directly: %s', files.wm));
                end
                wm_ts = extract_voxelwise_pca(files.func, wm_use, wm_n_components);
            else
                warning('WM file missing.');
                append_log(log_file, '[WARN] WM file missing.');
            end

            if ~isempty(files.csf) && exist(files.csf, 'file')
                csf_use = files.csf;
                if ~same_space(files.csf, files.func)
                    csf_use = reslice_roi_to_func_space(files.csf, files.func, run_outdir, 'csf_resliced', log_file);
                else
                    append_log(log_file, sprintf('[CSF] same-space used directly: %s', files.csf));
                end
                csf_ts = extract_voxelwise_pca(files.func, csf_use, csf_n_components);
            else
                warning('CSF file missing.');
                append_log(log_file, '[WARN] CSF file missing.');
            end

            % ===== build nuisance design =====
            X = ones(nScans,1);

            if ~isempty(motion)
                nUse = min(size(motion,1), nScans);
                motion = motion(1:nUse,:);
                if nUse < nScans
                    motion = [motion; zeros(nScans-nUse,size(motion,2))];
                end
                X = [X, motion];
            end

            if ~isempty(art_regs)
                nUse = min(size(art_regs,1), nScans);
                art_regs = art_regs(1:nUse,:);
                if nUse < nScans
                    art_regs = [art_regs; zeros(nScans-nUse,size(art_regs,2))];
                end
                X = [X, art_regs];
            end

            if ~isempty(wm_ts)
                X = [X, wm_ts];
            end

            if ~isempty(csf_ts)
                X = [X, csf_ts];
            end

            if use_dct_highpass
                Xhp = make_dct_basis(nScans, TR, highpass_cutoff_sec);
                X = [X, Xhp];
            end

            badX = any(~isfinite(X),1) | (std(X,0,1)==0);
            X = X(:, ~badX);

            % ===== denoise ROI time series =====
            Yclean = nan(size(Yraw));
            for c = 1:size(Yraw,2)
                y = Yraw(:,c);
                idx = isfinite(y) & all(isfinite(X),2);
                yc = nan(size(y));

                if sum(idx) > size(X,2) + 2
                    beta = X(idx,:) \ y(idx);
                    yc(idx) = y(idx) - X(idx,:)*beta;
                    if do_zscore_after_denoise
                        mu = mean(yc(idx));
                        sd = std(yc(idx));
                        if sd > 0
                            yc(idx) = (yc(idx) - mu) ./ sd;
                        end
                    end
                end
                Yclean(:,c) = yc;
            end

            clean_tbl = array2table(Yclean, 'VariableNames', matlab.lang.makeValidName(roi_names));
            writetable(clean_tbl, fullfile(run_outdir, 'roi_timeseries_clean.csv'));

            save(fullfile(run_outdir, 'run_workspace.mat'), ...
                'subj','run_name','files','TR','StartTime','roi_names','roi_files','roi_used_paths', ...
                'Yraw','Yclean','motion','art_regs','wm_ts','csf_ts','X');

            % ===== FC ALL =====
            disp('--- Raw FC ALL-level diagnostic ---');
            disp(size(Yclean));
            disp(sum(isfinite(Yclean)));
            disp('Calling compute_and_save_fc for ALL...');

            compute_and_save_fc(Yclean, roi_names, run_outdir, 'ALL', ...
                make_figs, min_valid_timepoints_fc, low_variance_threshold, clip_r_for_fisher, ...
                subj, run_name, log_file);

            % ===== BLOCK-WISE FC =====
            session_num = parse_run_number(run_name);
            Tsub = get_subject_condition_rows(Tcond, subj, session_num, has_real_subject);

            disp('--- CONDITION ROWS FOUND ---');
            try
                disp(Tsub(:, {'subject','session','condition','onset','duration'}));
            catch
                disp(Tsub);
            end

            if isempty(Tsub)
                warning('No condition rows found for %s session %d', subj, session_num);
                append_log(log_file, sprintf('[WARN] No condition rows found for %s session %d', subj, session_num));
            else
                cond_values = unique(Tsub.condition);
                cond_values = cond_values(~isnan(cond_values));

                for cc = 1:numel(cond_values)
                    cond_num = cond_values(cc);
                    rows = Tsub(Tsub.condition == cond_num, :);

                    block_idx = false(nScans,1);
                    for k = 1:height(rows)
                        onset_sec = rows.onset(k) - StartTime;
                        dur_sec   = rows.duration(k);

                        t0 = max(0, onset_sec);
                        t1 = max(t0, onset_sec + dur_sec);

                        scan_idx = ((0:nScans-1)' * TR >= t0) & ((0:nScans-1)' * TR < t1);
                        block_idx = block_idx | scan_idx;
                    end

                    fprintf('Block %d | selected scans = %d\n', cond_num, sum(block_idx));
                    append_log(log_file, sprintf('[BLOCK] %s %s Block-%02d selected_scans=%d', subj, run_name, cond_num, sum(block_idx)));

                    if sum(block_idx) < min_valid_timepoints_fc
                        warning('Too few scans for %s %s Block-%02d. Skipping.', subj, run_name, cond_num);
                        append_log(log_file, sprintf('[WARN] Too few scans for %s %s Block-%02d', subj, run_name, cond_num));
                        continue;
                    end

                    Yblk = Yclean(block_idx, :);
                    compute_and_save_fc(Yblk, roi_names, run_outdir, sprintf('Block-%02d', cond_num), ...
                        make_figs, min_valid_timepoints_fc, low_variance_threshold, clip_r_for_fisher, ...
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
fprintf('\n[FINISHED] All processing complete.\n');

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
function tf = same_space(file1, file2)
try
    V1 = spm_vol(file1);
    V2 = spm_vol(file2);

    if numel(V1) > 1
        V1 = V1(1);
    end
    if numel(V2) > 1
        V2 = V2(1);
    end

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
if isempty(ext)
    ext = '.nii';
end

out_roi = fullfile(outdir, [out_prefix ext]);

% Reuse an existing resliced mask when available.
if exist(out_roi, 'file')
    try
        if same_space(out_roi, func_file)
            append_log(log_file, sprintf('[RESLICE-SKIP] Existing resliced file reused: %s', out_roi));
            return;
        else
            delete(out_roi);
            append_log(log_file, sprintf('[RESLICE-REMAKE] Existing file removed due to space mismatch: %s', out_roi));
        end
    catch
        delete(out_roi);
        append_log(log_file, sprintf('[RESLICE-REMAKE] Existing file removed after validation failure: %s', out_roi));
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

if exist(tmp_copy, 'file')
    delete(tmp_copy);
end

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
    vals = Y(idx);
    X(t,:) = vals(:)';
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

[~,score,~] = pca(X, 'NumComponents', ncomp);
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
function session_num = parse_run_number(run_name)

tok = regexp(run_name, 'run_(\d+)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse session number from run name: %s', run_name);
end
session_num = str2double(tok{1});

end


%% ========================================================================
function Tsub = get_subject_condition_rows(Tcond, subj_label, session_num, has_real_subject)

subnum = sscanf(subj_label, 'sub-%d');
cand_rows = false(height(Tcond),1);

if ismember('subject', Tcond.Properties.VariableNames)
    subjcol = Tcond.subject;

    if isnumeric(subjcol)
        cand_rows = cand_rows | (subjcol == subnum);
    elseif iscell(subjcol)
        subjtxt = cell(size(subjcol));
        for i = 1:numel(subjcol)
            subjtxt{i} = char(string(subjcol{i}));
        end
        cand_rows = cand_rows | strcmp(subjtxt, subj_label) | ...
                               strcmp(subjtxt, sprintf('%d',subnum)) | ...
                               strcmp(subjtxt, sprintf('sub-%02d',subnum));
    elseif isstring(subjcol)
        subjtxt = cellstr(subjcol);
        cand_rows = cand_rows | strcmp(subjtxt, subj_label) | ...
                               strcmp(subjtxt, sprintf('%d',subnum)) | ...
                               strcmp(subjtxt, sprintf('sub-%02d',subnum));
    elseif ischar(subjcol)
        cand_rows = strcmp(cellstr(subjcol), subj_label);
    end
end

if has_real_subject
    rscol = Tcond.('real subject');
    if isnumeric(rscol)
        cand_rows = cand_rows | (rscol == subnum);
    elseif iscell(rscol)
        rstxt = cell(size(rscol));
        for i = 1:numel(rscol)
            rstxt{i} = char(string(rscol{i}));
        end
        cand_rows = cand_rows | strcmp(rstxt, subj_label) | ...
                               strcmp(rstxt, sprintf('%d',subnum)) | ...
                               strcmp(rstxt, sprintf('sub-%02d',subnum));
    elseif isstring(rscol)
        rstxt = cellstr(rscol);
        cand_rows = cand_rows | strcmp(rstxt, subj_label) | ...
                               strcmp(rstxt, sprintf('%d',subnum)) | ...
                               strcmp(rstxt, sprintf('sub-%02d',subnum));
    elseif ischar(rscol)
        cand_rows = cand_rows | strcmp(cellstr(rscol), subj_label);
    end
end

sess_ok = false(height(Tcond),1);
if isnumeric(Tcond.session)
    sess_ok = Tcond.session == session_num;
else
    sess_vals = Tcond.session;
    if isstring(sess_vals)
        sess_vals = cellstr(sess_vals);
    elseif ischar(sess_vals)
        sess_vals = cellstr(sess_vals);
    end

    if iscell(sess_vals)
        sess_ok = strcmp(sess_vals, sprintf('%d',session_num)) | ...
                  strcmp(sess_vals, sprintf('%.0f',session_num));
    end
end

Tsub = Tcond(cand_rows & sess_ok, :);

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

writecell(r_cell, fullfile(outdir, ['FC_' tag '_r.csv']));
writecell(z_cell, fullfile(outdir, ['FC_' tag '_z.csv']));

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
writetable(edge_tbl, fullfile(outdir, ['FC_' tag '_zvec.csv']));

save(fullfile(outdir, ['FC_' tag '.mat']), 'R','Z','roi_names','tag','subj','run_name');

if make_figs
    try
        f1 = figure('Visible','off','Color','w');
        imagesc(R);
        axis square;
        colorbar;
        caxis([-1 1]);
        title(['FC r - ' tag], 'Interpreter', 'none');
        set(gca,'XTick',1:nrois,'XTickLabel',roi_names,'XTickLabelRotation',45);
        set(gca,'YTick',1:nrois,'YTickLabel',roi_names);
        saveas(f1, fullfile(outdir, ['FC_' tag '_r_heatmap.png']));
        close(f1);

        f2 = figure('Visible','off','Color','w');
        imagesc(Z);
        axis square;
        colorbar;
        title(['FC z - ' tag], 'Interpreter', 'none');
        set(gca,'XTick',1:nrois,'XTickLabel',roi_names,'XTickLabelRotation',45);
        set(gca,'YTick',1:nrois,'YTickLabel',roi_names);
        saveas(f2, fullfile(outdir, ['FC_' tag '_z_heatmap.png']));
        close(f2);

        validz = z_vec(isfinite(z_vec));
        f3 = figure('Visible','off','Color','w');
        histogram(validz, 20);
        title(['FC z histogram - ' tag], 'Interpreter', 'none');
        xlabel('Fisher z');
        ylabel('Count');
        saveas(f3, fullfile(outdir, ['FC_' tag '_z_hist.png']));
        close(f3);
    catch figME
        warning('Figure saving failed for %s: %s', tag, figME.message);
        append_log(log_file, sprintf('[WARN] Figure saving failed for %s: %s', tag, figME.message));
    end
end

append_log(log_file, sprintf('[SAVED] %s | %s | %s', subj, run_name, tag));
fprintf('[SAVED] %s\n', fullfile(outdir, ['FC_' tag '.mat']));

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


%% ========================================================================
function export_group_level_fc(output_root, roi_names, log_file)

group_dir = fullfile(output_root, 'FC_group_exports');
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
edge_summary_rows = {};
edge_summary_header_done = false;

for rr = 1:numel(run_labels)
    run_name = run_labels{rr};
    run_num = parse_run_number(run_name);

    for tt = 1:numel(tag_labels)
        tag = tag_labels{tt};

        mats = {};
        subj_list = {};
        edge_table_this = [];

        for s = 1:numel(subdirs)
            subj = subdirs(s).name;
            matfile = fullfile(output_root, subj, run_name, ['FC_' tag '.mat']);

            if exist(matfile, 'file')
                S = load(matfile);

                if isfield(S, 'Z')
                    mats{end+1} = S.Z; %#ok<AGROW>
                    subj_list{end+1} = subj; %#ok<AGROW>

                    zvec_file = fullfile(output_root, subj, run_name, ['FC_' tag '_zvec.csv']);
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

        if nSub == 0
            continue;
        end

        nroi = size(mats{1},1);
        Zstack = nan(nroi, nroi, nSub);
        for k = 1:nSub
            Zstack(:,:,k) = mats{k};
        end

        Zmean = nanmean(Zstack, 3);
        Zmedian = nanmedian(Zstack, 3);
        Zmode = matrix_mode_3d(Zstack);

        write_matrix_csv_with_labels(fullfile(group_dir, sprintf('FC_group_R%d_%s_mean.csv', run_num, tag)), Zmean, roi_names);
        write_matrix_csv_with_labels(fullfile(group_dir, sprintf('FC_group_R%d_%s_median.csv', run_num, tag)), Zmedian, roi_names);
        write_matrix_csv_with_labels(fullfile(group_dir, sprintf('FC_group_R%d_%s_mode.csv', run_num, tag)), Zmode, roi_names);

        save(fullfile(group_dir, sprintf('FC_group_R%d_%s_stack.mat', run_num, tag)), ...
            'Zstack','Zmean','Zmedian','Zmode','roi_names','subj_list','run_name','tag');

        if ~isempty(edge_table_this)
            edge_names = unique(edge_table_this.edge, 'stable');
            summary_edge = table();

            edge_col = cell(numel(edge_names),1);
            mean_col = nan(numel(edge_names),1);
            median_col = nan(numel(edge_names),1);
            mode_col = nan(numel(edge_names),1);
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
                end
            end

            summary_edge.edge = edge_col;
            summary_edge.run = run_col;
            summary_edge.tag = tag_col;
            summary_edge.n = nsub_col;
            summary_edge.mean_z = mean_col;
            summary_edge.median_z = median_col;
            summary_edge.mode_z = mode_col;

            writetable(summary_edge, fullfile(group_dir, sprintf('FC_group_R%d_%s_edge_summary.csv', run_num, tag)));

            if ~edge_summary_header_done
                edge_summary_rows = summary_edge;
                edge_summary_header_done = true;
            else
                edge_summary_rows = [edge_summary_rows; summary_edge]; %#ok<AGROW>
            end
        end

        % allsubjects wide
        [wideTbl, headerCells] = make_allsubjects_wide_table(subj_list, mats, roi_names, run_name, tag);

        writetable(wideTbl, fullfile(group_dir, sprintf('FC_group_R%d_%s_allsubjects_wide.csv', run_num, tag)));

        if ~allsubject_header_done
            allsubject_rows = headerCells;
            allsubject_header_done = true;
        else
            allsubject_rows = [allsubject_rows; headerCells(2:end,:)]; %#ok<AGROW>
        end
    end
end

if edge_summary_header_done
    writetable(edge_summary_rows, fullfile(group_dir, 'FC_group_edge_summary.csv'));
end

if allsubject_header_done
    writecell(allsubject_rows, fullfile(group_dir, 'FC_allsubjects_wide.csv'));
end

inventory = build_fc_inventory(output_root);
writetable(inventory, fullfile(group_dir, 'FC_file_inventory.csv'));

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
            f = fullfile(output_root, subj, runs{r}, ['FC_' tags{t} '.mat']);

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
