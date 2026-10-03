function run_cropped_segmentwise_glm_3blocks(varargin)
% =========================================================================
% run_cropped_segmentwise_glm_3blocks
% =========================================================================
% Purpose
%   Build GLMs for block-specific DCM in a way that each block GLM sees ONLY
%   the scans belonging to that temporal segment.
%
% Difference from the primary full-run block-specific GLM:
%   Primary full-run specification:
%       sess.scans = all scans of the run
%       block identity was only represented by a block regressor
%
%   This version:
%       Block-01 GLM uses only scans from Block-01
%       Block-02 GLM uses only scans from Block-02
%       Block-03 GLM uses only scans from Block-03
%
% Output folder structure
%   <work_root>/glm_cropped/sub-01/GLM_run_1_block-01/SPM.mat
%   <work_root>/glm_cropped/sub-01/GLM_run_1_block-02/SPM.mat
%   <work_root>/glm_cropped/sub-01/GLM_run_1_block-03/SPM.mat
%   <work_root>/glm_cropped/sub-01/GLM_run_1_ALL/SPM.mat
%
% The generated GLMs are used by build_dcm_spm12_cropped_segments.m.
%
% -------------------------------------------------------------------------
% Optional name-value inputs:
%   'root_data'       : BIDS/preprocessed data root
%   'root_out'        : output GLM root
%   'csv_path'        : FC_blocks_COND.csv path
%   'TR_default'      : fallback TR if not found in NIfTI header. Default = 3.0
%   'hpf'             : high pass filter. Default = 128
%   'input_mode'      : 'events' or 'block'. Default = 'events'
%                       'events': use task events within each cropped segment
%                       'block' : one sustained input covering the segment
%   'overwrite'       : true/false. Default = true
% =========================================================================

% ------------------------- Defaults ---------------------------------------
cfg = project_paths();
P = inputParser;
P.addParameter('root_data', cfg.data_root, @ischar);
P.addParameter('root_out', cfg.cropped_glm_root, @ischar);
P.addParameter('csv_path', cfg.block_timing_csv, @ischar);
P.addParameter('TR_default', 3.0, @isnumeric);     % acquisition TR = 3 s; header overrides if present
P.addParameter('hpf',        128, @isnumeric);
P.addParameter('fmri_t',     16, @isnumeric);
P.addParameter('fmri_t0',    8, @isnumeric);
P.addParameter('input_mode', 'events', @(x) any(strcmpi(x, {'events','block'})));
P.addParameter('overwrite',  true, @islogical);
P.parse(varargin{:});
S = P.Results;

clc;
spm('defaults','FMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

fprintf('\n============================================================\n');
fprintf('CROPPED SEGMENT-WISE GLM BUILDER\n');
fprintf('root_data  : %s\n', S.root_data);
fprintf('root_out   : %s\n', S.root_out);
fprintf('csv_path   : %s\n', S.csv_path);
fprintf('input_mode : %s\n', S.input_mode);
fprintf('TR fallback: %.3f sec\n', S.TR_default);
fprintf('============================================================\n');

if ~exist(S.root_data, 'dir'), error('root_data not found: %s', S.root_data); end
if ~exist(S.csv_path,  'file'), error('csv_path not found: %s', S.csv_path); end
ensure_dir(S.root_out);

% ------------------------- Load block CSV --------------------------------
B = readtable(S.csv_path, 'VariableNamingRule','preserve');
required_cols = {'subject','session','condition','onset','duration'};
for i = 1:numel(required_cols)
    if ~ismember(required_cols{i}, B.Properties.VariableNames)
        error('Column "%s" not found in CSV file.', required_cols{i});
    end
end

% Use the explicit participant identifier column when present.
if ismember('real subject', B.Properties.VariableNames)
    subject_col = string(B.('real subject'));
else
    subject_col = string(B.subject);
end

% ------------------------- Subject loop ----------------------------------
subs = dir(fullfile(S.root_data, 'sub-*'));
subs = subs([subs.isdir]);

n_done = 0; n_fail = 0; n_skip = 0;
log_rows = {};

for si = 1:numel(subs)
    sub_id  = subs(si).name;
    funcdir = fullfile(S.root_data, sub_id, 'func');
    if ~isfolder(funcdir)
        warning('%s: func folder not found. Skipping.', sub_id);
        continue;
    end

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Subject: %s\n', sub_id);
    fprintf('------------------------------------------------------------\n');

    Bsub = B(strcmp(subject_col, sub_id), :);
    if isempty(Bsub)
        warning('%s: no block information in CSV. Skipping.', sub_id);
        n_skip = n_skip + 1;
        continue;
    end

    runs = dir(fullfile(funcdir, 'run_*'));
    runs = runs([runs.isdir]);
    [~, ix] = sort({runs.name}); runs = runs(ix);

    for ri = 1:numel(runs)
        run_name = runs(ri).name;      % e.g., run_1
        run_dir  = fullfile(funcdir, run_name);
        run_num  = extract_run_number(run_name);
        if isnan(run_num)
            warning('%s | %s: cannot parse run number. Skipping.', sub_id, run_name);
            continue;
        end

        fprintf('\n  >>> %s | %s\n', sub_id, run_name);

        scans_all = collect_scans_for_run(run_dir);
        if isempty(scans_all)
            warning('%s | %s: no BOLD scans found.', sub_id, run_name);
            n_skip = n_skip + 1; continue;
        end
        nvol_all = numel(scans_all);
        TR = detect_TR_from_scans(scans_all, S.TR_default);
        fprintf('      nvol = %d | TR = %.4f sec\n', nvol_all, TR);

        evt = dir(fullfile(run_dir, '*events.tsv'));
        if isempty(evt)
            warning('%s | %s: events.tsv not found. Skipping.', sub_id, run_name);
            n_skip = n_skip + 1; continue;
        end
        evt_path = fullfile(evt(1).folder, evt(1).name);
        T = read_events_tsv(evt_path);
        [on_all, du_all, tt_all] = normalize_events_columns(T); %#ok<ASGLU>
        if isempty(on_all)
            warning('%s | %s: event onsets not readable. Skipping.', sub_id, run_name);
            n_skip = n_skip + 1; continue;
        end

        % Rows for this subject/run from FC_blocks_COND.csv
        run_mask = to_double(Bsub.session) == run_num;
        Brun = Bsub(run_mask, :);
        if isempty(Brun)
            warning('%s | %s: no block rows for this run. Skipping.', sub_id, run_name);
            n_skip = n_skip + 1; continue;
        end
        Brun = sortrows(Brun, 'condition');

        % Motion regressors: crop rows to selected volume indices for each block.
        Rfile_all = find_motion_regressors(run_dir);

        % Save QA summary
        qa_dir = fullfile(S.root_out, sub_id, 'GLM_blocks_summaries');
        ensure_dir(qa_dir);
        writetable(Brun, fullfile(qa_dir, sprintf('%s_blocks_3merged_cropped.tsv', run_name)), ...
            'FileType','text','Delimiter','\t');

        % ===================== Cropped Block GLMs =========================
        for bi = 1:height(Brun)
            block_idx = to_double(Brun.condition(bi));
            if isnan(block_idx), block_idx = bi; end
            block_onset = to_double(Brun.onset(bi));
            block_dur   = to_double(Brun.duration(bi));
            block_end   = block_onset + block_dur;

            [idx_vol, segment_scan_start, segment_scan_end] = select_scans_for_segment(nvol_all, TR, block_onset, block_dur);
            if isempty(idx_vol) || numel(idx_vol) < 8
                warning('%s | %s | block-%02d: too few scans selected (%d). Skipping.', ...
                    sub_id, run_name, block_idx, numel(idx_vol));
                n_skip = n_skip + 1; continue;
            end
            scans_seg = scans_all(idx_vol);
            seg_duration = numel(scans_seg) * TR;

            glm_dir = fullfile(S.root_out, sub_id, sprintf('GLM_%s_block-%02d', run_name, block_idx));
            if S.overwrite
                reset_dir(glm_dir);
            else
                ensure_dir(glm_dir);
                if exist(fullfile(glm_dir,'SPM.mat'), 'file')
                    fprintf('      [SKIP] GLM exists: %s\n', glm_dir);
                    n_skip = n_skip + 1; continue;
                end
            end

            % Build regressors for cropped time series
            switch lower(S.input_mode)
                case 'events'
                    [seg_on, seg_du] = crop_events_to_segment(on_all, du_all, segment_scan_start, segment_scan_end);
                    if isempty(seg_on)
                        % Fallback: if no event found, use one sustained block regressor
                        seg_on = 0;
                        seg_du = seg_duration;
                        cond_name = sprintf('BLOCK_%02d_SEGMENT', block_idx);
                    else
                        cond_name = sprintf('TASK_EVENTS_BLOCK_%02d', block_idx);
                    end
                case 'block'
                    seg_on = 0;
                    seg_du = seg_duration;
                    cond_name = sprintf('BLOCK_%02d_SEGMENT', block_idx);
            end

            % Crop motion/confound rows to selected volumes.
            Rfile_seg = crop_motion_regressors(Rfile_all, idx_vol, glm_dir, sprintf('R_%s_block-%02d_cropped.txt', run_name, block_idx));

            try
                build_single_glm(glm_dir, scans_seg, TR, S.fmri_t, S.fmri_t0, S.hpf, ...
                    cond_name, seg_on, seg_du, Rfile_seg);
                n_done = n_done + 1;
                fprintf('      [OK] Block-%02d GLM | scans=%d | time %.2f-%.2f sec | input=%s\n', ...
                    block_idx, numel(scans_seg), segment_scan_start, segment_scan_end, cond_name);
                log_rows(end+1,:) = {sub_id, run_name, sprintf('Block-%02d',block_idx), numel(scans_seg), ...
                    TR, block_onset, block_end, segment_scan_start, segment_scan_end, cond_name}; %#ok<AGROW>
            catch ME
                n_fail = n_fail + 1;
                warning('%s | %s | block-%02d GLM failed: %s', sub_id, run_name, block_idx, ME.message);
            end
        end

        % ===================== ALL GLM, full run ==========================
        glm_dir = fullfile(S.root_out, sub_id, sprintf('GLM_%s_ALL', run_name));
        if S.overwrite
            reset_dir(glm_dir);
        else
            ensure_dir(glm_dir);
        end
        try
            Rfile_all_for_spm = prepare_full_motion_regressors(Rfile_all, nvol_all, glm_dir, sprintf('R_%s_ALL.txt', run_name));
            build_single_glm(glm_dir, scans_all, TR, S.fmri_t, S.fmri_t0, S.hpf, ...
                'ALL_TASK_EVENTS', on_all(:), du_all(:), Rfile_all_for_spm);
            n_done = n_done + 1;
            fprintf('      [OK] ALL GLM | scans=%d\n', nvol_all);
            log_rows(end+1,:) = {sub_id, run_name, 'ALL', nvol_all, TR, 0, nvol_all*TR, 0, nvol_all*TR, 'ALL_TASK_EVENTS'}; %#ok<AGROW>
        catch ME
            n_fail = n_fail + 1;
            warning('%s | %s | ALL GLM failed: %s', sub_id, run_name, ME.message);
        end
    end
end

% ------------------------- Write QA log ----------------------------------
if ~isempty(log_rows)
    Q = cell2table(log_rows, 'VariableNames', ...
        {'subject','run','model','n_scans','TR','csv_block_onset','csv_block_end', ...
         'actual_first_scan_time','actual_segment_end_time','input_name'});
    writetable(Q, fullfile(S.root_out, 'cropped_segmentwise_glm_QA.csv'));
end

fprintf('\n============================================================\n');
fprintf('CROPPED GLM SUMMARY\n');
fprintf('Successfully built : %d\n', n_done);
fprintf('Skipped            : %d\n', n_skip);
fprintf('Failed             : %d\n', n_fail);
fprintf('QA file            : %s\n', fullfile(S.root_out, 'cropped_segmentwise_glm_QA.csv'));
fprintf('============================================================\n');
end

%% =========================================================================
% GLM builder helper
% =========================================================================
function build_single_glm(glm_dir, scans, TR, fmri_t, fmri_t0, hpf, cond_name, onsets, durations, Rfile)
    mb = {};
    mb{1}.spm.stats.fmri_spec.dir = {glm_dir};
    mb{1}.spm.stats.fmri_spec.timing.units = 'secs';
    mb{1}.spm.stats.fmri_spec.timing.RT    = TR;
    mb{1}.spm.stats.fmri_spec.timing.fmri_t  = fmri_t;
    mb{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;

    mb{1}.spm.stats.fmri_spec.sess.scans = scans(:);
    mb{1}.spm.stats.fmri_spec.sess.cond(1).name     = cond_name;
    mb{1}.spm.stats.fmri_spec.sess.cond(1).onset    = onsets(:);
    mb{1}.spm.stats.fmri_spec.sess.cond(1).duration = durations(:);
    mb{1}.spm.stats.fmri_spec.sess.cond(1).tmod     = 0;
    mb{1}.spm.stats.fmri_spec.sess.cond(1).pmod     = struct('name',{},'param',{},'poly',{});
    mb{1}.spm.stats.fmri_spec.sess.cond(1).orth     = 0;

    if ~isempty(Rfile) && exist(Rfile, 'file')
        mb{1}.spm.stats.fmri_spec.sess.multi_reg = {Rfile};
    else
        mb{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    end
    mb{1}.spm.stats.fmri_spec.sess.hpf = hpf;

    mb{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    mb{1}.spm.stats.fmri_spec.volt   = 1;
    mb{1}.spm.stats.fmri_spec.global = 'None';
    mb{1}.spm.stats.fmri_spec.mask   = {''};
    mb{1}.spm.stats.fmri_spec.cvi    = 'AR(1)';

    mb{2}.spm.stats.fmri_est.spmmat = {fullfile(glm_dir, 'SPM.mat')};
    spm_jobman('run', mb);
end

%% =========================================================================
% Segment and events helpers
% =========================================================================
function [idx, seg_start, seg_end] = select_scans_for_segment(nvol, TR, onset, duration)
    vol_times = (0:nvol-1)' * TR;
    block_end = onset + duration;
    % Use acquisition onset times. This keeps the cropped time series as a
    % proper consecutive SPM session after relabelling first selected scan as t=0.
    idx = find(vol_times >= onset & vol_times < block_end);
    if isempty(idx)
        % fallback to nearest range
        idx1 = max(1, floor(onset/TR) + 1);
        idx2 = min(nvol, ceil(block_end/TR));
        idx = (idx1:idx2)';
    end
    seg_start = vol_times(idx(1));
    seg_end   = vol_times(idx(end)) + TR;
end

function [seg_on, seg_du] = crop_events_to_segment(on, du, seg_start, seg_end)
    ev_start = on(:);
    ev_end   = on(:) + du(:);
    keep = ev_end > seg_start & ev_start < seg_end;
    ev_start = ev_start(keep);
    ev_end   = ev_end(keep);
    if isempty(ev_start)
        seg_on = [];
        seg_du = [];
        return;
    end
    clipped_start = max(ev_start, seg_start);
    clipped_end   = min(ev_end, seg_end);
    seg_on = clipped_start - seg_start;
    seg_du = clipped_end - clipped_start;
    seg_du(seg_du <= 0 | isnan(seg_du)) = 0.35;
end

%% =========================================================================
% Motion/confound helpers
% =========================================================================
function Rfile_seg = crop_motion_regressors(Rfile_all, idx_vol, glm_dir, out_name)
    Rfile_seg = '';
    if isempty(Rfile_all) || ~exist(Rfile_all, 'file')
        return;
    end
    [~,~,ext] = fileparts(Rfile_all);
    if ~strcmpi(ext, '.txt')
        warning('Motion file is not a .txt numeric matrix; skipping crop: %s', Rfile_all);
        return;
    end
    try
        R = readmatrix(Rfile_all);
        if isempty(R) || size(R,1) < max(idx_vol)
            warning('Motion file rows do not match scans; skipping: %s', Rfile_all);
            return;
        end
        Rseg = R(idx_vol, :);
        Rfile_seg = fullfile(glm_dir, out_name);
        writematrix(Rseg, Rfile_seg, 'Delimiter', 'tab');
    catch ME
        warning('Could not crop motion regressors: %s', ME.message);
        Rfile_seg = '';
    end
end

function Rfile_out = prepare_full_motion_regressors(Rfile_all, nvol, glm_dir, out_name)
    Rfile_out = '';
    if isempty(Rfile_all) || ~exist(Rfile_all, 'file')
        return;
    end
    [~,~,ext] = fileparts(Rfile_all);
    if ~strcmpi(ext, '.txt')
        Rfile_out = Rfile_all;
        return;
    end
    try
        R = readmatrix(Rfile_all);
        if size(R,1) == nvol
            Rfile_out = Rfile_all;
        elseif size(R,1) > nvol
            Rfile_out = fullfile(glm_dir, out_name);
            writematrix(R(1:nvol,:), Rfile_out, 'Delimiter', 'tab');
        else
            warning('Full motion rows (%d) less than nvol (%d); ignoring motion.', size(R,1), nvol);
            Rfile_out = '';
        end
    catch
        Rfile_out = Rfile_all;
    end
end

%% =========================================================================
% Existing pipeline helpers, made robust
% =========================================================================
function ensure_dir(p)
    if ~exist(p,'dir'), mkdir(p); end
end

function reset_dir(p)
    if exist(p,'dir'), rmdir(p,'s'); end
    mkdir(p);
end

function scans = collect_scans_for_run(run_dir)
    scans = {};
    patt = {'sw*bold*.nii','dswausub*bold*.nii','*bold.nii','*bold*.nii'};
    for i = 1:numel(patt)
        L = dir(fullfile(run_dir, patt{i}));
        if ~isempty(L)
            if numel(L) == 1
                nii = fullfile(L(1).folder, L(1).name);
                V = spm_vol(nii);
                nvol = numel(V);
                scans = arrayfun(@(k) sprintf('%s,%d', nii, k), 1:nvol, 'UniformOutput', false)';
            else
                names = sort({L.name});
                scans = strcat(fullfile(run_dir, names(:)), ',1');
            end
            return;
        end
    end
end

function TR = detect_TR_from_scans(scans, TR_fallback)
    TR = TR_fallback;
    try
        f = regexprep(scans{1}, ',\d+$', '');
        V = spm_vol(f);
        if isfield(V,'private') && isfield(V.private,'timing') && ...
                isfield(V.private.timing,'tspace') && V.private.timing.tspace > 0
            TR = V.private.timing.tspace;
        elseif isfield(V,'pixdim') && numel(V.pixdim) >= 4 && V.pixdim(4) > 0
            TR = V.pixdim(4);
        end
    catch
        TR = TR_fallback;
    end
end

function T = read_events_tsv(p)
    fid = fopen(p,'r','n','UTF-8');
    raw = fread(fid,'*char')';
    fclose(fid);
    if ~isempty(raw) && numel(raw) >= 3 && all(double(raw(1:3)) == [239 187 191])
        raw = raw(4:end);
    end
    lines = regexp(raw,'\r?\n','split');
    lines = lines(~cellfun(@isempty, lines));
    hdr = regexp(lines{1},'\t','split');
    body = '';
    if numel(lines) > 1, body = strjoin(lines(2:end), '\n'); end
    data = textscan(body, repmat('%s',1,numel(hdr)), 'Delimiter','\t', 'CollectOutput', true);
    T = cell2table(data{1}, 'VariableNames', hdr);
end

function [on,du,tt] = normalize_events_columns(T)
    raw  = T.Properties.VariableNames;
    norm = lower(regexprep(raw,'[^a-z0-9]',''));
    f = @(alts) find(ismember(norm, alts), 1);
    iOn  = f({'onset','onsets','onsettime','onsetsec'});
    iDur = f({'duration','dur','durationsec','dursec'});
    iTT  = f({'trialtype','trialtype','condition','stimtype','category','type','stimulus'});
    on = []; du = []; tt = string([]);
    if isempty(iOn), return; end
    colOn = T.(raw{iOn});
    on = str2double(string(colOn));
    if ~isempty(iDur)
        du = str2double(string(T.(raw{iDur})));
    else
        du = 0.35 * ones(size(on));
    end
    du(isnan(du)) = 0.35;
    if ~isempty(iTT)
        tt = string(lower(strtrim(string(T.(raw{iTT})))));
    else
        tt = repmat("event", size(on));
    end
    valid = ~isnan(on) & ~isnan(du);
    on = on(valid); du = du(valid); tt = tt(valid);
end

function Rfile = find_motion_regressors(run_dir)
    Rfile = '';
    patt = {'rp_*.txt','rp_*txt','art_regression_outliers*.txt'};
    for i = 1:numel(patt)
        L = dir(fullfile(run_dir, patt{i}));
        if ~isempty(L)
            Rfile = fullfile(L(1).folder, L(1).name);
            return;
        end
    end
end

function run_num = extract_run_number(run_name)
    tok = regexp(run_name, 'run_(\d+)', 'tokens', 'once');
    if isempty(tok), run_num = NaN; else, run_num = str2double(tok{1}); end
end

function x = to_double(v)
    if isnumeric(v)
        x = double(v);
    elseif iscell(v)
        x = str2double(string(v{1}));
    else
        x = str2double(string(v));
    end
end
