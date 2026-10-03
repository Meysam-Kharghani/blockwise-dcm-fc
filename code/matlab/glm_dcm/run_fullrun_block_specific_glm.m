function run_fullrun_block_specific_glm
%% ================================================
% Full-run block-specific and run-level GLM construction
% - Per RUN: 3 GLMs (one per merged block)
% - Per RUN: 1 GLM pooling all events (ALL)
% - Uses CSV block definitions directly
% ================================================
clear; clc;
spm('defaults','FMRI');
spm_jobman('initcfg');

% ---------- Paths ----------
cfg = project_paths();
root_data = cfg.data_root;
root_out = cfg.glm_root;
csv_path = cfg.block_timing_csv;

% ---------- GLM defaults ----------
TR_default  = 3.0;
hpf         = 128;
fmri_t      = 16;
fmri_t0     = 8;

% ---------- Load block CSV ----------
B = readtable(csv_path, 'VariableNamingRule','preserve');

% Required columns in the block-timing table.
required_cols = {'subject','session','condition','onset','duration','real subject'};
for i = 1:numel(required_cols)
    if ~ismember(required_cols{i}, B.Properties.VariableNames)
        error('Column "%s" not found in CSV file.', required_cols{i});
    end
end

% String representation of the participant identifier.
real_sub_col = string(B.('real subject'));

% ---------- Loop subjects ----------
subs = dir(fullfile(root_data,'sub-*'));
for s = 1:numel(subs)
    if ~subs(s).isdir, continue; end

    sub_id  = subs(s).name;   % e.g., sub-01
    funcdir = fullfile(root_data, sub_id, 'func');
    if ~isfolder(funcdir)
        warning('%s: func not found. Skipping.', sub_id);
        continue;
    end

    fprintf('\n==============================\nProcessing %s\n==============================\n', sub_id);

    % Select block rows for the current participant.
    sub_mask = strcmp(real_sub_col, sub_id);
    Bsub = B(sub_mask, :);

    if isempty(Bsub)
        warning('%s: no block information found in CSV. Skipping.', sub_id);
        continue;
    end

    runs = dir(fullfile(funcdir,'run_*'));
    for r = 1:numel(runs)
        run_name = runs(r).name;   % e.g., run_1
        run_dir  = fullfile(funcdir, run_name);

        fprintf(' -> Parsing %s | %s ...\n', sub_id, run_name);

        % Parse the run number from the directory name.
        run_num = extract_run_number(run_name);
        if isnan(run_num)
            warning('%s | %s: could not parse run number. Skipping.', sub_id, run_name);
            continue;
        end

        % --- Collect scans ---
        scans = collect_scans_for_run(run_dir);
        if isempty(scans)
            warning('%s | %s: No BOLD scans found.', sub_id, run_name);
            continue;
        end

        % --- TR detection ---
        TR = detect_TR_from_scans(scans, TR_default);

        % --- Load events.tsv for ALL model ---
        evt = dir(fullfile(run_dir,'*events.tsv'));
        if isempty(evt)
            warning('%s | %s: events.tsv not found.', sub_id, run_name);
            continue;
        end
        evt_path = fullfile(evt(1).folder, evt(1).name);
        T = read_events_tsv(evt_path);

        [on, du, tt] = normalize_events_columns(T);
        if isempty(on)
            warning('%s | %s: could not read onset/duration/trial_type.', sub_id, run_name);
            continue;
        end
        du(isnan(du)) = 0.35;

        % --- Get 3 merged blocks from CSV for this subject/run ---
        run_mask = Bsub.session == run_num;
        Brun = Bsub(run_mask, :);

        if isempty(Brun)
            warning('%s | %s: no block rows found in CSV for this run.', sub_id, run_name);
            continue;
        end

        % Sort block rows by condition index.
        Brun = sortrows(Brun, 'condition');

        % Verify that the three expected blocks are present.
        conds = Brun.condition(:)';
        if numel(conds) ~= 3 || ~isequal(double(conds), [1 2 3])
            warning('%s | %s: expected 3 blocks with conditions [1 2 3], found: %s', ...
                sub_id, run_name, mat2str(double(conds)));
        end

        blocks = struct('condition',{},'block_onset',{},'block_duration',{});
        for k = 1:height(Brun)
            blocks(k).condition      = Brun.condition(k);
            blocks(k).block_onset    = Brun.onset(k);
            blocks(k).block_duration = Brun.duration(k);
        end

        % Save a summary TSV for QA
        out_tab_dir = fullfile(root_out, sub_id, 'GLM_blocks_summaries');
        ensure_dir(out_tab_dir);
        writetable(Brun, fullfile(out_tab_dir, sprintf('%s_blocks_3merged.tsv', run_name)), ...
            'FileType','text','Delimiter','\t');

        % --- Motion regressors (optional) ---
        Rfile = find_motion_regressors(run_dir);

        % ===== 1) Per-block GLMs: 3 per run =====
        for b = 1:numel(blocks)
            blk = blocks(b);

            glm_dir = fullfile(root_out, sub_id, sprintf('GLM_%s_block-%02d', run_name, b));
            reset_dir(glm_dir);

            mb = [];
            mb{1}.spm.stats.fmri_spec.dir = {glm_dir};
            mb{1}.spm.stats.fmri_spec.timing.units = 'secs';
            mb{1}.spm.stats.fmri_spec.timing.RT    = TR;
            mb{1}.spm.stats.fmri_spec.timing.fmri_t  = fmri_t;
            mb{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;

            mb{1}.spm.stats.fmri_spec.sess.scans = scans;

            mb{1}.spm.stats.fmri_spec.sess.cond(1).name     = sprintf('BLOCK_%02d', b);
            mb{1}.spm.stats.fmri_spec.sess.cond(1).onset    = blk.block_onset(:);
            mb{1}.spm.stats.fmri_spec.sess.cond(1).duration = blk.block_duration(:);
            mb{1}.spm.stats.fmri_spec.sess.cond(1).tmod     = 0;
            mb{1}.spm.stats.fmri_spec.sess.cond(1).pmod     = struct('name',{},'param',{},'poly',{});
            mb{1}.spm.stats.fmri_spec.sess.cond(1).orth     = 0;

            if ~isempty(Rfile)
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

            mb{2}.spm.stats.fmri_est.spmmat = {fullfile(glm_dir,'SPM.mat')};

            fprintf('    Building GLM (block %02d) ...\n', b);
            spm_jobman('run', mb);
        end

        % ===== 2) Run-wise GLM (ALL) =====
        glm_dir = fullfile(root_out, sub_id, sprintf('GLM_%s_ALL', run_name));
        reset_dir(glm_dir);

        mb = [];
        mb{1}.spm.stats.fmri_spec.dir = {glm_dir};
        mb{1}.spm.stats.fmri_spec.timing.units = 'secs';
        mb{1}.spm.stats.fmri_spec.timing.RT    = TR;
        mb{1}.spm.stats.fmri_spec.timing.fmri_t  = fmri_t;
        mb{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;

        mb{1}.spm.stats.fmri_spec.sess.scans = scans;

        mb{1}.spm.stats.fmri_spec.sess.cond(1).name     = 'ALL';
        mb{1}.spm.stats.fmri_spec.sess.cond(1).onset    = on(:);
        mb{1}.spm.stats.fmri_spec.sess.cond(1).duration = du(:);
        mb{1}.spm.stats.fmri_spec.sess.cond(1).tmod     = 0;
        mb{1}.spm.stats.fmri_spec.sess.cond(1).pmod     = struct('name',{},'param',{},'poly',{});
        mb{1}.spm.stats.fmri_spec.sess.cond(1).orth     = 0;

        if ~isempty(Rfile)
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

        mb{2}.spm.stats.fmri_est.spmmat = {fullfile(glm_dir,'SPM.mat')};

        fprintf('    Building GLM (ALL) ...\n');
        spm_jobman('run', mb);

        fprintf('    Completed %s | %s\n', sub_id, run_name);
    end
end

fprintf('\nCompleted all available participants.\n');
end


%% ===================== helpers =====================

function ensure_dir(p)
    if ~exist(p,'dir')
        mkdir(p);
    end
end

function reset_dir(p)
    if exist(p,'dir')
        rmdir(p,'s');
    end
    mkdir(p);
end

function scans = collect_scans_for_run(run_dir)
    scans = {};
    patt = {'sw*bold*.nii','dswausub*bold*.nii','*bold.nii','*bold*.nii'};
    for i = 1:numel(patt)
        L = dir(fullfile(run_dir, patt{i}));
        if ~isempty(L)
            if numel(L) == 1
                % 4D nifti
                nii = fullfile(L(1).folder, L(1).name);
                V = spm_vol(nii);
                nvol = numel(V);
                scans = arrayfun(@(k) sprintf('%s,%d', nii, k), 1:nvol, 'UniformOutput', false)';
            else
                % many 3D
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
        f = regexprep(scans{1},',\d+$','');
        V = spm_vol(f);
        if isfield(V,'private') && isfield(V.private,'timing') && ...
                isfield(V.private.timing,'tspace') && V.private.timing.tspace > 0
            TR = V.private.timing.tspace;
        elseif isfield(V,'pixdim') && numel(V.pixdim) >= 4 && V.pixdim(4) > 0
            TR = V.pixdim(4);
        end
    catch
    end
end

function T = read_events_tsv(p)
    fid = fopen(p,'r','n','UTF-8');
    raw = fread(fid,'*char')';
    fclose(fid);

    if ~isempty(raw) && strlength(raw) >= 3 && all(double(raw(1:3)) == [239 187 191])
        raw = raw(4:end);
    end

    lines = regexp(raw,'\r?\n','split');
    lines = lines(~cellfun(@isempty, lines));
    hdr   = regexp(lines{1},'\t','split');

    body = '';
    if numel(lines) > 1
        body = strjoin(lines(2:end), '\n');
    end

    data = textscan(body, repmat('%s',1,numel(hdr)), 'Delimiter','\t', 'CollectOutput', true);
    T = cell2table(data{1}, 'VariableNames', hdr);
end

function [on,du,tt] = normalize_events_columns(T)
    raw  = T.Properties.VariableNames;
    norm = lower(regexprep(raw,'[^a-z0-9]',''));
    f    = @(alts) find(ismember(norm, alts),1);

    iOn  = f({'onset','onsets','onsettime','onsetsec'});
    iDur = f({'duration','dur','durationsec','dursec'});
    iTT  = f({'trialtype','trial_type','condition','stimtype','category','type','stimulus'});

    on = [];
    du = [];
    tt = string([]);

    if isempty(iOn) || isempty(iTT)
        return;
    end

    colOn = T.(raw{iOn});
    if iscell(colOn) || isstring(colOn)
        on = str2double(string(colOn));
    else
        on = double(colOn);
    end

    if ~isempty(iDur)
        colDu = T.(raw{iDur});
        if iscell(colDu) || isstring(colDu)
            du = str2double(string(colDu));
        else
            du = double(colDu);
        end
    else
        du = 0.35 * ones(size(on));
    end
    du(isnan(du)) = 0.35;

    colTT = T.(raw{iTT});
    tt = string(lower(strtrim(string(colTT))));

    valid = ~isnan(on) & ~isnan(du) & tt ~= "";
    on = on(valid);
    du = du(valid);
    tt = tt(valid);
end

function Rfile = find_motion_regressors(run_dir)
    Rfile = '';
    patt = {'rp_*.txt','rp_*','art_regression_outliers*.mat'};
    for i = 1:numel(patt)
        L = dir(fullfile(run_dir,patt{i}));
        if ~isempty(L)
            Rfile = fullfile(L(1).folder, L(1).name);
            return;
        end
    end
end

function run_num = extract_run_number(run_name)
    tok = regexp(run_name, 'run_(\d+)', 'tokens', 'once');
    if isempty(tok)
        run_num = NaN;
    else
        run_num = str2double(tok{1});
    end
end
