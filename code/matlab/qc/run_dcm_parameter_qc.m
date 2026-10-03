function run_dcm_parameter_qc(root_dir, out_dir, varargin)
%RUN_DCM_PARAMETER_QC  Quality control and posterior edge prevalence for SPM fMRI DCMs.
%
% This script recursively finds estimated DCM_*.mat files, extracts model-level
% diagnostics and A-matrix posterior summaries, and writes publication-ready QC
% tables/figures.
%
% Example:
%   cfg = project_paths();
%   run_dcm_parameter_qc(cfg.glm_root, fullfile(cfg.work_root, 'qc'));
%
%       'roi_names', {'V1_L','lOTC','vOTC','IFG_L'}, ...
%       'pp_threshold', 0.95, ...
%       'abs_fallback_threshold', 1e-6, ...
%       'make_figures', true);
%
% Main outputs:
%   dcm_qc_model_summary.csv
%   dcm_qc_A_edges_long.csv
%   dcm_qc_A_edge_prevalence_by_run_block.csv
%   dcm_qc_warnings.txt
%   fig01_dcm_qc_overview.png/pdf
%   fig02_A_edge_prevalence_heatmaps.png/pdf
%   fig03_A_edge_uncertainty_heatmaps.png/pdf
%
% Notes:
%   - A(i,j) is interpreted as SOURCE j -> TARGET i.
%   - "present_pp" uses posterior sign probability when Cp is available:
%       max(P(theta > 0), P(theta < 0)) >= pp_threshold
%   - If Cp is unavailable, the script falls back to:
%       abs(Ep.A) > abs_fallback_threshold
%   - Self-connections are included but should be interpreted as DCM self-inhibition
%     / local gain terms, not ordinary directed edges between different ROIs.

% ---------------------------- inputs -------------------------------------
p = inputParser;
p.addRequired('root_dir', @(x)ischar(x) || isstring(x));
p.addRequired('out_dir',  @(x)ischar(x) || isstring(x));
p.addParameter('roi_names', {'V1_L','lOTC','vOTC','IFG_L'}, @(x)iscell(x) || isstring(x));
p.addParameter('pp_threshold', 0.95, @(x)isnumeric(x) && isscalar(x));
p.addParameter('abs_fallback_threshold', 1e-6, @(x)isnumeric(x) && isscalar(x));
p.addParameter('include_all_tag', true, @(x)islogical(x) && isscalar(x));
p.addParameter('make_figures', true, @(x)islogical(x) && isscalar(x));
p.addParameter('file_pattern', 'DCM*.mat', @(x)ischar(x) || isstring(x));
p.parse(root_dir, out_dir, varargin{:});

root_dir = char(p.Results.root_dir);
out_dir  = char(p.Results.out_dir);
roi_names_default = cellstr(p.Results.roi_names);
pp_threshold = p.Results.pp_threshold;
abs_fallback_threshold = p.Results.abs_fallback_threshold;
include_all_tag = p.Results.include_all_tag;
make_figures = p.Results.make_figures;
file_pattern = char(p.Results.file_pattern);

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

fprintf('\n=== DCM QC: searching for %s under:\n%s\n', file_pattern, root_dir);
files = dir(fullfile(root_dir, '**', file_pattern));
files = files(~[files.isdir]);

% Avoid accidentally reading group-level PEB/BMA files if they match pattern.
keep = true(numel(files),1);
for i = 1:numel(files)
    low = lower(fullfile(files(i).folder, files(i).name));
    if contains(low, 'peb') || contains(low, 'bma') || contains(low, 'qc')
        keep(i) = false;
    end
end
files = files(keep);

if isempty(files)
    error('No DCM files found. Check root_dir and file_pattern.');
end
fprintf('Found %d candidate DCM files.\n', numel(files));

model_rows = {};
edge_rows  = {};
warnings   = {};
model_idx  = 0;

for f = 1:numel(files)
    fpath = fullfile(files(f).folder, files(f).name);
    meta = parse_dcm_path_metadata(fpath);
    if ~include_all_tag && strcmpi(meta.tag, 'ALL')
        continue;
    end

    model_idx = model_idx + 1;
    status = 'OK';
    err_msg = '';

    try
        S = load(fpath);
        if isfield(S, 'DCM')
            DCM = S.DCM;
        else
            fn = fieldnames(S);
            DCM = S.(fn{1});
            warnings{end+1,1} = sprintf('Loaded first variable instead of DCM in %s', fpath); %#ok<AGROW>
        end

        if ~isfield(DCM, 'Ep') || ~isfield(DCM.Ep, 'A')
            error('Missing DCM.Ep.A');
        end

        A = DCM.Ep.A;
        nROI = size(A,1);
        roi_names = get_roi_names(DCM, roi_names_default, nROI);

        [ve_percent, r2_value, n_time, n_regions_data] = compute_variance_explained(DCM);
        free_energy = get_numeric_field(DCM, 'F');
        n_inputs = get_n_inputs(DCM);
        max_abs_A = max(abs(A(:)));
        mean_abs_A = mean(abs(A(:)), 'omitnan');
        n_nan_A = sum(isnan(A(:)));
        n_inf_A = sum(isinf(A(:)));

        [A_sd, A_ci_low, A_ci_high, A_pp_dir, A_present_pp, A_present_ci, method_used] = ...
            posterior_A_summary(DCM, A, pp_threshold, abs_fallback_threshold);

        if strcmp(method_used, 'abs_value_fallback')
            warnings{end+1,1} = sprintf('Cp unavailable/unusable; used abs-value fallback for %s', fpath); %#ok<AGROW>
        end

        if any(~isfinite(A(:)))
            status = 'BAD_A_NONFINITE';
        end

        model_rows(end+1,:) = {model_idx, meta.subject, meta.run, meta.tag, fpath, status, err_msg, ...
            free_energy, ve_percent, r2_value, n_time, n_regions_data, n_inputs, nROI, ...
            max_abs_A, mean_abs_A, n_nan_A, n_inf_A, method_used}; %#ok<AGROW>

        for target_i = 1:nROI
            for source_j = 1:nROI
                source = roi_names{source_j};
                target = roi_names{target_i};
                is_self = source_j == target_i;
                edge_name = sprintf('%s -> %s', source, target);

                edge_rows(end+1,:) = {model_idx, meta.subject, meta.run, meta.tag, source, target, edge_name, is_self, ...
                    A(target_i, source_j), A_sd(target_i, source_j), A_ci_low(target_i, source_j), A_ci_high(target_i, source_j), ...
                    A_pp_dir(target_i, source_j), A_present_pp(target_i, source_j), A_present_ci(target_i, source_j), ...
                    free_energy, ve_percent, status, fpath, method_used}; %#ok<AGROW>
            end
        end

    catch ME
        status = 'FAILED_LOAD_OR_QC';
        err_msg = ME.message;
        warnings{end+1,1} = sprintf('FAILED: %s | %s', fpath, err_msg); %#ok<AGROW>
        model_rows(end+1,:) = {model_idx, meta.subject, meta.run, meta.tag, fpath, status, err_msg, ...
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'failed'}; %#ok<AGROW>
    end
end

% ---------------------------- tables -------------------------------------
model_varnames = {'model_id','subject','run','tag','file_path','status','error_message', ...
    'free_energy_F','variance_explained_percent','R2','n_timepoints','n_regions_in_data','n_inputs','n_roi', ...
    'max_abs_A','mean_abs_A','n_nan_A','n_inf_A','presence_method'};
model_tbl = cell2table(model_rows, 'VariableNames', model_varnames);

edge_varnames = {'model_id','subject','run','tag','source','target','edge','is_self', ...
    'A_value','A_sd','A_ci_low','A_ci_high','A_pp_direction','present_pp','present_ci', ...
    'free_energy_F','variance_explained_percent','model_status','file_path','presence_method'};
if isempty(edge_rows)
    edge_tbl = cell2table(cell(0,numel(edge_varnames)), 'VariableNames', edge_varnames);
else
    edge_tbl = cell2table(edge_rows, 'VariableNames', edge_varnames);
end

model_csv = fullfile(out_dir, 'dcm_qc_model_summary.csv');
edge_csv  = fullfile(out_dir, 'dcm_qc_A_edges_long.csv');
writetable(model_tbl, model_csv);
writetable(edge_tbl, edge_csv);

summary_tbl = summarize_edges(edge_tbl);
summary_csv = fullfile(out_dir, 'dcm_qc_A_edge_prevalence_by_run_block.csv');
writetable(summary_tbl, summary_csv);

% Subject-level model availability matrix.
avail_tbl = summarize_model_availability(model_tbl);
writetable(avail_tbl, fullfile(out_dir, 'dcm_qc_model_availability.csv'));

% Warning report.
warn_file = fullfile(out_dir, 'dcm_qc_warnings.txt');
fid = fopen(warn_file, 'w');
fprintf(fid, 'DCM QC warning report\n');
fprintf(fid, 'Root: %s\n', root_dir);
fprintf(fid, 'Generated: %s\n\n', datestr(now));
if isempty(warnings)
    fprintf(fid, 'No warnings.\n');
else
    for i = 1:numel(warnings)
        fprintf(fid, '%03d. %s\n', i, warnings{i});
    end
end
fclose(fid);

fprintf('\nWrote:\n  %s\n  %s\n  %s\n  %s\n', model_csv, edge_csv, summary_csv, warn_file);

% ---------------------------- figures ------------------------------------
if make_figures
    try
        make_qc_figures(model_tbl, edge_tbl, summary_tbl, roi_names_default, out_dir);
    catch ME
        warning('Could not make figures: %s', ME.message);
    end
end

fprintf('\nDone.\n');
end

% ========================================================================
% local functions
% ========================================================================
function meta = parse_dcm_path_metadata(fpath)
    low = lower(fpath);
    [~, fname, ~] = fileparts(fpath);

    sub = regexp(fpath, 'sub[-_]?\d+', 'match', 'once', 'ignorecase');
    if isempty(sub)
        sub = regexp(fname, 's\d+', 'match', 'once', 'ignorecase');
    end
    if isempty(sub); sub = 'unknown_subject'; end

    run_tok = regexp(low, 'run[-_ ]?0?([12])', 'tokens', 'once');
    if isempty(run_tok)
        run = 'unknown_run';
    else
        run = sprintf('run_%s', run_tok{1});
    end

    if contains(low, 'block-01') || contains(low, 'block_01') || contains(low, 'block01')
        tag = 'Block-01';
    elseif contains(low, 'block-02') || contains(low, 'block_02') || contains(low, 'block02')
        tag = 'Block-02';
    elseif contains(low, 'block-03') || contains(low, 'block_03') || contains(low, 'block03')
        tag = 'Block-03';
    elseif contains(low, 'all')
        tag = 'ALL';
    else
        tag = 'unknown_tag';
    end

    meta = struct('subject', sub, 'run', run, 'tag', tag);
end

function roi_names = get_roi_names(DCM, roi_names_default, nROI)
    roi_names = {};
    try
        if isfield(DCM, 'Y') && isfield(DCM.Y, 'name') && numel(DCM.Y.name) == nROI
            roi_names = cellstr(DCM.Y.name);
        end
    catch
        roi_names = {};
    end
    if isempty(roi_names)
        if numel(roi_names_default) == nROI
            roi_names = roi_names_default;
        else
            roi_names = arrayfun(@(x)sprintf('ROI_%02d', x), 1:nROI, 'UniformOutput', false);
        end
    end
    roi_names = matlab.lang.makeValidName(roi_names, 'ReplacementStyle', 'delete');
end

function val = get_numeric_field(S, fieldname)
    if isfield(S, fieldname) && isnumeric(S.(fieldname)) && isscalar(S.(fieldname))
        val = S.(fieldname);
    else
        val = NaN;
    end
end

function n_inputs = get_n_inputs(DCM)
    n_inputs = NaN;
    try
        if isfield(DCM, 'U') && isfield(DCM.U, 'u')
            n_inputs = size(DCM.U.u, 2);
        elseif isfield(DCM, 'U') && isfield(DCM.U, 'name')
            n_inputs = numel(DCM.U.name);
        end
    catch
        n_inputs = NaN;
    end
end

function [ve_percent, r2_value, n_time, n_regions_data] = compute_variance_explained(DCM)
    ve_percent = NaN;
    r2_value = NaN;
    n_time = NaN;
    n_regions_data = NaN;

    y_obs = [];
    y_pred = [];

    try
        if isfield(DCM, 'Y') && isfield(DCM.Y, 'y')
            y_obs = DCM.Y.y;
        end
        if isfield(DCM, 'y')
            y_pred = DCM.y;
        elseif isfield(DCM, 'H')
            y_pred = DCM.H;
        end
    catch
        return;
    end

    if isempty(y_obs)
        return;
    end
    n_time = size(y_obs, 1);
    n_regions_data = size(y_obs, 2);

    if isempty(y_pred) || ~isequal(size(y_obs), size(y_pred))
        return;
    end

    yy = y_obs(:);
    yh = y_pred(:);
    ok = isfinite(yy) & isfinite(yh);
    yy = yy(ok);
    yh = yh(ok);

    if numel(yy) < 3 || var(yy) == 0
        return;
    end

    resid = yy - yh;
    sse = sum(resid.^2);
    sst = sum((yy - mean(yy)).^2);
    r2_value = 1 - sse / sst;
    ve_percent = 100 * r2_value;
end

function [A_sd, A_ci_low, A_ci_high, A_pp_dir, A_present_pp, A_present_ci, method_used] = ...
    posterior_A_summary(DCM, A, pp_threshold, abs_fallback_threshold)

    A_sd = NaN(size(A));
    A_ci_low = NaN(size(A));
    A_ci_high = NaN(size(A));
    A_pp_dir = NaN(size(A));
    A_present_pp = false(size(A));
    A_present_ci = false(size(A));
    method_used = 'abs_value_fallback';

    if isfield(DCM, 'Cp') && ~isempty(DCM.Cp) && isfield(DCM, 'Ep') && isfield(DCM.Ep, 'A')
        try
            idx_A = get_spm_field_indices_A(DCM.Ep);
            cp_diag = diag(DCM.Cp);
            if numel(idx_A) == numel(A) && max(idx_A) <= numel(cp_diag)
                A_var = reshape(cp_diag(idx_A), size(A));
                A_sd = sqrt(max(A_var, 0));
                A_ci_low = A - 1.96 .* A_sd;
                A_ci_high = A + 1.96 .* A_sd;

                z = A ./ A_sd;
                z(A_sd == 0 & A > 0) = Inf;
                z(A_sd == 0 & A < 0) = -Inf;
                z(A_sd == 0 & A == 0) = 0;

                % Normal CDF using erf to avoid requiring Statistics Toolbox.
                p_pos = 1 - normal_cdf(0, A, A_sd); % P(theta > 0)
                p_neg = normal_cdf(0, A, A_sd);     % P(theta < 0)
                A_pp_dir = max(p_pos, p_neg);
                A_pp_dir(~isfinite(A_pp_dir)) = NaN;

                A_present_pp = A_pp_dir >= pp_threshold;
                A_present_ci = (A_ci_low > 0 & A_ci_high > 0) | (A_ci_low < 0 & A_ci_high < 0);
                method_used = 'posterior_probability_from_Cp';
                return;
            end
        catch
            % fall through to abs-value fallback
        end
    end

    % Fallback: less rigorous. Use only as descriptive threshold.
    A_present_pp = abs(A) > abs_fallback_threshold;
    A_present_ci = A_present_pp;
end

function p = normal_cdf(x, mu, sigma)
    p = NaN(size(mu));
    ok = isfinite(mu) & isfinite(sigma) & sigma > 0;
    p(ok) = 0.5 .* (1 + erf((x - mu(ok)) ./ (sigma(ok) .* sqrt(2))));
    p(sigma == 0 & mu < x) = 1;
    p(sigma == 0 & mu > x) = 0;
    p(sigma == 0 & mu == x) = 0.5;
end

function idx_A = get_spm_field_indices_A(Ep)
    if exist('spm_fieldindices', 'file') == 2
        idx_A = spm_fieldindices(Ep, 'A');
    else
        % Minimal fallback matching MATLAB field order. Prefer SPM function if available.
        [~, paths] = local_vec_with_paths(Ep, '');
        idx_A = find(startsWith(paths, '.A('));
    end
end

function [v, paths] = local_vec_with_paths(x, prefix)
    v = [];
    paths = strings(0,1);
    if isstruct(x)
        fn = fieldnames(x);
        for k = 1:numel(fn)
            [vk, pk] = local_vec_with_paths(x.(fn{k}), [prefix '.' fn{k}]);
            v = [v; vk]; %#ok<AGROW>
            paths = [paths; pk]; %#ok<AGROW>
        end
    elseif isnumeric(x) || islogical(x)
        xv = x(:);
        v = xv;
        paths = strings(numel(xv),1);
        for ii = 1:numel(xv)
            paths(ii) = sprintf('%s(%d)', prefix, ii);
        end
    elseif iscell(x)
        for k = 1:numel(x)
            [vk, pk] = local_vec_with_paths(x{k}, sprintf('%s{%d}', prefix, k));
            v = [v; vk]; %#ok<AGROW>
            paths = [paths; pk]; %#ok<AGROW>
        end
    end
end

function summary_tbl = summarize_edges(edge_tbl)
    if isempty(edge_tbl) || height(edge_tbl) == 0
        summary_tbl = table();
        return;
    end

    ok = strcmp(edge_tbl.model_status, 'OK') | strcmp(edge_tbl.model_status, 'BAD_A_NONFINITE');
    E = edge_tbl(ok,:);
    key_tbl = E(:, {'run','tag','source','target','edge','is_self'});
    [keys, ~, g] = unique(key_tbl, 'rows');

    nG = height(keys);
    n_models = zeros(nG,1);
    n_present_pp = zeros(nG,1);
    n_present_ci = zeros(nG,1);
    prevalence_pp = nan(nG,1);
    prevalence_ci = nan(nG,1);
    mean_all = nan(nG,1);
    mean_present_pp = nan(nG,1);
    median_all = nan(nG,1);
    sd_all = nan(nG,1);
    sem_all = nan(nG,1);
    mean_sd_posterior = nan(nG,1);
    mean_pp_direction = nan(nG,1);
    n_subjects = zeros(nG,1);

    for k = 1:nG
        ix = g == k;
        vals = E.A_value(ix);
        pres_pp = logical(E.present_pp(ix));
        pres_ci = logical(E.present_ci(ix));
        pp = E.A_pp_direction(ix);
        psd = E.A_sd(ix);
        subs = E.subject(ix);

        n_models(k) = sum(ix);
        n_subjects(k) = numel(unique(subs));
        n_present_pp(k) = sum(pres_pp);
        n_present_ci(k) = sum(pres_ci);
        prevalence_pp(k) = n_present_pp(k) ./ n_models(k);
        prevalence_ci(k) = n_present_ci(k) ./ n_models(k);
        mean_all(k) = mean(vals, 'omitnan');
        median_all(k) = median(vals, 'omitnan');
        sd_all(k) = std(vals, 'omitnan');
        sem_all(k) = sd_all(k) ./ sqrt(sum(isfinite(vals)));
        if any(pres_pp)
            mean_present_pp(k) = mean(vals(pres_pp), 'omitnan');
        end
        mean_sd_posterior(k) = mean(psd, 'omitnan');
        mean_pp_direction(k) = mean(pp, 'omitnan');
    end

    summary_tbl = [keys, table(n_models, n_subjects, n_present_pp, prevalence_pp, n_present_ci, prevalence_ci, ...
        mean_all, mean_present_pp, median_all, sd_all, sem_all, mean_sd_posterior, mean_pp_direction)];
end

function avail_tbl = summarize_model_availability(model_tbl)
    if isempty(model_tbl) || height(model_tbl) == 0
        avail_tbl = table();
        return;
    end
    key_tbl = model_tbl(:, {'subject','run','tag'});
    [keys, ~, g] = unique(key_tbl, 'rows');
    n_files = accumarray(g, 1);
    n_ok = accumarray(g, strcmp(model_tbl.status, 'OK'));
    avail_tbl = [keys, table(n_files, n_ok)];
end

function make_qc_figures(model_tbl, edge_tbl, summary_tbl, roi_names_default, out_dir)
    fprintf('Making QC figures...\n');
    make_qc_overview_figure(model_tbl, out_dir);
    make_prevalence_heatmaps(summary_tbl, roi_names_default, out_dir);
    make_uncertainty_heatmaps(summary_tbl, roi_names_default, out_dir);
end

function make_qc_overview_figure(model_tbl, out_dir)
    tags_order = {'ALL','Block-01','Block-02','Block-03'};
    runs_order = {'run_1','run_2'};

    f = figure('Color','w','Position',[100 100 1500 650]);
    tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

    nexttile;
    hold on;
    plot_grouped_points(model_tbl, 'variance_explained_percent', runs_order, tags_order);
    ylabel('Variance explained (%)');
    title('DCM fit quality');
    grid on;

    nexttile;
    hold on;
    plot_grouped_points(model_tbl, 'free_energy_F', runs_order, tags_order);
    ylabel('Free energy F');
    title('Model evidence proxy');
    grid on;

    export_figure(f, fullfile(out_dir, 'fig01_dcm_qc_overview'));
    close(f);
end

function plot_grouped_points(tbl, value_col, runs_order, tags_order)
    xlabels = {};
    xpos = 0;
    colors = lines(numel(runs_order));
    for r = 1:numel(runs_order)
        for t = 1:numel(tags_order)
            xpos = xpos + 1;
            ix = strcmp(tbl.run, runs_order{r}) & strcmp(tbl.tag, tags_order{t});
            vals = tbl.(value_col)(ix);
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                jitter = (rand(size(vals))-0.5)*0.18;
                scatter(xpos + jitter, vals, 22, colors(r,:), 'filled', 'MarkerFaceAlpha',0.45);
                plot([xpos-0.25 xpos+0.25], [median(vals) median(vals)], 'k-', 'LineWidth',2);
            end
            xlabels{end+1} = sprintf('%s\n%s', runs_order{r}, tags_order{t}); %#ok<AGROW>
        end
        xpos = xpos + 0.7;
    end
    xlim([0 xpos+0.5]);
    xticks(1:(numel(runs_order)*numel(tags_order)+numel(runs_order)-1));
    xticks(find(~cellfun(@isempty,xlabels)));
    xticklabels(xlabels);
    xtickangle(45);
end

function make_prevalence_heatmaps(summary_tbl, roi_names_default, out_dir)
    if isempty(summary_tbl); return; end
    runs = unique(summary_tbl.run, 'stable');
    tags = {'Block-01','Block-02','Block-03'};
    if any(strcmp(summary_tbl.tag, 'ALL'))
        tags = [{'ALL'}, tags];
    end
    roi = get_roi_order_from_summary(summary_tbl, roi_names_default);

    nPanels = numel(runs) * numel(tags);
    f = figure('Color','w','Position',[100 100 1900 950]);
    tiledlayout(numel(runs), numel(tags), 'Padding','compact', 'TileSpacing','compact');

    for r = 1:numel(runs)
        for t = 1:numel(tags)
            nexttile;
            ix = strcmp(summary_tbl.run, runs{r}) & strcmp(summary_tbl.tag, tags{t});
            M = nan(numel(roi));
            Txt = strings(numel(roi));
            for a = 1:height(summary_tbl(ix,:))
                row = summary_tbl(ix,:);
                src = row.source{a}; tgt = row.target{a};
                jj = find(strcmp(roi, src));
                ii = find(strcmp(roi, tgt));
                if ~isempty(ii) && ~isempty(jj)
                    M(ii,jj) = row.prevalence_pp(a);
                    Txt(ii,jj) = sprintf('%d/%d\n%.3g', row.n_present_pp(a), row.n_models(a), row.mean_all(a));
                end
            end
            imagesc(M, [0 1]);
            axis image;
            colormap(gca, parula);
            cb = colorbar; cb.Label.String = 'Posterior prevalence';
            title(sprintf('%s | %s', runs{r}, tags{t}), 'Interpreter','none');
            xticks(1:numel(roi)); yticks(1:numel(roi));
            xticklabels(roi); yticklabels(roi); xtickangle(45);
            xlabel('Source'); ylabel('Target');
            add_cell_text(M, Txt, 0, 1);
        end
    end
    sgtitle('A-matrix edge prevalence: n present / N and mean A');
    export_figure(f, fullfile(out_dir, 'fig02_A_edge_prevalence_heatmaps'));
    close(f);
end

function make_uncertainty_heatmaps(summary_tbl, roi_names_default, out_dir)
    if isempty(summary_tbl); return; end
    runs = unique(summary_tbl.run, 'stable');
    tags = {'Block-01','Block-02','Block-03'};
    if any(strcmp(summary_tbl.tag, 'ALL'))
        tags = [{'ALL'}, tags];
    end
    roi = get_roi_order_from_summary(summary_tbl, roi_names_default);

    f = figure('Color','w','Position',[100 100 1900 950]);
    tiledlayout(numel(runs), numel(tags), 'Padding','compact', 'TileSpacing','compact');

    all_unc = summary_tbl.mean_sd_posterior;
    vmax = max(all_unc(isfinite(all_unc)));
    if isempty(vmax) || ~isfinite(vmax) || vmax == 0; vmax = 1; end

    for r = 1:numel(runs)
        for t = 1:numel(tags)
            nexttile;
            ix = strcmp(summary_tbl.run, runs{r}) & strcmp(summary_tbl.tag, tags{t});
            M = nan(numel(roi));
            Txt = strings(numel(roi));
            row = summary_tbl(ix,:);
            for a = 1:height(row)
                src = row.source{a}; tgt = row.target{a};
                jj = find(strcmp(roi, src));
                ii = find(strcmp(roi, tgt));
                if ~isempty(ii) && ~isempty(jj)
                    M(ii,jj) = row.mean_sd_posterior(a);
                    Txt(ii,jj) = sprintf('sd %.3g\npp %.2f', row.mean_sd_posterior(a), row.mean_pp_direction(a));
                end
            end
            imagesc(M, [0 vmax]);
            axis image;
            colormap(gca, turbo_if_available());
            cb = colorbar; cb.Label.String = 'Mean posterior SD';
            title(sprintf('%s | %s', runs{r}, tags{t}), 'Interpreter','none');
            xticks(1:numel(roi)); yticks(1:numel(roi));
            xticklabels(roi); yticklabels(roi); xtickangle(45);
            xlabel('Source'); ylabel('Target');
            add_cell_text(M, Txt, 0, vmax);
        end
    end
    sgtitle('A-matrix posterior uncertainty: mean posterior SD and mean directional posterior probability');
    export_figure(f, fullfile(out_dir, 'fig03_A_edge_uncertainty_heatmaps'));
    close(f);
end

function cmap = turbo_if_available()
    try
        cmap = turbo;
    catch
        cmap = parula;
    end
end

function roi = get_roi_order_from_summary(summary_tbl, roi_names_default)
    all_roi = unique([summary_tbl.source; summary_tbl.target], 'stable');
    roi = roi_names_default(:);
    roi = roi(ismember(roi, all_roi));
    extras = all_roi(~ismember(all_roi, roi));
    roi = [roi; extras];
end

function add_cell_text(M, Txt, vmin, vmax)
    for ii = 1:size(M,1)
        for jj = 1:size(M,2)
            if isfinite(M(ii,jj))
                val = M(ii,jj);
                normval = (val - vmin) / (vmax - vmin + eps);
                if normval > 0.55
                    txt_color = 'w';
                else
                    txt_color = 'k';
                end
                text(jj, ii, Txt(ii,jj), 'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', 'FontSize',8, 'Color',txt_color, ...
                    'FontWeight','bold');
            end
        end
    end
end

function export_figure(f, basepath)
    try
        exportgraphics(f, [basepath '.png'], 'Resolution', 300);
        exportgraphics(f, [basepath '.pdf'], 'ContentType','vector');
    catch
        saveas(f, [basepath '.png']);
        saveas(f, [basepath '.pdf']);
    end
end
