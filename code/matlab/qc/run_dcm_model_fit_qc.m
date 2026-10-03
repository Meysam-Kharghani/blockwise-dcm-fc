function run_dcm_model_fit_qc(dcm_root, out_root, varargin)
%RUN_DCM_MODEL_FIT_QC Export free energy and variance explained for SPM DCMs.
%
% This QC script recursively searches for DCM*.mat files, loads each DCM, and
% exports model-level fit metrics that are useful for reporting and checking
% block-specific DCM analyses:
%   - Free energy / negative variational free energy approximation: DCM.F
%   - Variance explained / R2, preferably from DCM.R2 when present
%   - Fallback R2 estimated from observed and predicted BOLD responses
%   - Number of time points, ROIs, inputs, and basic data dimensions
%
% Example:
%   cfg = project_paths();
%   run_dcm_model_fit_qc(cfg.glm_root, fullfile(cfg.work_root, 'qc'));
%
%       'make_figures', true);
%
% If SPM is on the MATLAB path, the script can try to regenerate predicted
% responses using DCM.M.IS when stored predictions are not present.
%
% Outputs:
%   dcm_model_fit_long.csv
%   dcm_model_fit_group_summary.csv
%   dcm_model_fit_warnings.txt
%   fig01_r2_heatmap.png/.pdf
%   fig02_free_energy_heatmap.png/.pdf
%   fig03_r2_subject_scatter.png/.pdf
%
% Notes:
%   - R2 is exported as fraction and percent.
%   - If DCM.R2 already exists, that value is used first.
%   - If R2 cannot be computed, the script writes NaN and logs a warning.

    if nargin < 1 || isempty(dcm_root)
        dcm_root = pwd;
    end
    if nargin < 2 || isempty(out_root)
        out_root = fullfile(pwd, 'DCM_ModelFit_QC');
    end

    p = inputParser;
    p.addParameter('spm_path', '', @(x) ischar(x) || isstring(x));
    p.addParameter('make_figures', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('include_all', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('file_pattern', 'DCM*.mat', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    opt = p.Results;

    if ~isempty(opt.spm_path) && exist(opt.spm_path, 'dir')
        addpath(opt.spm_path);
    end

    if ~exist(out_root, 'dir')
        mkdir(out_root);
    end

    files = recursive_dir(dcm_root, char(opt.file_pattern));
    if isempty(files)
        error('No DCM files were found under: %s', dcm_root);
    end

    warnings = {};
    rows = struct([]);
    row_i = 0;

    fprintf('Found %d candidate DCM files.\n', numel(files));

    for i = 1:numel(files)
        fpath = fullfile(files(i).folder, files(i).name);
        meta = parse_dcm_path_metadata(fpath);

        if ~opt.include_all && strcmpi(meta.tag, 'ALL')
            continue;
        end

        row_i = row_i + 1;
        rows(row_i).file = string(fpath);
        rows(row_i).subject = string(meta.subject);
        rows(row_i).run = string(meta.run);
        rows(row_i).tag = string(meta.tag);
        rows(row_i).status = "OK";
        rows(row_i).error_message = "";
        rows(row_i).free_energy = NaN;
        rows(row_i).r2_fraction = NaN;
        rows(row_i).r2_percent = NaN;
        rows(row_i).r2_source = "missing";
        rows(row_i).r2_min_roi = NaN;
        rows(row_i).r2_max_roi = NaN;
        rows(row_i).n_timepoints = NaN;
        rows(row_i).n_rois = NaN;
        rows(row_i).n_inputs = NaN;
        rows(row_i).has_Ep = false;
        rows(row_i).has_Cp = false;
        rows(row_i).has_Pp = false;
        rows(row_i).has_R2_field = false;
        rows(row_i).has_F_field = false;

        try
            DCM = load_dcm_struct(fpath);

            rows(row_i).has_Ep = isfield(DCM, 'Ep');
            rows(row_i).has_Cp = isfield(DCM, 'Cp');
            rows(row_i).has_Pp = isfield(DCM, 'Pp');
            rows(row_i).has_R2_field = isfield(DCM, 'R2');
            rows(row_i).has_F_field = isfield(DCM, 'F');

            if isfield(DCM, 'F') && ~isempty(DCM.F)
                rows(row_i).free_energy = safe_scalar(DCM.F);
            end

            [nt, nr] = get_observed_dims(DCM);
            rows(row_i).n_timepoints = nt;
            rows(row_i).n_rois = nr;
            rows(row_i).n_inputs = get_n_inputs(DCM);

            [r2_frac, r2_source, r2_per_roi] = get_dcm_r2(DCM);
            rows(row_i).r2_fraction = r2_frac;
            rows(row_i).r2_percent = 100 .* r2_frac;
            rows(row_i).r2_source = string(r2_source);
            if ~isempty(r2_per_roi) && all(isfinite(r2_per_roi))
                rows(row_i).r2_min_roi = min(r2_per_roi);
                rows(row_i).r2_max_roi = max(r2_per_roi);
            end

            if ~isfinite(rows(row_i).r2_fraction)
                warnings{end+1,1} = sprintf('R2 could not be computed: %s', fpath); %#ok<AGROW>
            end
            if ~isfinite(rows(row_i).free_energy)
                warnings{end+1,1} = sprintf('Free energy missing: %s', fpath); %#ok<AGROW>
            end

        catch ME
            rows(row_i).status = "ERROR";
            rows(row_i).error_message = string(ME.message);
            warnings{end+1,1} = sprintf('ERROR in %s: %s', fpath, ME.message); %#ok<AGROW>
        end
    end

    T = struct2table(rows);
    T = sortrows(T, {'subject','run','tag','file'});

    out_long = fullfile(out_root, 'dcm_model_fit_long.csv');
    writetable(T, out_long);

    G = group_fit_summary(T);
    out_group = fullfile(out_root, 'dcm_model_fit_group_summary.csv');
    writetable(G, out_group);

    warn_file = fullfile(out_root, 'dcm_model_fit_warnings.txt');
    fid = fopen(warn_file, 'w');
    if fid > 0
        if isempty(warnings)
            fprintf(fid, 'No warnings.\n');
        else
            for i = 1:numel(warnings)
                fprintf(fid, '%s\n', warnings{i});
            end
        end
        fclose(fid);
    end

    fprintf('\nSaved:\n  %s\n  %s\n  %s\n', out_long, out_group, warn_file);

    if opt.make_figures
        make_fit_figures(T, out_root);
    end
end

%% ------------------------------------------------------------------------
function files = recursive_dir(root_dir, pattern)
    % Uses ** when available; falls back to manual recursion otherwise.
    try
        files = dir(fullfile(root_dir, '**', pattern));
        files = files(~[files.isdir]);
    catch
        files = manual_recursive_dir(root_dir, pattern);
    end

    % Some MATLAB versions return empty for ** if root has unusual chars.
    if isempty(files)
        files = manual_recursive_dir(root_dir, pattern);
    end
end

function files = manual_recursive_dir(root_dir, pattern)
    files = dir(fullfile(root_dir, pattern));
    files = files(~[files.isdir]);
    d = dir(root_dir);
    d = d([d.isdir]);
    for i = 1:numel(d)
        name = d(i).name;
        if strcmp(name, '.') || strcmp(name, '..')
            continue;
        end
        files = [files; manual_recursive_dir(fullfile(root_dir, name), pattern)]; %#ok<AGROW>
    end
end

function DCM = load_dcm_struct(fpath)
    S = load(fpath);
    if isfield(S, 'DCM')
        DCM = S.DCM;
        return;
    end

    fn = fieldnames(S);
    for i = 1:numel(fn)
        x = S.(fn{i});
        if isstruct(x) && (isfield(x, 'Ep') || isfield(x, 'M') || isfield(x, 'Y'))
            DCM = x;
            return;
        end
    end
    error('No DCM-like structure found in MAT file.');
end

function meta = parse_dcm_path_metadata(fpath)
    s = char(fpath);

    subject = regexp(s, '(sub[-_]?\d+)', 'match', 'once');
    if isempty(subject)
        subject = regexp(s, '(?i)(subject[-_]?\d+)', 'match', 'once');
    end
    if isempty(subject), subject = 'unknown_subject'; end
    subject = regexprep(subject, '_', '-');

    r = regexp(s, '(?i)run[-_ ]?0?(\d+)', 'tokens', 'once');
    if isempty(r)
        r = regexp(s, '(?i)R0?(\d+)', 'tokens', 'once');
    end
    if isempty(r)
        run = 'unknown_run';
    else
        run = sprintf('run_%d', str2double(r{1}));
    end

    b = regexp(s, '(?i)block[-_ ]?0?(\d+)', 'tokens', 'once');
    if ~isempty(b)
        tag = sprintf('Block-%02d', str2double(b{1}));
    elseif ~isempty(regexp(s, '(?i)(^|[/\\_ -])ALL($|[/\\_ -])', 'once'))
        tag = 'ALL';
    else
        tag = 'unknown_tag';
    end

    meta.subject = subject;
    meta.run = run;
    meta.tag = tag;
end

function v = safe_scalar(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = x(1);
    end
end

function [nt, nr] = get_observed_dims(DCM)
    Y = [];
    if isfield(DCM, 'Y') && isfield(DCM.Y, 'y')
        Y = DCM.Y.y;
    elseif isfield(DCM, 'xY') && isfield(DCM.xY, 'y')
        Y = DCM.xY.y;
    end

    if isempty(Y)
        nt = NaN;
        nr = NaN;
    else
        [nt, nr] = size(Y);
    end
end

function nU = get_n_inputs(DCM)
    nU = NaN;
    if isfield(DCM, 'U') && isfield(DCM.U, 'u') && ~isempty(DCM.U.u)
        nU = size(DCM.U.u, 2);
    elseif isfield(DCM, 'U') && isfield(DCM.U, 'name')
        nU = numel(DCM.U.name);
    end
end

function [r2_frac, source, r2_per_roi] = get_dcm_r2(DCM)
    r2_frac = NaN;
    source = 'missing';
    r2_per_roi = [];

    % 1) Prefer DCM.R2 if SPM saved it.
    if isfield(DCM, 'R2') && ~isempty(DCM.R2)
        r2_raw = double(DCM.R2(:));
        r2_raw = r2_raw(isfinite(r2_raw));
        if ~isempty(r2_raw)
            % DCM.R2 may be stored as fraction or percent depending on code/SPM.
            if nanmedian(r2_raw) > 1.5
                r2_per_roi = r2_raw ./ 100;
            else
                r2_per_roi = r2_raw;
            end
            r2_frac = mean(r2_per_roi, 'omitnan');
            source = 'DCM.R2';
            return;
        end
    end

    % 2) Try stored prediction fields.
    [Y, X0] = get_observed_y_and_confounds(DCM);
    if isempty(Y)
        return;
    end

    Yhat = get_stored_prediction(DCM, size(Y));
    if isempty(Yhat)
        % 3) Try regenerating prediction using the SPM integration scheme.
        Yhat = try_generate_prediction(DCM, size(Y));
    end

    if isempty(Yhat)
        return;
    end

    [r2_frac, r2_per_roi] = compute_r2_from_prediction(Y, Yhat, X0);
    if isfinite(r2_frac)
        source = 'computed_from_prediction';
    end
end

function [Y, X0] = get_observed_y_and_confounds(DCM)
    Y = [];
    X0 = [];
    if isfield(DCM, 'Y')
        if isfield(DCM.Y, 'y')
            Y = double(DCM.Y.y);
        end
        if isfield(DCM.Y, 'X0')
            X0 = double(DCM.Y.X0);
        end
    elseif isfield(DCM, 'xY') && isfield(DCM.xY, 'y')
        Y = double(DCM.xY.y);
    end
end

function Yhat = get_stored_prediction(DCM, target_size)
    Yhat = [];
    candidates = {};

    if isfield(DCM, 'y'), candidates{end+1} = DCM.y; end %#ok<AGROW>
    if isfield(DCM, 'yp'), candidates{end+1} = DCM.yp; end %#ok<AGROW>
    if isfield(DCM, 'yhat'), candidates{end+1} = DCM.yhat; end %#ok<AGROW>
    if isfield(DCM, 'predicted'), candidates{end+1} = DCM.predicted; end %#ok<AGROW>
    if isfield(DCM, 'Y')
        if isfield(DCM.Y, 'yhat'), candidates{end+1} = DCM.Y.yhat; end %#ok<AGROW>
        if isfield(DCM.Y, 'yp'), candidates{end+1} = DCM.Y.yp; end %#ok<AGROW>
        if isfield(DCM.Y, 'pred'), candidates{end+1} = DCM.Y.pred; end %#ok<AGROW>
    end

    for i = 1:numel(candidates)
        y = candidates{i};
        y = normalize_prediction_array(y);
        y = align_matrix_to_size(y, target_size);
        if ~isempty(y)
            Yhat = y;
            return;
        end
    end
end

function y = normalize_prediction_array(y)
    if isempty(y)
        y = [];
        return;
    end
    if iscell(y)
        if numel(y) == 1
            y = y{1};
        else
            try
                y = cat(2, y{:});
            catch
                y = [];
                return;
            end
        end
    end
    if ~isnumeric(y)
        y = [];
        return;
    end
    y = double(y);
    if ndims(y) > 2
        y = squeeze(y);
    end
end

function y = align_matrix_to_size(y, target_size)
    if isempty(y)
        return;
    end
    if isequal(size(y), target_size)
        return;
    end
    if isequal(size(y'), target_size)
        y = y';
        return;
    end
    % If rows match but columns do not, reject. Do not silently crop ROIs.
    y = [];
end

function Yhat = try_generate_prediction(DCM, target_size)
    Yhat = [];
    try
        if isfield(DCM, 'M') && isfield(DCM, 'U') && isfield(DCM, 'Ep') && ...
                isfield(DCM.M, 'IS') && ~isempty(DCM.M.IS)
            y = feval(DCM.M.IS, DCM.Ep, DCM.M, DCM.U);
            y = normalize_prediction_array(y);
            y = align_matrix_to_size(y, target_size);
            if ~isempty(y)
                Yhat = y;
            end
        end
    catch
        Yhat = [];
    end
end

function [r2_global, r2_per_roi] = compute_r2_from_prediction(Y, Yhat, X0)
    Y = double(Y);
    Yhat = double(Yhat);

    if isempty(Y) || isempty(Yhat) || ~isequal(size(Y), size(Yhat))
        r2_global = NaN;
        r2_per_roi = [];
        return;
    end

    n = size(Y, 1);
    if ~isempty(X0) && size(X0, 1) == n && rank(X0) > 0
        R = eye(n) - X0 * pinv(X0);
        Y = R * Y;
        Yhat = R * Yhat;
    else
        Y = bsxfun(@minus, Y, mean(Y, 1, 'omitnan'));
        Yhat = bsxfun(@minus, Yhat, mean(Yhat, 1, 'omitnan'));
    end

    residual = Y - Yhat;

    rss_roi = sum(residual.^2, 1, 'omitnan');
    tss_roi = sum(bsxfun(@minus, Y, mean(Y, 1, 'omitnan')).^2, 1, 'omitnan');

    r2_per_roi = 1 - (rss_roi ./ tss_roi);
    r2_per_roi(~isfinite(r2_per_roi)) = NaN;

    rss = sum(rss_roi, 'omitnan');
    tss = sum(tss_roi, 'omitnan');
    r2_global = 1 - rss ./ tss;
    if ~isfinite(r2_global)
        r2_global = mean(r2_per_roi, 'omitnan');
    end
end

function G = group_fit_summary(T)
    ok = strcmp(T.status, 'OK');
    T2 = T(ok, :);

    if isempty(T2)
        G = table();
        return;
    end

    groups = unique(T2(:, {'run','tag'}), 'rows');

    run = strings(height(groups),1);
    tag = strings(height(groups),1);
    n_models = zeros(height(groups),1);
    n_ok = zeros(height(groups),1);
    mean_free_energy = NaN(height(groups),1);
    sd_free_energy = NaN(height(groups),1);
    mean_r2_percent = NaN(height(groups),1);
    sd_r2_percent = NaN(height(groups),1);
    median_r2_percent = NaN(height(groups),1);
    min_r2_percent = NaN(height(groups),1);
    max_r2_percent = NaN(height(groups),1);
    n_missing_r2 = zeros(height(groups),1);

    for i = 1:height(groups)
        idx = strcmp(T2.run, groups.run(i)) & strcmp(T2.tag, groups.tag(i));
        x = T2(idx,:);
        run(i) = groups.run(i);
        tag(i) = groups.tag(i);
        n_models(i) = height(x);
        n_ok(i) = sum(strcmp(x.status, 'OK'));
        mean_free_energy(i) = mean(x.free_energy, 'omitnan');
        sd_free_energy(i) = std(x.free_energy, 'omitnan');
        mean_r2_percent(i) = mean(x.r2_percent, 'omitnan');
        sd_r2_percent(i) = std(x.r2_percent, 'omitnan');
        median_r2_percent(i) = median(x.r2_percent, 'omitnan');
        min_r2_percent(i) = min(x.r2_percent, [], 'omitnan');
        max_r2_percent(i) = max(x.r2_percent, [], 'omitnan');
        n_missing_r2(i) = sum(~isfinite(x.r2_percent));
    end

    G = table(run, tag, n_models, n_ok, mean_free_energy, sd_free_energy, ...
        mean_r2_percent, sd_r2_percent, median_r2_percent, ...
        min_r2_percent, max_r2_percent, n_missing_r2);

    G = sort_by_run_tag(G);
end

function T = sort_by_run_tag(T)
    run_order = run_order_num(T.run);
    tag_order = tag_order_num(T.tag);
    [~, idx] = sortrows([run_order, tag_order]);
    T = T(idx,:);
end

function x = run_order_num(run)
    x = NaN(numel(run),1);
    for i = 1:numel(run)
        tok = regexp(char(run(i)), '(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            x(i) = str2double(tok{1});
        else
            x(i) = 999;
        end
    end
end

function x = tag_order_num(tag)
    x = NaN(numel(tag),1);
    for i = 1:numel(tag)
        s = char(tag(i));
        if strcmpi(s, 'ALL')
            x(i) = 0;
        else
            tok = regexp(s, '(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                x(i) = str2double(tok{1});
            else
                x(i) = 999;
            end
        end
    end
end

function make_fit_figures(T, out_root)
    ok = strcmp(T.status, 'OK');
    T = T(ok, :);
    if isempty(T)
        warning('No OK rows for plotting.');
        return;
    end

    runs = unique(T.run, 'stable');
    tags = {'ALL','Block-01','Block-02','Block-03'};
    tags = tags(ismember(tags, cellstr(T.tag)));

    R2 = make_run_tag_matrix(T, runs, tags, 'r2_percent');
    F = make_run_tag_matrix(T, runs, tags, 'free_energy');

    % Figure 1: R2 heatmap
    fig = figure('Color','w','Units','pixels','Position',[100 100 1200 450]);
    imagesc(R2);
    axis equal tight;
    colorbar;
    title('DCM variance explained (R^2, %)');
    set(gca, 'XTick', 1:numel(tags), 'XTickLabel', tags, ...
        'YTick', 1:numel(runs), 'YTickLabel', cellstr(runs), 'FontSize', 12);
    xlabel('Model segment');
    ylabel('Run');
    annotate_matrix(R2, '%.1f%%');
    save_figure(fig, out_root, 'fig01_r2_heatmap');

    % Figure 2: Free energy heatmap
    fig = figure('Color','w','Units','pixels','Position',[100 100 1200 450]);
    imagesc(F);
    axis equal tight;
    colorbar;
    title('DCM free energy');
    set(gca, 'XTick', 1:numel(tags), 'XTickLabel', tags, ...
        'YTick', 1:numel(runs), 'YTickLabel', cellstr(runs), 'FontSize', 12);
    xlabel('Model segment');
    ylabel('Run');
    annotate_matrix(F, '%.1f');
    save_figure(fig, out_root, 'fig02_free_energy_heatmap');

    % Figure 3: subject-level R2 scatter
    fig = figure('Color','w','Units','pixels','Position',[100 100 1400 600]);
    hold on;
    xlabels = {};
    xpos = 0;
    rng(1);
    for r = 1:numel(runs)
        for t = 1:numel(tags)
            xpos = xpos + 1;
            idx = strcmp(T.run, runs(r)) & strcmp(T.tag, tags{t});
            y = T.r2_percent(idx);
            jitter = (rand(size(y)) - 0.5) * 0.25;
            scatter(xpos + jitter, y, 24, 'filled', 'MarkerFaceAlpha', 0.55);
            mu = mean(y, 'omitnan');
            se = std(y, 'omitnan') ./ sqrt(sum(isfinite(y)));
            errorbar(xpos, mu, se, 'k', 'LineWidth', 1.8, 'CapSize', 10);
            xlabels{end+1} = sprintf('%s\n%s', char(runs(r)), tags{t}); %#ok<AGROW>
        end
    end
    yline(0, '--');
    ylabel('Variance explained (R^2, %)');
    title('Subject-level DCM model fit');
    set(gca, 'XTick', 1:numel(xlabels), 'XTickLabel', xlabels, 'FontSize', 10);
    grid on;
    box on;
    save_figure(fig, out_root, 'fig03_r2_subject_scatter');
end

function M = make_run_tag_matrix(T, runs, tags, fieldname)
    M = NaN(numel(runs), numel(tags));
    for r = 1:numel(runs)
        for t = 1:numel(tags)
            idx = strcmp(T.run, runs(r)) & strcmp(T.tag, tags{t});
            M(r,t) = mean(T.(fieldname)(idx), 'omitnan');
        end
    end
end

function annotate_matrix(M, fmt)
    clim = caxis;
    for i = 1:size(M,1)
        for j = 1:size(M,2)
            val = M(i,j);
            if ~isfinite(val), txt = 'NA'; else, txt = sprintf(fmt, val); end
            if isfinite(val)
                normv = (val - clim(1)) / max(eps, (clim(2)-clim(1)));
                if normv > 0.55
                    c = 'w';
                else
                    c = 'k';
                end
            else
                c = 'k';
            end
            text(j, i, txt, 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', 'FontWeight','bold', 'Color', c);
        end
    end
end

function save_figure(fig, out_root, stem)
    png_path = fullfile(out_root, [stem '.png']);
    pdf_path = fullfile(out_root, [stem '.pdf']);
    try
        exportgraphics(fig, png_path, 'Resolution', 300);
        exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    catch
        print(fig, png_path, '-dpng', '-r300');
        print(fig, pdf_path, '-dpdf', '-painters');
    end
    fprintf('Saved figure: %s\n', png_path);
    close(fig);
end
