function run_peb_group_analysis()
% =========================================================================
% run_peb_group_analysis.m
%
% Purpose
%   Bayesian group analysis of block-specific DCM A-matrix parameters.
%
%   The script searches the configured GLM root for estimated DCMs:
%       GLM_run_1_block-01 / DCM_*.mat
%       GLM_run_1_block-02 / DCM_*.mat
%       GLM_run_1_block-03 / DCM_*.mat
%       GLM_run_2_block-01 / DCM_*.mat
%       GLM_run_2_block-02 / DCM_*.mat
%       GLM_run_2_block-03 / DCM_*.mat
%
%   It then runs PEB on the A-matrix parameters using a second-level design:
%       1) Mean connectivity
%       2) Late effect: Block-03 > mean(Block-01, Block-02)
%       3) Mid effect:  Block-02 > Block-01
%       4) Run effect:  Run-02 > Run-01
%       5) Late x Run interaction
%       6) Mid x Run interaction
%
%   Optional participant indicators can be added to control for repeated
%   observations from the same subject.
%
% Requirements
%   - SPM12 and the estimated DCM_*.mat files are required.
%   - The analysis targets endogenous/effective connectivity in DCM.Ep.A.
%
% Outputs
%   out_root/peb_design_table.csv
%   out_root/PEB_A_main.mat
%   out_root/BMA_A_main.mat
%   out_root/PEB_A_subjectFixed.mat          optional
%   out_root/PEB_BMA_A_flat_table.csv
%   out_root/PEB_BMA_A_<effect>_matrix.csv
%   out_root/PEB_BMA_A_<effect>_posterior_probability.csv
%   out_root/PEB_BMA_A_<effect>_heatmap.png/.pdf
% =========================================================================

clear; clc;

% ============================ PATH CONFIGURATION ==========================
cfg = project_paths();
glm_root = cfg.glm_root;
out_root = cfg.peb_root;
if ~isempty(cfg.spm12_dir) && exist(cfg.spm12_dir, 'dir')
    addpath(cfg.spm12_dir);
end

% DCM filename pattern inside each GLM folder.
dcm_pattern = 'DCM_*_4ROI.mat';

% ROI order MUST match the order used when building the DCMs.
% ROI order must match the first-level DCM specification.
roi_names = {'vOTC','lOTC','V1_L','IFG_L'};

% Which DCM parameters to analyse.
peb_field = {'A'};

% Run the main PEB/BMA without participant indicators.
% The primary PEB/BMA is estimated without participant-indicator nuisance columns.
run_main_peb = true;

% Also run a sensitivity PEB with participant indicators.
% This sensitivity model includes participant indicators because each participant contributes six DCMs.
run_subject_fixed_sensitivity = true;

% Run BMR/BMA for the participant-indicator PEB too?
% Disabled by default because participant-indicator nuisance terms substantially expand the model-reduction space.
run_bma_for_subject_fixed = false;

% Only complete subjects are kept: subject must have 2 runs x 3 blocks = 6 DCMs.
require_complete_subjects = true;

% Posterior probability threshold used only for adding stars to figures.
pp_threshold = 0.95;

% Figure settings.
make_figures = true;
fig_format = {'png','pdf'};

% ============================== START SPM ================================
if ~exist(out_root, 'dir'); mkdir(out_root); end

spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

fprintf('\n============================================================\n');
fprintf('Block-specific DCM PEB analysis\n');
fprintf('GLM root: %s\n', glm_root);
fprintf('Output  : %s\n', out_root);
fprintf('============================================================\n');

% ============================ COLLECT DCMS ===============================
meta = collect_segmentwise_dcms(glm_root, dcm_pattern, require_complete_subjects);

if isempty(meta)
    error('No complete DCM set found. Check glm_root and dcm_pattern.');
end

% Build GCM cell array of DCM filenames, one row per DCM observation.
GCM = meta.dcm_path(:);

% Build design matrix.
[X_base, Xnames_base, design_table] = build_segment_design(meta);

% Save design table for transparency.
design_csv = fullfile(out_root, 'peb_design_table.csv');
writetable(design_table, design_csv);
fprintf('Saved design table: %s\n', design_csv);

% Load DCMs into memory. If memory is tight, SPM can also use file paths in many
% versions, but loading is safer across SPM12 installations.
fprintf('\nLoading %d DCMs...\n', numel(GCM));
GCM_loaded = spm_dcm_load(GCM);

% ============================== MAIN PEB =================================
if run_main_peb
    fprintf('\n============================================================\n');
    fprintf('Running MAIN PEB/BMA: block + run effects, no subject dummies\n');
    fprintf('============================================================\n');

    M = struct();
    M.X = X_base;
    M.Xnames = Xnames_base;
    M.Q = 'all';

    [PEB_main, GCM_main] = spm_dcm_peb(GCM_loaded, M, peb_field);
    BMA_main = spm_dcm_peb_bmc(PEB_main);

    save(fullfile(out_root, 'PEB_A_main.mat'), 'PEB_main', 'GCM_main', 'M', 'meta', 'design_table', '-v7.3');
    save(fullfile(out_root, 'BMA_A_main.mat'), 'BMA_main', 'PEB_main', 'M', 'meta', 'design_table', '-v7.3');

    export_peb_bma_outputs(BMA_main, Xnames_base, roi_names, out_root, 'PEB_BMA_A_main', pp_threshold, make_figures, fig_format);
end

% ======================= PARTICIPANT-INDICATOR SENSITIVITY ========================
if run_subject_fixed_sensitivity
    fprintf('\n============================================================\n');
    fprintf('Running SENSITIVITY PEB: block + run effects + participant indicators\n');
    fprintf('============================================================\n');

    [X_sf, Xnames_sf] = add_subject_fixed_effects(X_base, Xnames_base, meta.subject);

    Msf = struct();
    Msf.X = X_sf;
    Msf.Xnames = Xnames_sf;
    Msf.Q = 'all';

    [PEB_sf, GCM_sf] = spm_dcm_peb(GCM_loaded, Msf, peb_field);
    save(fullfile(out_root, 'PEB_A_subjectFixed.mat'), 'PEB_sf', 'GCM_sf', 'Msf', 'meta', 'design_table', '-v7.3');

    if run_bma_for_subject_fixed
        BMA_sf = spm_dcm_peb_bmc(PEB_sf);
        save(fullfile(out_root, 'BMA_A_subjectFixed.mat'), 'BMA_sf', 'PEB_sf', 'Msf', 'meta', 'design_table', '-v7.3');
        export_peb_bma_outputs(BMA_sf, Xnames_sf, roi_names, out_root, 'PEB_BMA_A_subjectFixed', pp_threshold, make_figures, fig_format);
    else
        % Export posterior summaries from PEB directly, without BMA.
        export_peb_bma_outputs(PEB_sf, Xnames_sf, roi_names, out_root, 'PEB_A_subjectFixed_noBMA', pp_threshold, make_figures, fig_format);
    end
end

fprintf('\n============================================================\n');
fprintf('DONE. Results saved in:\n%s\n', out_root);
fprintf('============================================================\n');

end

% =========================================================================
% Collect block-specific DCM files
% =========================================================================
function meta = collect_segmentwise_dcms(glm_root, dcm_pattern, require_complete_subjects)

block_labels = {'block-01','block-02','block-03'};
run_nums = [1 2];

sub_dirs = dir(fullfile(glm_root, 'sub-*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);
[~, idx] = sort({sub_dirs.name});
sub_dirs = sub_dirs(idx);

rows = {};

for s = 1:numel(sub_dirs)
    sub_id = sub_dirs(s).name;
    sub_path = fullfile(sub_dirs(s).folder, sub_id);

    subject_rows = {};
    ok_count = 0;

    for r = run_nums
        for b = 1:numel(block_labels)
            block = block_labels{b};
            glm_name = sprintf('GLM_run_%d_%s', r, block);
            glm_dir = fullfile(sub_path, glm_name);
            dcm_files = dir(fullfile(glm_dir, dcm_pattern));

            if isempty(dcm_files)
                fprintf('[MISSING] %s | run %d | %s | no DCM found in %s\n', sub_id, r, block, glm_dir);
                continue;
            elseif numel(dcm_files) > 1
                fprintf('[WARN] Multiple DCMs found for %s %s; using first: %s\n', sub_id, glm_name, dcm_files(1).name);
            end

            dcm_path = fullfile(dcm_files(1).folder, dcm_files(1).name);
            subject_rows(end+1,:) = {sub_id, r, block, b, dcm_path}; %#ok<AGROW>
            ok_count = ok_count + 1;
        end
    end

    if require_complete_subjects && ok_count ~= 6
        fprintf('[EXCLUDE] %s has %d/6 segment DCMs.\n', sub_id, ok_count);
    else
        rows = [rows; subject_rows]; %#ok<AGROW>
    end
end

if isempty(rows)
    meta = [];
    return;
end

meta = table();
meta.subject = rows(:,1);
meta.run = cell2mat(rows(:,2));
meta.block = rows(:,3);
meta.block_num = cell2mat(rows(:,4));
meta.dcm_path = rows(:,5);

% Stable sort: subject, run, block.
[~, idx] = sortrows([double(categorical(meta.subject)), meta.run, meta.block_num]);
meta = meta(idx,:);

fprintf('\nCollected %d DCM observations from %d subjects.\n', height(meta), numel(unique(meta.subject)));
fprintf('Expected if complete: subjects x 2 runs x 3 blocks.\n');

end

% =========================================================================
% Build block-specific design matrix
% =========================================================================
function [X, Xnames, design_table] = build_segment_design(meta)

n = height(meta);

% Contrast codes:
% Late: Block-03 > mean(Block-01, Block-02)
%   block 1 = -0.5, block 2 = -0.5, block 3 = +1
late = zeros(n,1);
late(meta.block_num == 1) = -0.5;
late(meta.block_num == 2) = -0.5;
late(meta.block_num == 3) =  1.0;

% Mid: Block-02 > Block-01
%   block 1 = -1, block 2 = +1, block 3 = 0
mid = zeros(n,1);
mid(meta.block_num == 1) = -1;
mid(meta.block_num == 2) =  1;
mid(meta.block_num == 3) =  0;

% Run: Run-02 > Run-01
run_eff = zeros(n,1);
run_eff(meta.run == 1) = -0.5;
run_eff(meta.run == 2) =  0.5;

late_by_run = late .* run_eff;
mid_by_run  = mid  .* run_eff;

X = [ones(n,1), late, mid, run_eff, late_by_run, mid_by_run];
Xnames = {'Mean', 'Late_B3_gt_mean_B1B2', 'Mid_B2_gt_B1', 'Run2_gt_Run1', 'Late_x_Run', 'Mid_x_Run'};

% Rank check.
r = rank(X);
if r < size(X,2)
    warning('Design matrix is rank deficient: rank %d but %d columns.', r, size(X,2));
end

% Save readable design table.
design_table = meta;
design_table.Mean = X(:,1);
design_table.Late_B3_gt_mean_B1B2 = X(:,2);
design_table.Mid_B2_gt_B1 = X(:,3);
design_table.Run2_gt_Run1 = X(:,4);
design_table.Late_x_Run = X(:,5);
design_table.Mid_x_Run = X(:,6);

end

% =========================================================================
% Add participant indicators as nuisance columns
% =========================================================================
function [Xsf, Xnames_sf] = add_subject_fixed_effects(X, Xnames, subjects)

subs = unique(subjects, 'stable');
Xsf = X;
Xnames_sf = Xnames;

% Use N-1 dummy columns to avoid collinearity with intercept.
for i = 2:numel(subs)
    dummy = double(strcmp(subjects, subs{i}));
    Xsf = [Xsf, dummy]; %#ok<AGROW>
    clean_sub = regexprep(subs{i}, '[^a-zA-Z0-9_]', '_');
    Xnames_sf{end+1} = sprintf('Subj_%s', clean_sub); %#ok<AGROW>
end

r = rank(Xsf);
if r < size(Xsf,2)
    warning('Participant-indicator design matrix is rank deficient: rank %d but %d columns.', r, size(Xsf,2));
end

fprintf('Participant-indicator design: %d rows x %d columns; rank = %d.\n', size(Xsf,1), size(Xsf,2), r);

end

% =========================================================================
% Export PEB/BMA outputs into CSVs and figures
% =========================================================================
function export_peb_bma_outputs(B, Xnames, roi_names, out_root, prefix, pp_threshold, make_figures, fig_format)

if ~exist(out_root, 'dir'); mkdir(out_root); end

% Extract posterior mean/covariance/probability.
Ep = B.Ep;
Cp = B.Cp;

if isfield(B, 'Pp')
    Pp = B.Pp;
else
    Pp = [];
end

if isfield(B, 'Pnames')
    Pnames = B.Pnames;
elseif isfield(B, 'Pind') && isfield(B, 'M') && isfield(B.M, 'pE')
    Pnames = make_generic_pnames(numel(Ep));
else
    Pnames = make_generic_pnames(size(Ep,1));
end

% PEB/BMA may store Ep as matrix [parameters x covariates] or vector.
[n_param, n_cov, Ep_mat, SD_mat, Pp_mat, param_names] = normalise_peb_arrays(Ep, Cp, Pp, Pnames, Xnames);

% Flat table.
flat = table();
param_col = {};
effect_col = {};
mean_col = [];
sd_col = [];
pp_col = [];
source_col = {};
target_col = {};
row_col = [];
col_col = [];

for c = 1:n_cov
    for p = 1:n_param
        pname = param_names{p};
        [src, tgt, irow, icol] = parse_A_name(pname, roi_names);
        param_col{end+1,1} = pname; %#ok<AGROW>
        effect_col{end+1,1} = Xnames{c}; %#ok<AGROW>
        mean_col(end+1,1) = Ep_mat(p,c); %#ok<AGROW>
        sd_col(end+1,1) = SD_mat(p,c); %#ok<AGROW>
        pp_col(end+1,1) = Pp_mat(p,c); %#ok<AGROW>
        source_col{end+1,1} = src; %#ok<AGROW>
        target_col{end+1,1} = tgt; %#ok<AGROW>
        row_col(end+1,1) = irow; %#ok<AGROW>
        col_col(end+1,1) = icol; %#ok<AGROW>
    end
end

flat.parameter = param_col;
flat.effect = effect_col;
flat.posterior_mean = mean_col;
flat.posterior_sd = sd_col;
flat.posterior_probability = pp_col;
flat.source = source_col;
flat.target = target_col;
flat.target_row = row_col;
flat.source_col = col_col;

flat_csv = fullfile(out_root, sprintf('%s_flat_table.csv', prefix));
writetable(flat, flat_csv);
fprintf('Saved flat PEB/BMA table: %s\n', flat_csv);

% Export one matrix per effect for A parameters.
for c = 1:n_cov
    effect_name = sanitize_filename(Xnames{c});
    Amean = nan(numel(roi_names));
    App = nan(numel(roi_names));
    Asd = nan(numel(roi_names));

    for p = 1:n_param
        [~, ~, irow, icol] = parse_A_name(param_names{p}, roi_names);
        if ~isnan(irow) && ~isnan(icol)
            Amean(irow, icol) = Ep_mat(p,c);
            App(irow, icol) = Pp_mat(p,c);
            Asd(irow, icol) = SD_mat(p,c);
        end
    end

    write_matrix_csv(fullfile(out_root, sprintf('%s_%s_matrix.csv', prefix, effect_name)), Amean, roi_names);
    write_matrix_csv(fullfile(out_root, sprintf('%s_%s_posterior_probability.csv', prefix, effect_name)), App, roi_names);
    write_matrix_csv(fullfile(out_root, sprintf('%s_%s_posterior_sd.csv', prefix, effect_name)), Asd, roi_names);

    if make_figures && any(~isnan(Amean(:)))
        figfile_base = fullfile(out_root, sprintf('%s_%s_heatmap', prefix, effect_name));
        plot_A_effect_heatmap(Amean, App, roi_names, Xnames{c}, figfile_base, pp_threshold, fig_format);
    end
end

end

% =========================================================================
% Normalise Ep/Cp/Pp to matrices [parameters x covariates]
% =========================================================================
function [n_param, n_cov, Ep_mat, SD_mat, Pp_mat, param_names] = normalise_peb_arrays(Ep, Cp, Pp, Pnames, Xnames)

n_cov = numel(Xnames);

if ismatrix(Ep) && size(Ep,2) == n_cov
    Ep_mat = Ep;
    n_param = size(Ep_mat,1);
    param_names = get_param_names(Pnames, n_param);

    SD_mat = nan(size(Ep_mat));
    if ~isempty(Cp)
        if isequal(size(Cp), size(Ep_mat))
            SD_mat = sqrt(abs(Cp));
        elseif all(size(Cp) == [numel(Ep_mat), numel(Ep_mat)])
            v = sqrt(abs(diag(Cp)));
            SD_mat = reshape(v, size(Ep_mat));
        elseif all(size(Cp) == [n_param, n_param])
            sd = sqrt(abs(diag(Cp)));
            SD_mat = repmat(sd, 1, n_cov);
        end
    end

    if isempty(Pp)
        Pp_mat = nan(size(Ep_mat));
    elseif isequal(size(Pp), size(Ep_mat))
        Pp_mat = Pp;
    elseif numel(Pp) == numel(Ep_mat)
        Pp_mat = reshape(Pp, size(Ep_mat));
    else
        Pp_mat = nan(size(Ep_mat));
    end

else
    % Vector form. Assume parameters vary fastest within covariates.
    Ep_vec = Ep(:);
    if mod(numel(Ep_vec), n_cov) ~= 0
        warning('Cannot infer [parameter x covariate] shape cleanly. Exporting as one covariate.');
        n_cov = 1;
        Xnames = {'Effect'}; %#ok<NASGU>
    end
    n_param = numel(Ep_vec) / n_cov;
    Ep_mat = reshape(Ep_vec, n_param, n_cov);
    param_names = get_param_names(Pnames, n_param);

    SD_mat = nan(size(Ep_mat));
    if ~isempty(Cp) && all(size(Cp) == [numel(Ep_vec), numel(Ep_vec)])
        SD_mat = reshape(sqrt(abs(diag(Cp))), n_param, n_cov);
    end

    if isempty(Pp)
        Pp_mat = nan(size(Ep_mat));
    elseif numel(Pp) == numel(Ep_vec)
        Pp_mat = reshape(Pp(:), n_param, n_cov);
    else
        Pp_mat = nan(size(Ep_mat));
    end
end

end

% =========================================================================
% Parameter names helpers
% =========================================================================
function names = get_param_names(Pnames, n_param)

if iscell(Pnames) && numel(Pnames) >= n_param
    names = Pnames(1:n_param);
else
    names = make_generic_pnames(n_param);
end

end

function names = make_generic_pnames(n)

names = cell(n,1);
for i = 1:n
    names{i} = sprintf('Param_%03d', i);
end

end

% =========================================================================
% Parse A parameter name into source/target indices.
% SPM versions use slightly different naming conventions, so several regexes
% are tried. DCM A matrices are row = target, column = source.
% =========================================================================
function [source, target, irow, icol] = parse_A_name(pname, roi_names)

source = '';
target = '';
irow = NaN;
icol = NaN;

patterns = { ...
    'A\((\d+),(\d+)\)', ...
    'A\{\d+\}\((\d+),(\d+)\)', ...
    'A\s*\((\d+)\s*,\s*(\d+)\)' ...
};

for k = 1:numel(patterns)
    tok = regexp(pname, patterns{k}, 'tokens', 'once');
    if ~isempty(tok)
        irow = str2double(tok{1});
        icol = str2double(tok{2});
        break;
    end
end

% Fallback: if generic parameter names are used, leave indices as NaN.
if ~isnan(irow) && ~isnan(icol) && irow >= 1 && irow <= numel(roi_names) && icol >= 1 && icol <= numel(roi_names)
    target = roi_names{irow};
    source = roi_names{icol};
end

end

% =========================================================================
% Write matrix CSV with row/column labels
% =========================================================================
function write_matrix_csv(filename, M, roi_names)

T = array2table(M, 'VariableNames', matlab.lang.makeValidName(roi_names), 'RowNames', roi_names);
writetable(T, filename, 'WriteRowNames', true);

end

% =========================================================================
% Heatmap for PEB effects
% =========================================================================
function plot_A_effect_heatmap(Amean, App, roi_names, effect_title, figfile_base, pp_threshold, fig_format)

n = numel(roi_names);

fig = figure('Color','w', 'Position', [100 100 1100 850]);
imagesc(Amean);
axis image;
set(gca, 'XTick', 1:n, 'XTickLabel', roi_names, 'YTick', 1:n, 'YTickLabel', roi_names, 'FontSize', 12);
xlabel('Source region');
ylabel('Target region');
title(sprintf('PEB/BMA A-matrix effect: %s', strrep(effect_title, '_', '\_')), 'FontSize', 14, 'FontWeight', 'bold');

% Symmetric color limits around zero.
mx = max(abs(Amean(:)), [], 'omitnan');
if isempty(mx) || isnan(mx) || mx == 0; mx = 1; end
caxis([-mx mx]);
colormap(redblue_local(256));
cb = colorbar;
ylabel(cb, 'Posterior mean');

% Cell labels.
for i = 1:n
    for j = 1:n
        val = Amean(i,j);
        if isnan(val); continue; end
        pp = App(i,j);
        star = '';
        if ~isnan(pp) && pp >= pp_threshold
            star = '*';
        end
        txt = sprintf('%.3f%s\nP=%.2f', val, star, pp);
        % Dynamic text color for readability.
        if abs(val) > 0.55 * mx
            txt_color = [1 1 1];
        else
            txt_color = [0 0 0];
        end
        text(j, i, txt, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 9, 'Color', txt_color, 'FontWeight','bold');
    end
end

% Grid lines.
for k = 0.5:1:(n+0.5)
    line([0.5 n+0.5], [k k], 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
    line([k k], [0.5 n+0.5], 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
end

for f = 1:numel(fig_format)
    ext = fig_format{f};
    out_file = sprintf('%s.%s', figfile_base, ext);
    if strcmpi(ext, 'pdf')
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, out_file, '-dpdf', '-painters');
    else
        print(fig, out_file, ['-d' ext], '-r300');
    end
end

close(fig);

end

% =========================================================================
% Red-blue colormap without requiring external toolboxes
% =========================================================================
function cmap = redblue_local(m)

if nargin < 1; m = 256; end
x = linspace(-1, 1, m)';
r = zeros(m,1); g = zeros(m,1); b = zeros(m,1);

% Negative: blue -> white
neg = x < 0;
t = (x(neg) + 1);  % 0..1
r(neg) = t;
g(neg) = t;
b(neg) = 1;

% Positive: white -> red
pos = x >= 0;
t = x(pos);  % 0..1
r(pos) = 1;
g(pos) = 1 - t;
b(pos) = 1 - t;

cmap = [r g b];

end

% =========================================================================
% Safe filename helper
% =========================================================================
function s = sanitize_filename(s)

s = regexprep(s, '[^a-zA-Z0-9_\-]', '_');
s = regexprep(s, '_+', '_');

end
