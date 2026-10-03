function check_cropped_dcm_timepoints(varargin)
% =========================================================================
% check_cropped_dcm_timepoints
% =========================================================================
% Quality-assurance script for cropped block-specific DCMs.
% Verifies that Block-01/02/03 DCMs contain fewer time points than ALL.
%
% Default root is obtained from project_paths.m.
% =========================================================================

cfg = project_paths();
P = inputParser;
P.addParameter('glm_root', cfg.cropped_glm_root, @ischar);
P.parse(varargin{:});
S = P.Results;

if ~exist(S.glm_root, 'dir'), error('GLM root not found: %s', S.glm_root); end

D = dir(fullfile(S.glm_root, 'sub-*', 'GLM_*', 'DCM_*.mat'));
rows = {};
for i = 1:numel(D)
    p = fullfile(D(i).folder, D(i).name);
    try
        L = load(p, 'DCM');
        DCM = L.DCM;
        nY = NaN;
        if isfield(DCM, 'Y') && isfield(DCM.Y, 'y')
            nY = size(DCM.Y.y, 1);
        end
        [sub_id, run_id, model] = parse_path(p);
        rows(end+1,:) = {sub_id, run_id, model, nY, p, 'OK', ''}; %#ok<AGROW>
    catch ME
        [sub_id, run_id, model] = parse_path(p);
        rows(end+1,:) = {sub_id, run_id, model, NaN, p, 'FAIL', ME.message}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'subject','run','model','n_timepoints','path','status','message'});
out_csv = fullfile(S.glm_root, 'cropped_dcm_timepoint_QA.csv');
writetable(T, out_csv);

fprintf('\nSaved QA table: %s\n', out_csv);

if ~isempty(T)
    G = groupsummary(T, {'run','model'}, {'mean','median','min','max'}, 'n_timepoints');
    disp(G);
    writetable(G, fullfile(S.glm_root, 'cropped_dcm_timepoint_QA_summary.csv'));
end
end

function [sub_id, run_id, model] = parse_path(p)
    sub_id = ''; run_id = ''; model = '';
    parts = regexp(p, filesep, 'split');
    for i = 1:numel(parts)
        if startsWith(parts{i}, 'sub-'), sub_id = parts{i}; end
        tok = regexp(parts{i}, '^GLM_run_(\d+)_(.*)$', 'tokens', 'once');
        if ~isempty(tok)
            run_id = sprintf('run_%s', tok{1});
            model = tok{2};
        end
    end
end
