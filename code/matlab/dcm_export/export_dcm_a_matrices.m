%% export_dcm_a_matrices.m
% Export subject-level DCM A-matrix estimates to a wide table.
% The DCM processing root and export directory are defined by project_paths.m.

clear; clc;

cfg = project_paths();
glm_root = cfg.glm_root;
out_dir = cfg.dcm_export_root;
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

roiNames = {'vOTC','lOTC','V1_L','IFG_L'};
nROI = numel(roiNames);
nEdges = nROI * nROI;

% Directed edge labels: sender_to_receiver
AedgeLabels = strings(1, nEdges);
k = 0;
for r = 1:nROI
    for c = 1:nROI
        k = k + 1;
        AedgeLabels(k) = "A_" + roiNames{c} + "_to_" + roiNames{r};
    end
end

Subj      = strings(0,1);
Run       = zeros(0,1);
Condition = strings(0,1);
BlockIdx  = zeros(0,1);
IsALL     = false(0,1);
Avals     = zeros(0, nEdges);
FilePath  = strings(0,1);
Status    = strings(0,1);

subDirs = dir(fullfile(glm_root, 'sub-*'));
subDirs = subDirs([subDirs.isdir]);

fprintf('Found %d subject folders in %s\n', numel(subDirs), glm_root);

for s = 1:numel(subDirs)
    subName = string(subDirs(s).name);
    subPath  = fullfile(subDirs(s).folder, subDirs(s).name);

    % Search only DCM mat files
    dcmFiles = dir(fullfile(subPath, '**', 'DCM*.mat'));
    if isempty(dcmFiles)
        dcmFiles = dir(fullfile(subPath, '**', 'dcm*.mat'));
    end

    if isempty(dcmFiles)
        warning('No DCM*.mat files found under %s', subPath);
        continue;
    end

    for f = 1:numel(dcmFiles)
        dcmFile = fullfile(dcmFiles(f).folder, dcmFiles(f).name);

        [runNum, condLabel, isAll, blockIdx] = parse_dcm_path(dcmFile);

        try
            S = load(dcmFile);
            DCM = extract_dcm_from_loaded_struct(S);

            if isempty(DCM) || ~isfield(DCM, 'Ep') || ~isfield(DCM.Ep, 'A')
                error('No DCM.Ep.A found');
            end

            A = double(DCM.Ep.A);
            if any(size(A) ~= [nROI nROI])
                error('A matrix is not %dx%d', nROI, nROI);
            end

            vec = reshape(A', 1, []);  % row-major export for readability

            Subj(end+1,1)      = subName; %#ok<SAGROW>
            Run(end+1,1)       = runNum; %#ok<SAGROW>
            Condition(end+1,1) = condLabel; %#ok<SAGROW>
            BlockIdx(end+1,1)   = blockIdx; %#ok<SAGROW>
            IsALL(end+1,1)      = isAll; %#ok<SAGROW>
            Avals(end+1,:)      = vec; %#ok<SAGROW>
            FilePath(end+1,1)   = string(dcmFile); %#ok<SAGROW>
            Status(end+1,1)     = "OK"; %#ok<SAGROW>

        catch ME
            warning('Failed: %s\n%s', dcmFile, ME.message);

            Subj(end+1,1)      = subName; %#ok<SAGROW>
            Run(end+1,1)       = runNum; %#ok<SAGROW>
            Condition(end+1,1) = condLabel; %#ok<SAGROW>
            BlockIdx(end+1,1)   = blockIdx; %#ok<SAGROW>
            IsALL(end+1,1)      = isAll; %#ok<SAGROW>
            Avals(end+1,:)      = NaN(1, nEdges); %#ok<SAGROW>
            FilePath(end+1,1)   = string(dcmFile); %#ok<SAGROW>
            Status(end+1,1)     = "FAILED"; %#ok<SAGROW>
        end
    end
end

if isempty(Subj)
    error('No DCM files were successfully read.');
end

% Build table
varNames = [{'Subject','Run','Condition','BlockIndex','IsALL'}, cellstr(AedgeLabels), {'FilePath','Status'}];
T = table(Subj, Run, Condition, BlockIdx, IsALL, ...
    Avals(:,1), Avals(:,2), Avals(:,3), Avals(:,4), ...
    Avals(:,5), Avals(:,6), Avals(:,7), Avals(:,8), ...
    Avals(:,9), Avals(:,10), Avals(:,11), Avals(:,12), ...
    Avals(:,13), Avals(:,14), Avals(:,15), Avals(:,16), ...
    FilePath, Status, ...
    'VariableNames', varNames);

csvFile  = fullfile(out_dir, 'DCM_A_allsubjects_wide.csv');
xlsxFile = fullfile(out_dir, 'DCM_A_allsubjects_wide.xlsx');
matFile  = fullfile(out_dir, 'DCM_A_allsubjects_wide.mat');

writetable(T, csvFile);
writetable(T, xlsxFile);
save(matFile, 'T', 'roiNames', 'AedgeLabels');

fprintf('\nDone.\nSaved:\n%s\n%s\n%s\n', csvFile, xlsxFile, matFile);

%% -------- local functions --------
function [runNum, condLabel, isAll, blockIdx] = parse_dcm_path(dcmFile)
    p = lower(string(dcmFile));

    tokRun = regexp(p, 'run[_-]?(\d+)', 'tokens', 'once');
    if isempty(tokRun)
        runNum = NaN;
    else
        runNum = str2double(tokRun{1});
    end

    if contains(p, 'all')
        condLabel = "ALL";
        isAll = true;
        blockIdx = 0;
        return;
    end

    tokBlk = regexp(p, 'block[-_]?0*(\d+)', 'tokens', 'once');
    if ~isempty(tokBlk)
        blockIdx = str2double(tokBlk{1});
        condLabel = "Block-" + sprintf('%02d', blockIdx);
        isAll = false;
    else
        blockIdx = -1;
        condLabel = "Unknown";
        isAll = false;
    end
end

function DCM = extract_dcm_from_loaded_struct(S)
    DCM = [];
    if isfield(S, 'DCM')
        DCM = S.DCM;
        return;
    end
    fn = fieldnames(S);
    for i = 1:numel(fn)
        if isstruct(S.(fn{i})) && isfield(S.(fn{i}), 'Ep') && isfield(S.(fn{i}).Ep, 'A')
            DCM = S.(fn{i});
            return;
        end
    end
end
