%% summarize_dcm_a_matrices.m
% ------------------------------------------------------------
% Group summary for DCM A-matrices
% Input: DCM_A_allsubjects_wide.csv from cfg.dcm_export_root.
%
% Outputs:
%   DCM_group_edge_summary.csv
%   DCM_group_edge_summary.xlsx
%   DCM_group_matrices.mat
%   Per condition:
%       DCM_group_R{run}_{condition}_mean.csv
%       DCM_group_R{run}_{condition}_median.csv
%       DCM_group_R{run}_{condition}_mode.csv
%       PNG heatmaps and histograms
% ------------------------------------------------------------

clear; clc;

cfg = project_paths();
glm_root = cfg.glm_root;
in_dir = cfg.dcm_export_root;
wideFile = fullfile(in_dir, 'DCM_A_allsubjects_wide.csv');

if ~exist(wideFile, 'file')
    error('Cannot find DCM wide table: %s', wideFile);
end

T = readtable(wideFile, 'VariableNamingRule', 'preserve');

roiNames = {'vOTC','lOTC','V1_L','IFG_L'};
nROI = numel(roiNames);

% Directed edge order must match export script
edgeVars = strings(1, nROI*nROI);
k = 0;
for r = 1:nROI
    for c = 1:nROI
        k = k + 1;
        edgeVars(k) = "A_" + roiNames{c} + "_to_" + roiNames{r};
    end
end
edgeVars = cellstr(edgeVars);
nEdges = numel(edgeVars);

conditionOrder = ["ALL","Block-01","Block-02","Block-03"];
runOrder = sort(unique(T.Run));

summaryRun   = zeros(0,1);
summaryCond  = strings(0,1);
summaryEdge  = strings(0,1);
Nvals        = zeros(0,1);
MeanVals     = zeros(0,1);
MedianVals   = zeros(0,1);
SDVals       = zeros(0,1);
MinVals      = zeros(0,1);
MaxVals      = zeros(0,1);
ModeVals     = zeros(0,1);
ModeCenter   = zeros(0,1);

groupMats = struct();

for rr = 1:numel(runOrder)
    runNum = runOrder(rr);

    for cc = 1:numel(conditionOrder)
        cond = conditionOrder(cc);

        idx = T.Run == runNum & string(T.Condition) == cond & T.Status == "OK";
        if ~any(idx)
            continue;
        end

        X = T{idx, edgeVars};
        X = double(X);

        % Collect summary stats per edge
        for e = 1:nEdges
            vals = X(:,e);
            vals = vals(isfinite(vals));

            [mMode, mCenter] = hist_mode(vals, 20);

            summaryRun(end+1,1)   = runNum; %#ok<SAGROW>
            summaryCond(end+1,1)  = cond; %#ok<SAGROW>
            summaryEdge(end+1,1)  = string(edgeVars{e}); %#ok<SAGROW>
            Nvals(end+1,1)        = numel(vals); %#ok<SAGROW>
            MeanVals(end+1,1)     = mean(vals, 'omitnan'); %#ok<SAGROW>
            MedianVals(end+1,1)   = median(vals, 'omitnan'); %#ok<SAGROW>
            SDVals(end+1,1)       = std(vals, 'omitnan'); %#ok<SAGROW>
            MinVals(end+1,1)      = min(vals); %#ok<SAGROW>
            MaxVals(end+1,1)      = max(vals); %#ok<SAGROW>
            ModeVals(end+1,1)     = mMode; %#ok<SAGROW>
            ModeCenter(end+1,1)   = mCenter; %#ok<SAGROW>
        end

        % Build mean/median/mode matrices
        meanMat   = zeros(nROI);
        medianMat = zeros(nROI);
        modeMat   = zeros(nROI);

        meanVec   = zeros(nEdges,1);
        medianVec = zeros(nEdges,1);
        modeVec   = zeros(nEdges,1);

        for e = 1:nEdges
            vals = X(:,e);
            vals = vals(isfinite(vals));
            [mMode, ~] = hist_mode(vals, 20);

            meanVec(e)   = mean(vals, 'omitnan');
            medianVec(e) = median(vals, 'omitnan');
            modeVec(e)   = mMode;
        end

        meanMat   = vec_to_dcm_matrix(meanVec, nROI);
        medianMat = vec_to_dcm_matrix(medianVec, nROI);
        modeMat   = vec_to_dcm_matrix(modeVec, nROI);

        groupKey = sprintf('R%d_%s', runNum, cond);
        groupKeySafe = matlab.lang.makeValidName(groupKey);

        groupMats.(groupKeySafe).mean   = meanMat;
        groupMats.(groupKeySafe).median = medianMat;
        groupMats.(groupKeySafe).mode   = modeMat;
        groupMats.(groupKeySafe).nSubj  = sum(idx);
        groupMats.(groupKeySafe).edges  = edgeVars;

        % Save matrices
        writematrix(meanMat,   fullfile(in_dir, sprintf('DCM_group_%s_mean.csv',   groupKey)));
        writematrix(medianMat, fullfile(in_dir, sprintf('DCM_group_%s_median.csv', groupKey)));
        writematrix(modeMat,   fullfile(in_dir, sprintf('DCM_group_%s_mode.csv',   groupKey)));

        % Plots
        plot_dcm_condition(groupKey, meanMat, medianMat, modeMat, X, edgeVars, roiNames, in_dir);
    end
end

summaryT = table(summaryRun, summaryCond, summaryEdge, Nvals, MeanVals, MedianVals, SDVals, MinVals, MaxVals, ModeVals, ModeCenter, ...
    'VariableNames', {'Run','Condition','Edge','N','Mean','Median','SD','Min','Max','ModeHist','ModeBinCenter'});

summaryCsv  = fullfile(in_dir, 'DCM_group_edge_summary.csv');
summaryXlsx = fullfile(in_dir, 'DCM_group_edge_summary.xlsx');
matFile     = fullfile(in_dir, 'DCM_group_matrices.mat');

writetable(summaryT, summaryCsv);
writetable(summaryT, summaryXlsx);
save(matFile, 'groupMats', 'summaryT', 'roiNames', 'edgeVars', 'nROI');

fprintf('\nDone.\nSaved:\n%s\n%s\n%s\n', summaryCsv, summaryXlsx, matFile);

%% -------- local functions --------
function M = vec_to_dcm_matrix(vec, nROI)
    M = reshape(vec, [nROI, nROI])';
end

function [modeVal, binCenter] = hist_mode(x, nBins)
    x = x(isfinite(x));
    if isempty(x)
        modeVal = NaN;
        binCenter = NaN;
        return;
    end
    if numel(unique(x)) == 1
        modeVal = x(1);
        binCenter = x(1);
        return;
    end
    [cnt, edges] = histcounts(x, nBins);
    [~, idx] = max(cnt);
    if isempty(idx) || idx < 1
        modeVal = NaN;
        binCenter = NaN;
    else
        binCenter = mean([edges(idx), edges(idx+1)]);
        modeVal = binCenter;
    end
end

function plot_dcm_condition(groupKey, meanMat, medianMat, modeMat, X, edgeVars, roiNames, out_dir)
    nEdges = numel(edgeVars);

    % Heatmaps
    fig1 = figure('Color','w','Visible','off','Position',[100 100 1200 400]);
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

    nexttile;
    imagesc(meanMat); axis image; colorbar;
    title([groupKey ' - Mean']);
    set(gca, 'XTick', 1:numel(roiNames), 'XTickLabel', roiNames, ...
             'YTick', 1:numel(roiNames), 'YTickLabel', roiNames);

    nexttile;
    imagesc(medianMat); axis image; colorbar;
    title([groupKey ' - Median']);
    set(gca, 'XTick', 1:numel(roiNames), 'XTickLabel', roiNames, ...
             'YTick', 1:numel(roiNames), 'YTickLabel', roiNames);

    nexttile;
    imagesc(modeMat); axis image; colorbar;
    title([groupKey ' - Mode']);
    set(gca, 'XTick', 1:numel(roiNames), 'XTickLabel', roiNames, ...
             'YTick', 1:numel(roiNames), 'YTickLabel', roiNames);

    exportgraphics(fig1, fullfile(out_dir, sprintf('DCM_group_%s_heatmaps.png', groupKey)), 'Resolution', 300);
    close(fig1);

    % Histograms
    fig2 = figure('Color','w','Visible','off','Position',[100 100 1600 1200]);
    tiledlayout(4,4,'Padding','compact','TileSpacing','compact');

    for e = 1:nEdges
        nexttile;
        vals = X(:,e);
        vals = vals(isfinite(vals));
        histogram(vals, 20);
        grid on;
        title(strrep(edgeVars{e}, '_', '\_'));
        xlabel('A');
        ylabel('Count');
    end

    exportgraphics(fig2, fullfile(out_dir, sprintf('DCM_group_%s_histograms.png', groupKey)), 'Resolution', 300);
    close(fig2);
end
