%% ============================================================
%  Build CONN-ready CSV (event-count based super-cycles)
%  Condition column = INTEGER (CONN compliant)
%  Block duration extended by +1s to cover final trial
%% ============================================================

clear; clc;

cfg = project_paths();
rootDir = cfg.data_root;
outCSV = cfg.block_timing_csv;

stimCats = { ...
    'Words', ...
    'Objects', ...
    'Scrambled objects', ...
    'Consonant strings' };

nSC = 6;
stimPerSC = 64;

blocks = { [1 2], [3 4], [5 6] };

% CONN condition identifiers (integer values).
% 1 = FC_Block1
% 2 = FC_Block2
% 3 = FC_Block3
blockIDs = [1 2 3];

blockEndPadding = 1.0;   % seconds

%% ---- Open CSV
fid = fopen(outCSV,'w');
fprintf(fid,'subject,session,condition,onset,duration\n');

subDirs = dir(fullfile(rootDir,'sub-*'));
subDirs = subDirs([subDirs.isdir]);

for s = 1:numel(subDirs)

    subj = subDirs(s).name;
    funcDir = fullfile(rootDir,subj,'func');
    if ~exist(funcDir,'dir'), continue; end

    runDirs = dir(fullfile(funcDir,'run_*'));
    runDirs = runDirs([runDirs.isdir]);

    for r = 1:numel(runDirs)

        runName = runDirs(r).name;
        session = sscanf(runName,'run_%d');

        evt = dir(fullfile(funcDir,runName,'*_events.tsv'));
        if isempty(evt), continue; end

        T = readtable(fullfile(evt(1).folder,evt(1).name), ...
            'FileType','text','Delimiter','\t','ReadVariableNames',true);

        onsets = T{:,1};
        types  = string(T.trial_type);

        isStim = ismember(types,stimCats);
        stimOnsets = sort(onsets(isStim));

        if numel(stimOnsets) ~= nSC*stimPerSC
            error('%s %s: expected %d stimulus, found %d', ...
                subj,runName,nSC*stimPerSC,numel(stimOnsets));
        end

        %% ---- Build super-cycles
        SC = cell(nSC,1);
        for sc = 1:nSC
            idx = (sc-1)*stimPerSC + (1:stimPerSC);
            SC{sc} = stimOnsets(idx);
        end

        %% ---- Build blocks
        for b = 1:3
            scIdx   = blocks{b};

            onsetB  = SC{scIdx(1)}(1);
            offsetB = SC{scIdx(2)}(end);

            durB = (offsetB + blockEndPadding) - onsetB;

            % Write integer condition identifiers for CONN compatibility.
            fprintf(fid,'%s,%d,%d,%.3f,%.3f\n', ...
                subj, session, blockIDs(b), onsetB, durB);
        end
    end
end

fclose(fid);
disp('CONN-compatible block timing table written successfully.');
