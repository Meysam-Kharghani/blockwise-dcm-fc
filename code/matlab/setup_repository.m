function cfg = setup_repository()
%SETUP_REPOSITORY Add repository MATLAB code and configured dependencies.

matlab_root = fileparts(mfilename('fullpath'));
addpath(genpath(matlab_root));

cfg = project_paths();
if ~isempty(cfg.spm12_dir) && exist(cfg.spm12_dir, 'dir')
    addpath(cfg.spm12_dir);
end
if ~isempty(cfg.conn_dir) && exist(cfg.conn_dir, 'dir')
    addpath(cfg.conn_dir);
end

fprintf('Repository MATLAB paths configured.\n');
fprintf('Repository root: %s\n', cfg.repo_root);
fprintf('Work root:       %s\n', cfg.work_root);
end
