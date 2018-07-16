clear
close all

% What the script does:
% - performs longitudinal VBM preprocessing using one subject with 2 available timepoints
%   in subfolder mpm_dir (see below)
% - details on pipeline steps are described in ./trajectory/spm_traj_pre_vbm.m
%   and supplementary information
% - pipeline outputs are normalized smoothed tissue segments and
%   intermediate files written in pre_dir (see below)
% - precalculated preprocessed data (using study-templates) are 
%   provided in ./mpm_spatially_processed/VBM

% warning: code not fully tested for all options and different application scenarios
% for any questions and requests please contact gabriel.ziegler@dzne.de

cwd = pwd;

% run demo longitudinal vbm preprocessing for example subject
spm_dir  = spm('dir');
mpm_dir  = fullfile(cwd,'mpm_example_data','mpm_generated_nativespace');
pre_dir  = fullfile(cwd,'mpm_example_data','mpm_spatially_processed_new');
cov_dir  = fullfile(cwd,'covariate');

% get structural files
% note: here a synthesized T1 weighting with maximal G/W contrast based on PD and 1/R1 image 
% is used to drive deformations over timepoints and subjects (it is not 1/R1).
filesT1syn  = spm_select('FPList',fullfile(mpm_dir,'T1syn'),'^NSPN.*_T1syn_head.nii');

% initialize spm 
spm('PET');

% job specification 
outdir = {fullfile(pre_dir,'VBM')};mkdir(outdir{1});
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.outdir           = outdir;
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.imagesStructural = cellstr(filesT1syn);
                                                                    
% infos on subjects and scans; convention PP=[id,age at scan in years]
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.ppfile        = {fullfile(cov_dir,'PP','PP.mat')};
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.segmalg       = 1; % use cat segmentation 

% shoot/dartel, 1mm/1.5mm, self-generated or cat default templates 
% (use cat default shooting in example, paper results based on whole-group study-wise template)
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.tmplfile      = {fullfile(spm_dir,'toolbox','cat12','templates_1.50mm','Template_0_IXI555_MNI152_GS.nii')};

% further parameters 
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.bias          = 1; % use globally intensity normalized and bias corrected images (to prevent global scanner shift effects)
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.skullstrip    = 2; % 0 - no skullstrip, 1 - skullstrip, 2 - additioal csf stripping 
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.estimatenoise = 2; % 0 - use fixed given sd,  1 - use whole image cf. stripped, 2 - white matter sd estimate
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.noisesd       = 1; % noise sd value used for longit reg. or just additional factor if estimatenoise is false
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.wparam        = [0 1 40 10 8]; % within-subject deformation regularization parameters (adapt to plasticity/pathology/development!)
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.space         = 2; % normalize to mni 
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.modulate      = 1; % use between-subject modulation (to include baseline differences on top of within-subject changes)
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.fwhm          = [6];
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbm.subjects      = 1; % which subjects (provided in PP.mat) to be processed (e.g. 1,2, 1:n) 

try
        spm_jobman('run', matlabbatch);
catch
        disp('*************crash preprocessing*********');
end
