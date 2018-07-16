clear
close all

% What the script does:
% - Performs longitudinal VBQ preprocessing using one subject with 2 available timepoints
%   in subfolder mpm_dir (see below)
% - Details on pipeline steps are described in ./trajectory/spm_traj_pre_vbq.m
%   and supplementary information
% - Pipeline outputs are normalized qMRI parameters (MT, R1, R2* and PD)
%   using tissue weighted smoothing and intermediate files all written in pre_dir (see below)
% - precalculated preprocessed data (using study-templates) are 
%   provided in ./mpm_spatially_processed/VBQ

% warning: code not fully tested for all options and different application scenarios
% for any questions and requests please contact gabriel.ziegler@dzne.de

cwd = pwd;

% run demo longitudinal vbm preprocessing for example subject
spm_dir  = spm('dir');
mpm_dir  = fullfile(cwd,'mpm_example_data','mpm_generated_nativespace');
pre_dir  = fullfile(cwd,'mpm_example_data','mpm_spatially_processed_new');
cov_dir  = fullfile(cwd,'covariate');

% structural images used to drive registration (native space from all timepoints) 
filesT1    = spm_select('FPList',fullfile(mpm_dir,'T1syn'),'^NSPN.*_T1syn_head.nii'); 
% note: this is a synthesized T1 weighting with maximal G/W
% contrast based on PD and 1/R1 image (it is not 1/R1).

% qMRI: Multi Parameter Maps (MPM)
filesPD    = spm_select('FPList',fullfile(mpm_dir,'PD'),'^NSPN.*_PD_head.nii');
filesMT    = spm_select('FPList',fullfile(mpm_dir,'MT'),'^NSPN.*_MT_head.nii');
filesR1    = spm_select('FPList',fullfile(mpm_dir,'R1'),'^NSPN.*_R1_head.nii');
filesR2s   = spm_select('FPList',fullfile(mpm_dir,'R2s'),'^NSPN.*_R2s_head.nii');
% note: The MPM multi-echo protocol was introduced in Weiskopf and Helms 2008 
% and Weiskopf et al. 2013 for the estimation of 
% - longitudinal and transverse relaxometry measures R1 = 1/T1 and R2* = 1/T2*, 
% - proton density PD and the 
% - magnetization transfer MT saturation
% The MPM protocol typically involves 6 to 8 echoes for each of the PD-, T1- and MT-weighted acquisitions 
% in a FLASH sequence (referred to as T1w, PDw and MTw echoes, respectively).
% The MPM map generation (in native space) has been performed using the hMRI open-source
% toolbox for SPM (Tabelow et al., in prep).

% initialize spm 
spm('PET');

% job specification 
outdir = {fullfile(pre_dir,'VBQ')};mkdir(outdir{1});
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.outdir    = outdir;
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.imagesT1  = cellstr(filesT1);
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.imagesA   = cellstr(filesPD);
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.imagesMT  = cellstr(filesMT);
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.imagesR1  = cellstr(filesR1);
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.imagesR2s = cellstr(filesR2s);
            
% person period information about all files above in the sample; PP=[id,age at scan]
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.ppfile     = {fullfile(cov_dir,'PP','PP.mat')};
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.segmalg    = 1;

% shoot/dartel, 1mm/1.5mm, self-generated or cat default templates 
% (use cat default shooting in example, paper results based on whole-group study-wise template)
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.tmplfile   = {fullfile(spm_dir,'toolbox','cat12','templates_1.50mm','Template_0_IXI555_MNI152_GS.nii')};

% further parameters for processing
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.bias          = 1;   % use bias corrected, intensity normalized images for longitudinal registration
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.skullstrip    = 2;   % 0 - no stripping, 1 - skullstripping, 2 - also csf is stripped 
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.estimatenoise = 2;   % 0 - use fixed sd, 1 - use whole image to estimate with JA's method cf stripped, 2 - white matter estimate using plain sd
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.noisesd       = 1;   % noise sd value used for longit reg. or just additional factor if estimatenoise is false
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.wparam        = [0 1 40 10 8]; % modified regularization parameters (adapt to plasticity/pathology/development) 
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.space         = 2;   % mni is fine
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.fwhm          = [7]; % depends on data quality and cortex/subcortex focus; 3-7 mm usually works
matlabbatch{1}.spm.tools.trajectory.traj_pre.traj_pre_vbq.vbq_proc_choice.traj_pre_vbq_T1.subjects      = 1;   % which subjects to be processed (can be used for parallel/cluster processing) 

try
        spm_jobman('run', matlabbatch);
catch
        disp('*************crash preprocessing*********');
end
