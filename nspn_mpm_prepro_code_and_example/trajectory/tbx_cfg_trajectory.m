function trajectory = tbx_cfg_trajectory
% MATLABBATCH Configuration file for toolbox 'Trajectory'

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','trajectory')); end

%--------------------------------------------------------------------------
% Working Directory
%--------------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory for SPM.mat';
dir.help    = {'Select a directory where the SPM.mat file containing the specified design matrix will be written.'};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

%--------------------------------------------------------------------------
% Output Directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory for Subject Folders with Preprocessed Files';
outdir.help    = {'Select an existing directory where the preprocessed data will be written.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

%--------------------------------------------------------------------------
% Scans 
%--------------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Preprocessed Scans';
scans.help    = {'Select a set of registered preprocessed images. They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter  = {'image','mesh'};
scans.ufilter = '.*';
scans.num     = [1 Inf];

%--------------------------------------------------------------------------
% Structural Images 
%--------------------------------------------------------------------------
imagesStructural         = cfg_files;
imagesStructural.tag     = 'imagesStructural';
imagesStructural.name    = 'Structural Images';
imagesStructural.help    = {'Select a set of T1 weighted images which the structural normalization will be based on. Alternatively, MT images obtained from biophysical modeling using SPM VBQ toolbox can be also used. Compativility with SPM longitudinal toolbox by J. Ashburner is required.'};
imagesStructural.filter  = {'image','mesh'};
imagesStructural.ufilter = '.*';
imagesStructural.num     = [1 Inf];

%--------------------------------------------------------------------------
% T1 Images 
%--------------------------------------------------------------------------
imagesT1         = cfg_files;
imagesT1.tag     = 'imagesT1';
imagesT1.name    = 'T1 Images';
imagesT1.help    = {'Select a set of raw T1 weighted images.'};
imagesT1.filter  = {'image','mesh'};
imagesT1.ufilter = '.*';
imagesT1.num     = [1 Inf];

%--------------------------------------------------------------------------
% A Images 
%--------------------------------------------------------------------------
imagesA         = cfg_files;
imagesA.tag     = 'imagesA';
imagesA.name    = 'A Images';
imagesA.help    = {'Select a set of A images obtained from first stage of VBQ preporcessing.'};
imagesA.filter  = {'image','mesh'};
imagesA.ufilter = '.*';
imagesA.num     = [1 Inf];

%--------------------------------------------------------------------------
% MT Images 
%--------------------------------------------------------------------------
imagesMT         = cfg_files;
imagesMT.tag     = 'imagesMT';
imagesMT.name    = 'MT Images';
imagesMT.help    = {'Select a set of MT images obtained from first stage of VBQ preporcessing.'};
imagesMT.filter  = {'image','mesh'};
imagesMT.ufilter = '.*';
imagesMT.num     = [1 Inf];

%--------------------------------------------------------------------------
% R1 Images 
%--------------------------------------------------------------------------
imagesR1         = cfg_files;
imagesR1.tag     = 'imagesR1';
imagesR1.name    = 'R1 Images';
imagesR1.help    = {'Select a set of R1 images obtained from first stage of VBQ preporcessing.'};
imagesR1.filter  = {'image','mesh'};
imagesR1.ufilter = '.*';
imagesR1.num     = [1 Inf];

%--------------------------------------------------------------------------
% R2* Images 
%--------------------------------------------------------------------------
imagesR2s         = cfg_files;
imagesR2s.tag     = 'imagesR2s';
imagesR2s.name    = 'R2* Images';
imagesR2s.help    = {'Select a set of R2* images obtained from first stage of VBQ preporcessing.'};
imagesR2s.filter  = {'image','mesh'};
imagesR2s.ufilter = '.*';
imagesR2s.num     = [1 Inf];

%--------------------------------------------------------------------------
% Subject Info
%--------------------------------------------------------------------------
ppfile         = cfg_files;
ppfile.tag     = 'ppfile';
ppfile.name    = 'Sample Info about Subjects Age at Multiple Scans';
ppfile.help    = {'Here one needs to provide basic person period (PP) information about the analyzed sample. Please choose a mat file which contains a matlab variable named PP. The variable PP should contain a matrix with two columns: the ID number (first column) and the AGE of the subjects at date of the scan in years (second column). Every row of PP therefore describes one scan of every subject. Consequently, the number of rows of PP must correspond to the total numbers of scans (raw/preprocessed) included in the analysis. A simple example of four scans from two elderly subjects would be PP=[1,61.1;1,62.5;2,56.0;2,57.4]. Please not, that the order of scans chosen as input (to preprocessing/analysis) must exactly correspond to the information provided in the matrix PP. Although the files and PP can be entered in any arbitrary but corresponding order, for reasons of simplicity in subsequent analysis steps it is recommended to rename raw files to be preprocessed according to an arbitray numeric ID, e.g. projectX_rid_001_time_scan1.nii, projectX_rid_001_time_scan2.nii, projectX_rid_002_time_scan1.nii, projectX_rid_002_time_scan2.nii, etc.'};    
ppfile.filter  = 'mat';
ppfile.ufilter = 'PP*\.mat$';
ppfile.num     = [1 1];

%--------------------------------------------------------------------------
% Segmentation/Normalization Algorithm
%--------------------------------------------------------------------------
segmalg         = cfg_menu;
segmalg.tag     = 'segmalg';
segmalg.name    = 'Segmentation/Normalization Algorithm';
segmalg.val     = {1};
segmalg.help    = {'Chose algorithm for segmentation/normalization.'};
segmalg.labels  = {
             'SPM12 segmentation and geodesic shooting (not yet)'
             'CAT12 Toolbox segmentation and Dartel (default)'
}';
segmalg.values  = {0 1};

%--------------------------------------------------------------------------
% Segmentation/Normalization Modality 
%--------------------------------------------------------------------------
segmod         = cfg_menu;
segmod.tag     = 'segmod';
segmod.name    = 'Segmentation/Normalization Modality';
segmod.val     = {0};
segmod.help    = {'Missing Description'};
segmod.labels  = {
             'MT (default)'
             'R1'
}';
segmod.values  = {0 1};

%--------------------------------------------------------------------------
% Template File(s)
%--------------------------------------------------------------------------
tmplfile         = cfg_files;
tmplfile.tag     = 'tmplfile';
tmplfile.name    = 'Series of Templates generated using Run Shoot/Dartel (create templates)';
tmplfile.help    = {'Select a (the series of) sample average shaped template(s) which was obtained using SPMs geodesic shooting or dartel toolbox, i.e. Run Shoot/Dartel (create template). In order to obtain such a template one needs to use SPM12 segment with Dartel export files to generate rc1* & rc2* images wich can work with shooting/dartel. It is recommended that the template was generated using either whole sample (for smaller samples) with all follow ups or a considerable sized subsample randomized across subjects and scans.'};
tmplfile.filter = 'nifti';
tmplfile.ufilter = '.*';
tmplfile.num     = [1 Inf];

%--------------------------------------------------------------------------
% Template File(s)
%--------------------------------------------------------------------------
tmplfile2         = cfg_files;
tmplfile2.tag     = 'tmplfile';
tmplfile2.name    = 'Series of Templates generated using Dartel (create templates)';
tmplfile2.help    = {'Select a the first of a series of average shaped template(s) which was obtained using SPMs dartel toolbox and after segmentaion with CAT12, i.e. Dartel (create template) on tissue segments affine registered with the MNI (CAT12 dartel export affine). In order to obtain such a template one needs to use CAT12 segment with Dartel export affine to generate rp1*_affine & rp2*_affine images wich can work with DARTEL create templates. It is recommended that the template was generated using either whole sample (for smaller samples) with all follow ups or a considerable sized subsample randomized across subjects and scans.'};
tmplfile2.filter = 'nifti';
tmplfile2.ufilter = '.*';
tmplfile2.num     = [1 Inf];

%--------------------------------------------------------------------------
% Use global intensity normalized/denoised/biascorr. structural images for
% registration
%--------------------------------------------------------------------------
bias         = cfg_menu;
bias.tag     = 'bias';
bias.name    = 'Use Intensity Normalized/Denoised/Biascorr. Images for Longitudinal Registration';
bias.val     = {0};
bias.help    = {'This option should be used when structural images for registration are obtained e.g. with scanner updates/changes across timepoints.'};
bias.labels  = {
             'No correction (default)'
             'Intensity normalize images'
}';
bias.values  = {0 1};

%--------------------------------------------------------------------------
% Skull-Stripping Before Registration
%--------------------------------------------------------------------------
skullstrip         = cfg_menu;
skullstrip.tag     = 'skullstrip';
skullstrip.name    = 'Use skullstripped images for longitudinal registration';
skullstrip.val     = {0};
skullstrip.help    = {'This option should be used when structural images for registration are obtained in any form from MPM/FLASH protocols rather than classic MPRAGE.'};
skullstrip.labels  = {
             'No skull-stripping of images (default)'
             'Register skull-stripped images'
             'Remove backround, skull and csf'
}';
skullstrip.values  = {0 1 2};

%--------------------------------------------------------------------------
% Noise Estimate (Within-Subject)
%--------------------------------------------------------------------------
estimatenoise         = cfg_menu;
estimatenoise.tag     = 'estimatenoise';
estimatenoise.name    = 'Use Noise Estimate';
estimatenoise.val     = {1};
estimatenoise.help    = {'For classic MPRAGE T1-weighted scans, the noise can be estimated using the whole image. For MPM/FLASH protocols this is expected to be biased due to the severe noise outside the brain. It would then be better use a reasonable fixed values (from experience).'};
estimatenoise.labels  = {
            'Specify fixed noise value (no estimation)'
            'Estimate noise within the whole image (default)'
            'Estimate noise within white matter'
}';
estimatenoise.values  = {0 1 2};

%--------------------------------------------------------------------------
% Fixed Noise value
%--------------------------------------------------------------------------
noisesd         = cfg_entry;
noisesd.tag     = 'noisesd';
noisesd.name    = 'Fixed Noise SD for Longitudinal Registration';
noisesd.val     = {[5]};
noisesd.num     = [1 1];
noisesd.help    = {'Specify fixed noise standard deviation used for longitudinal registration. This value is only used if specify fixed noise option was also chosen in menu. Rough points to play with can be manually estimated using function noise=spm_noise_estimate(Scans). However, these are expected to be too high for MPM based images.'};

%--------------------------------------------------------------------------
% Warping Regularization (Within-subject Registration) 
%--------------------------------------------------------------------------
wparam         = cfg_entry;
wparam.tag     = 'wparam';
wparam.name    = 'Regulatization for Within-Subject Deformations';
wparam.val     = {[0 1 80 25 80]};
wparam.num     = [1 5];
wparam.help    = {'Specify regularization parameters used for within-subject symmetric diffeomorphic registration. Note that this is an expert option for cases when different regularization is intended, but defaults should work fairly good for many applications scenarios. Please note that the regularization is linked to the entered timing in years of the follow up scans in the person period variable PP. Entering ages in months or study time in months might result in inappropriate regularization of the warps.'};

%--------------------------------------------------------------------------
% Output Space 
%--------------------------------------------------------------------------
space         = cfg_menu;
space.tag     = 'space';
space.name    = 'Output Space';
space.val     = {2};
space.help    = {'Please chose if output images (deformations or modulated tissue segments) should be produced in native space, template space (provided by geodesic shooting template) or the MNI space. Native space applies only John Ashburners Longitudinal toolbox and Segments unsing SPM12 segmentation. This option is useful for checking if longitudinal registration works for provided structural images. Segmented images of the within subject template (c1, c2 etc.) and their resliced versions (rc1, rc2) for Template generation are also produced. When chosing the MNI space, an affine registration is used to map the provided geodesic shooting templates to the TPM.nii of SPM.'};
space.labels  = {'Template space'
                 'MNI (default)'
       };
space.values  = {1 2};

%--------------------------------------------------------------------------
% Modulation (Between-Subjects)
%--------------------------------------------------------------------------
modulate         = cfg_menu;
modulate.tag     = 'modulate';
modulate.name    = 'Between-Subjects Modulation';
modulate.val     = {0};
modulate.help    = {'Although longitudinal analysis often exploits real within-subject changes, this option can be used to additionally include variability due to age-related differences across subjects. This option is suitable if e.g. age-related effects are analyzed and the design is less balanced with respect to age. Then age-related differences across subjects obtained from the between-subjects normalization included in terms of their respective jacobian determinant of the shooting deformations. Note that, this increases the between-subject variability which can improve sensitivity to detect structural change in proper mixed-effects modeling. Otherwise, in context of fixed-effects analyses this might render the analysis less sensitive if offsets/intercepts are not modelled explicitly.'};
modulate.labels  = {
             'No  Modulation (analyze only within-subject changes, default)'
             'Use Modulation (analyze both within-and between-subject variability)'
}';
modulate.values  = {0 1};

%--------------------------------------------------------------------------
% Voxel-Based Morphometry Smoothing Filter Size
%--------------------------------------------------------------------------
fwhm_vbm         = cfg_entry;
fwhm_vbm.tag     = 'fwhm';
fwhm_vbm.name    = 'FWHM for Voxel-Based Gaussian Filtering';
fwhm_vbm.val     = {[6]};
fwhm_vbm.strtype = 'e';
fwhm_vbm.num     = [1 Inf];
fwhm_vbm.help    = {'Specify the full-width at half maximum (FWHM) of the Gaussian blurring kernel in mm. Multiple entries result in multiple outputs with corresponding filter size applied in x, y and z direction. Note that default is 0 but any ``modulated'' data will show aliasing (see eg Wikipedia), which occurs because of the way the warped images are generated. In case of VBQ preprocessing additional tissue weighted filtering with FWHM kernel size is applied to preserve values within tissue classes.'};

%--------------------------------------------------------------------------
% Voxel-Based Quantification Smoothing Filter Size
%--------------------------------------------------------------------------
fwhm_vbq         = cfg_entry;
fwhm_vbq.tag     = 'fwhm';
fwhm_vbq.name    = 'FWHM for Voxel-Based Tissue-Weighted Smoothing (Draganski, 2011)';
fwhm_vbq.val     = {[5]};
fwhm_vbq.strtype = 'e';
fwhm_vbq.num     = [1 Inf];
fwhm_vbq.help    = {'Specify the full-width at half maximum (FWHM) of the Gaussian blurring kernel in mm used for tissue weighted smoothing. Multiple entries result in multiple outputs with corresponding filter size isotropicall applied in x, y and z direction. In case of VBQ preprocessing tissue weighted filtering is applied to preserve the actual quantitative values within tissue classes.'};

%--------------------------------------------------------------------------
% Surface-Based Smoothing Filter Size
%--------------------------------------------------------------------------
fwhm_sbq         = cfg_entry;
fwhm_sbq.tag     = 'fwhm';
fwhm_sbq.name    = 'FWHM of Kernel for Surface-Based Gaussian Filtering';
fwhm_sbq.val     = {[10]};
fwhm_sbq.strtype = 'e';
fwhm_sbq.num     = [1 Inf];
fwhm_sbq.help    = {'Surface based smoothing kernel size. Multiple entries result in multiple outputs with corresponding filter size.'};

%--------------------------------------------------------------------------
% Subjects
%--------------------------------------------------------------------------
subjects         = cfg_entry;
subjects.tag     = 'subjects';
subjects.name    = 'Subject Index';
subjects.val     = {1};
subjects.strtype = 'e';
subjects.num     = [1 Inf];
subjects.help    = {'Specify subjects indices to preprocess. This could be used to process only single subjects or all serially, e.g. 1 2 3 4 ... '};

%--------------------------------------------------------------------------
% Rois
%--------------------------------------------------------------------------
rois         = cfg_entry;
rois.tag     = 'rois';
rois.name    = 'Chose Roi Numbers to include';
rois.val     = {1};
rois.strtype = 'e';
rois.num     = [1 Inf];
rois.help    = {'Missing description'};

%--------------------------------------------------------------------------
% Atlasnumber for ROI selection 
%--------------------------------------------------------------------------
atlasnumber         = cfg_entry;
atlasnumber.tag     = 'atlasnumber';
atlasnumber.name    = 'Chose the Number of the Atlasfile (if multiple provided) to chose ROIs from';
atlasnumber.val     = {1};
atlasnumber.strtype = 'e';
atlasnumber.num     = [1 Inf];
atlasnumber.help    = {'Missing description'};

%--------------------------------------------------------------------------
% System Name
%--------------------------------------------------------------------------
sysname         = cfg_entry;
sysname.tag     = 'sysname';
sysname.name    = 'System Name';
sysname.help    = {'Name the Set of ROIs'};
sysname.strtype = 's';
sysname.num     = [1 Inf];

%--------------------------------------------------------------------------
% Trajectory Preprocessing using VBM 
%--------------------------------------------------------------------------
traj_pre_vbm         = cfg_exbranch;
traj_pre_vbm.tag     = 'traj_pre_vbm';
traj_pre_vbm.name    = 'Longitudinal preprocessing for Voxel-Based Morphometry (VBM)';
traj_pre_vbm.val     = {outdir imagesStructural ppfile segmalg tmplfile bias skullstrip estimatenoise noisesd wparam space modulate fwhm_vbm subjects};
traj_pre_vbm.help    = {'This calls SPM functions to preprocess MR images (of multiple subjects and timepoints) for subsequent longitudinal analysis of Voxel-Based Morphometry (VBM) ie. normalized segment images for tissue classes are produced. The routine runs a subject by subject loop consisting of a (1) symmetric diffeomorphic registration of multiple scans per subject (Ashburner & Ridgway, 2013) (2) segmentation of the obtained midpoint average image into grey and white matter tissue classes and modulation of these segments using the within-subject jacobian determinants (3) geodesic shooting of each individuals average image to an existing template (obtained from Shoot Tools (create templates) (Ashburner, 2013)) (4) uses SPM deformation utilities to push all output segment images from each individuals midpoint average space into the template space and (5) combines modulations from within- and cf between-subjects normalization, (6) applies Gaussian smoothing, and finally (7) uses an affine registration of the shooting template to provide registered output segment images in MNI space.'};
traj_pre_vbm.prog    = @spm_traj_pre_vbm;

%--------------------------------------------------------------------------
% VBQ Trajectory Preprocessing using additional T1 
%--------------------------------------------------------------------------
traj_pre_vbq_T1         = cfg_exbranch;
traj_pre_vbq_T1.tag     = 'traj_pre_vbq_T1';
traj_pre_vbq_T1.name    = 'Longitudinal Processing of MPM Maps With Additional T1-weighted Scans';
traj_pre_vbq_T1.val     = {outdir imagesT1 imagesA imagesMT imagesR1 imagesR2s ppfile segmalg tmplfile bias skullstrip estimatenoise noisesd wparam space fwhm_vbq subjects};
traj_pre_vbq_T1.help    = {'This calls SPM functions to preprocess multiple quantitative images (of multiple subjects and timepoints) for subsequent longitudinal analysis of Voxel-Based Quantification (VBQ). The routine runs a subject by subject pipeline consisting of a (0) a coregistration of the structural scans (typically T1w, T1 or MT) and quantitative maps (eg. MT, R1, etc) (1) symmetric diffeomorphic registration of multiple structural scans per subject (2) segmentation of the obtained midpoint average (3) geodesic shooting of the midpoint average image to an existing template (obtained from Shoot (create templates)) (4) uses Deformation utilities to warp quantitative maps from native space to individuals midpoint space further into the template space (given by the Shoot Template) (5) registration of template space with MNI space (6) applies a tissue-weighted smoothing procedure (see Draganski, 2011) to all scans/timepoints/paramters to preserve values within tissue classes of grey and white matter regions.'};
traj_pre_vbq_T1.prog    = @spm_traj_pre_vbq;

%--------------------------------------------------------------------------
% VBQ Trajectory Preprocessing using no additional T1 
%--------------------------------------------------------------------------
traj_pre_vbq_noT1         = cfg_exbranch;
traj_pre_vbq_noT1.tag     = 'traj_pre_vbq_noT1';
traj_pre_vbq_noT1.name    = 'Longitudinal Processing of MPM Maps Without Additional T1-weighted Scans';
traj_pre_vbq_noT1.val     = {outdir imagesA imagesMT imagesR1 imagesR2s ppfile segmod segmalg tmplfile bias skullstrip estimatenoise noisesd wparam space fwhm_vbq subjects};
traj_pre_vbq_noT1.help    = {'This calls SPM functions to preprocess multiple quantitative images (of multiple subjects and timepoints) for subsequent longitudinal analysis of Voxel-Based Quantification (VBQ). The routine runs a subject by subject pipeline consisting of a (0) a coregistration of the structural scans (typically T1w, T1 or MT) and quantitative maps (eg. MT, R1, etc) (1) symmetric diffeomorphic registration of multiple structural scans per subject (2) segmentation of the obtained midpoint average (3) geodesic shooting of the midpoint average image to an existing template (obtained from Shoot (create templates)) (4) uses Deformation utilities to warp quantitative maps from native space to individuals midpoint space further into the template space (given by the Shoot Template) (5) registration of template space with MNI space (6) applies a tissue-weighted smoothing procedure (see Draganski, 2011) to all scans/timepoints/paramters to preserve values within tissue classes of grey and white matter regions.'};
traj_pre_vbq_noT1.prog    = @spm_traj_pre_vbq;

% ---------------------------------------------------------------------
% Choice VBQ Dataspecification with additional T1 or without
% ---------------------------------------------------------------------
vbq_proc_choice             = cfg_choice;
vbq_proc_choice.tag         = 'vbq_proc_choice';
vbq_proc_choice.name        = 'Longitudinal VBQ Preprocessing Methods';
vbq_proc_choice.help        = {'You have the option to create B1 corrected parameter maps estimated from dual flip angle FLASH experiment.'};
vbq_proc_choice.values      = {traj_pre_vbq_noT1 traj_pre_vbq_T1};


%--------------------------------------------------------------------------
% VBQ Trajectory Preprocessing using MPMs (no additional T1) 
%--------------------------------------------------------------------------
traj_pre_vbq         = cfg_exbranch;
traj_pre_vbq.tag     = 'traj_pre_vbq';
traj_pre_vbq.name    = 'Longitudinal Preprocessing for Voxel-Based Quantification (VBQ)';
traj_pre_vbq.val     = {vbq_proc_choice};
traj_pre_vbq.help    = {'This calls SPM functions to preprocess multiple quantitative images (of multiple subjects and timepoints) for subsequent longitudinal analysis of Voxel-Based Quantification (VBQ). The routine runs a subject by subject pipeline consisting of a (0) a coregistration of the structural scans (typically T1w, T1 or MT) and quantitative maps (eg. MT, R1, etc) (1) symmetric diffeomorphic registration of multiple structural scans per subject (2) segmentation of the obtained midpoint average (3) geodesic shooting of the midpoint average image to an existing template (obtained from Shoot (create templates)) (4) uses Deformation utilities to warp quantitative maps from native space to individuals midpoint space further into the template space (given by the Shoot Template) (5) registration of template space with MNI space (6) applies a tissue-weighted smoothing procedure (see Draganski, 2011) to all scans/timepoints/paramters to preserve values within tissue classes of grey and white matter regions.'};
traj_pre_vbq.prog    = @spm_traj_pre_vbq;

%--------------------------------------------------------------------------
% Trajectory Preprocessing
%--------------------------------------------------------------------------
traj_pre         = cfg_choice;
traj_pre.tag     = 'traj_pre';
traj_pre.name    = 'Trajectory Preprocessing';
traj_pre.help    = {'Missing description'};
traj_pre.values  = {traj_pre_vbm traj_pre_vbq};

%--------------------------------------------------------------------------
% Trajectory Tools
%--------------------------------------------------------------------------
trajectory         = cfg_choice;
trajectory.tag     = 'trajectory';
trajectory.name    = 'Trajectory Tools';
trajectory.help    = {'Missing description'};
trajectory.values  = {traj_pre};

%======================================================================
function fout=insert_pref(f,p)
fout=strcat(spm_str_manip(f,'h'),filesep,p,spm_str_manip(f,'t'));
%======================================================================
function chk = check_entry(job)
n1 = numel(job.images);
chk = '';
n2 = sum(~cellfun('isempty',regexp(spm_str_manip(job.images,'t'),'(^c1|^c2).*.nii')));
if n1 ~= 2,
    chk = 'Wrong input - should be c1 and c2';
end
if n2 ~= 2,
    chk = 'Wrong input - should be c1 and c2';
end



function c = unlimit(c)
try
    if isa(c, 'cfg_files')
        c.num = [0 Inf];
    end
catch e
end
try
    for i=1:numel(c.val)
        c.val{i} = unlimit(c.val{i});
    end
catch e
end

function expr=cfg_expr(c, varargin) %#ok<INUSL>
    expr = 'c';
    for i=1:size(varargin,2)
        if strcmp(class(varargin{i}), 'double')
            expr = [expr '.val{' num2str(varargin{i}) '}']; %#ok<AGROW>
        else
            v = eval([expr ';']);
            for j=1:size(v.val,2)
                if strcmp(v.val{j}.tag, varargin{i})
                    break
                end
            end
            expr = [expr '.val{' num2str(j) '}']; %#ok<AGROW>
        end
    end
    expr = expr(2:end);
    
    
function expr=cfg_expr_values(c, varargin) %#ok<INUSL>
    expr = 'c';
    for i=1:size(varargin,2)
        if strcmp(class(varargin{i}), 'double')
            expr = [expr '.values{' num2str(varargin{i}) '}']; %#ok<AGROW>
        else
            v = eval([expr ';']);
            for j=1:size(v.values,2)
                if strcmp(v.values{j}.tag, varargin{i})
                    break
                end
            end
            expr = [expr '.values{' num2str(j) '}']; %#ok<AGROW>
        end
    end
    expr = expr(2:end);

function c=cfg_set_val(c, varargin)
    expr = ['c' cfg_expr(c, varargin{1:end-1})];
    eval([expr '.val={varargin{end}};']);




