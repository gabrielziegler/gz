function [SPM] = spm_traj_pre_vbq(job)
% SPM preprocessing for longitudinal Voxel-based Quantification (VBQ) 
fprintf('%-40s: %30s\n','Longitudinal VBQ preprocessing started',spm('time'))  

% SPM preprocessing steps (please download and install from
%     https://www.fil.ion.ucl.ac.uk/spm/software/download/)
% (1) prepares structural images (MT, T1-MPRAGE or synthesized T1) 
%     for longitudnal registration in terms of skull stripping (if required), 
%     gloabl intensity normalization (if appropriate for certain projects)
%     can be used
%     note: code can handle available T1 e.g. MPRAGE (coregistered with MPMs), 
%     a synthesized T1 image (based on PD and 1/R1 map or just the MT image at it is) 
% (2) longitudinal registration of the structural scans using symmetric 
%     diffeomorphic registration (Ashburner & Ridgway, 2013)
% (3) warping of all MPMs from all timepoints in the midpoint average space
%     from longitudinal registration (using SPMs deformation utilities)
% (4) segmentation of midpoint average image in gray/white/csf tissue classes 
%     and normalization into template space using 
%     Computational Anatomy Toolbox for SPM (please download and install from 
%     http://www.neuro.uni-jena.de/cat/index.html); 
%     study-wise template (Ashburner, 2007, 2011) from DARTEL or
%     geodesic shooting needs to be provided or default cat template
% (5) apply between-subjects deformation to obtain spatially normalized
%     MPMs from all timepoints and midpoint average segments 
% (6) tissue weighted smoothing (Draganski et al., 2011 describes the 
%     applied smoothing procedure to account for spatial registration errors
%     and preserve quantiative values within each tissue class; one output per gray, white, csf tissue class)

% warning: code not fully tested for all options and applicatino scenarios
% for any questions and requests please contact gabriel.ziegler@dzne.de

% last update: 12/07/2018
% the output image resolution is set to 1mm 

% processes the following MPM modalities/maps 
MPM_Lab={'A','MT','R1','R2s','T1'};

% parameters for symmetric diffeomorphic registration
if isfield(job.vbq_proc_choice,'traj_pre_vbq_T1') % additional T1 is available
    longit_use_T1=1;
    job=job.vbq_proc_choice.traj_pre_vbq_T1;
    job.segmod = 3;
    % parameters for coregistration of structural T1 scans and MPMs
    % first segmenting and skull stripping the R1 data, performed before
    % coregistering the T1 to R1 and finally the header of the T1 is updated; 
    % otherwise (set coregt1=0), the header of T1w was already updated from the VBQ toolbox to coregister with 
    % all other MPM images (A,MT,R1,R2s and PDw)
    coreg_use_mod              = 4;
    estcoreg.eoptions.cost_fun = 'nmi';
    estcoreg.eoptions.sep      = [4 2];
    estcoreg.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    estcoreg.eoptions.fwhm     = [4 4];
else % no T1 available, use MT instead
    longit_use_T1=0;
    job=job.vbq_proc_choice.traj_pre_vbq_noT1;
    MPM_Lab=MPM_Lab(1:4);
end

% parameters for longitudinal registration
usemedian      = 0;                  % we might use median (John's default) instead of mean (might be more suitable for asymmetric sampling)
excl_cross     = 0;                  % ignore cross-sectional scans (for observational studies) or push them through same pipeline
intnorm        = job.bias;           % uses cats global intensity normalized/bias corr./denoised images for longitudinal registration
                                     % this makes sense in case of undesired scanner changes over time 
longit_use_mod = job.segmod+2;       % which input modality used for registration
skullstrip     = job.skullstrip;     % 0 - no skullstripping; 1 - removes only skull and background; 2 - background, skull and csf masked out 
estnoise       = job.estimatenoise;  % allows to obtain noise estimate for longit registration within white matter only
                                     % 0 - use given sd; 1 - use whole image to estimate; 2 - white matter estimate
wparam0        = job.wparam;         % regularization parameters (might prefer to use less than JA defaults designed for dementia etc.)
noisesd        = job.noisesd;
rigid_within   = any(~isfinite(wparam0));
if rigid_within
    wparam0    = Inf*wparam0;        % rigid warps (only recommended for plasticity studies)
end
bparam         = [0 0 1e8];          % 1e6 default for MPRAGE; 1e8 for MPMs since less non-intensity corr. required
ord            = [3 3 3 0 0 0];

% use log transformed jacobian determinants
logjac  = 0; % (not yet)

% parameters for shoot/dartel normalization  
% to simplify checks of normalization accuracy
pushavg = 1;    % this pushes average to the output space                         
pushjac = 0;

% parameters for segmentation 
segmalg   = job.segmalg;
ngaus     = [2,3,2,3,4,2];
outtiss   = [1,1,1,0,0,0];
modulate  = 1; % needed due to tws 
switch segmalg
    case 0 % spm12 segmentation (not yet)
     
    case 1 % cat12 segmentation (works mostly very well)
      
        % not parallel 
        catjobss.nproc = 0;
        
        % bias
        catjobss.opts    = struct('biasstr',0.25,'affreg','mni','samp',3);
       
        % segmentation
        catjobss.extopts.segmentation = struct('APP',1070,'NCstr',-Inf,'LASstr',0.5,'gcutstr',0.25,'cleanupstr',0.5,'WMHCstr',0,'WMHC',0,'restypes',struct('best',[0.5 0.3]));
               
        % registration
        catjobss.extopts.registration.darteltpm     = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_1_IXI555_MNI152.nii')};
        catjobss.extopts.registration.shootingtpm = {job.tmplfile{1}};
        catjobss.extopts.registration.regstr      = 0.5; % optimized shooting default
        
        % output resolution
        catjobss.extopts.vox = 1;
        
        % surface
        catjobss.surface = struct('pbtres',0.5,'scale_cortex',0.7,'add_parahipp',0.1,'close_parahipp',0);
        
        % admin  
        catjobss.admin = struct('ignoreErrors',1,'verb',1,'print',1);
        nonlinonly       = 0;  % (not yet) 0 - affine and nonlin modulation or 1 - nonlinear modulation only 
      
        % output of skull stripping run on each scan
        catjobss.output  = struct(...
            'GM',struct('native',0,'warped',0,'mod',0,'dartel',0),...
            'WM',struct('native',0,'warped',0,'mod',0,'dartel',0),...
            'CSF',struct('native',0,'warped',0,'mod',0,'dartel',0),...
            'atlas',struct('native',0,'dartel',0),...
            'label',struct('native',1,'warped',0,'dartel',0),...
            'bias',struct('native',intnorm,'warped',0,'dartel',0),...
            'las',struct('native',0,'warped',0,'dartel',0),...
            'jacobian',struct('warped',0),'warps',[0 0]); 
        % note: the global intensity corrected image (created in /mri/m*.nii) created here might be used
        catjobss.output.surface  = 0;
        atlases.neuromorphometrics = 0;atlases.lpba40 = 0;atlases.cobra = 0;atlases.hammers = 0;atlases.ibsr = 0;atlases.aal = 0;atlases.mori = 0;atlases.anatomy = 0;
        catjobss.output.ROImenu.atlases = atlases;
      
        % for midpoint segmentation and normalization
        catjob = catjobss;
        catjob.output  = struct(...
            'GM',struct('native',1,'warped',outtiss(1)*(1-modulate),'mod',(nonlinonly+1)*outtiss(1)*modulate,'dartel',0),...
            'WM',struct('native',1,'warped',outtiss(2)*(1-modulate),'mod',(nonlinonly+1)*outtiss(2)*modulate,'dartel',0),...
            'CSF',struct('native',1,'warped',outtiss(3)*(1-modulate),'mod',(nonlinonly+1)*outtiss(3)*modulate,'dartel',0),...
            'atlas',struct('native',0,'dartel',0),...
            'label',struct('native',1,'warped',0,'dartel',0),...
            'bias',struct('native',0,'warped',pushavg,'dartel',0),...
            'las',struct('native',0,'warped',0,'dartel',0),...
            'jacobian',struct('warped',pushjac),'warps',[1 0]);
        catjob.output.surface = 0; % use SBQ pipeline instead
        catjob.output.ROImenu.atlases = atlases;clear atlases;     
end

% size of tissue-weighted smoothing applied at the end
fwhm = job.fwhm;

% main output directory
d   = spm_file(job.outdir{1},'cpath');
if ~exist(d,'dir')
    sts = mkdir(d);
    if ~sts, error('Error creating output directory "%s".',d); end
end
outdir = job.outdir{1};
   
% check number of files
eval(['job.imagesStructural = job.images' MPM_Lab{longit_use_mod} ';']);
nscan = numel(job.imagesStructural);
for mo=1:numel(MPM_Lab)
    if eval(['nscan~=numel(job.images' MPM_Lab{mo} ')'])
        error('wrong number of files');
    end
    mpmdir{mo}   = spm_fileparts(char(eval(['job.images' MPM_Lab{mo} '(1)'])));
end
datadir = mpmdir{longit_use_mod};
cd(datadir);

% load person period file and get within subject timing
% note: this contains information about the all input files (all subjects)
% since the code can be run on one folder with all files from a study
% here the selection of the scans of a single subject is made
% files should be renamed containing an ID number 
% convention: PP=[id,age at scan]; this can be used during analysis as well
pp       = load(job.ppfile{1},'PP');
try   PP = pp.PP;catch PP = pp.pp;end
clearvars pp;

% use ages at observations to calculate timing of groupwise registration
PP(:,3) = gz_longit(PP(:,2),PP,'min');
rid     = unique(PP(:,1));
nsub    = length(rid);

% tissue probability map for all segmentations 
tpm     = fullfile(spm('Dir'),'tpm','TPM.nii');
if segmalg
    catjob.opts.tpm   = {tpm};
    catjobss.opts.tpm = {tpm};
end
Ntpm    = nifti(tpm);ntiss=size(Ntpm.dat,4);
Mmni    = spm_get_space(tpm);
for t=1:ntiss
     tissue(t).tpm    = {[tpm ',' num2str(t)]};
     tissue(t).ngaus  = ngaus(t);
     tissue(t).warped = [0 0];
end
TPM_Lab = {'gray','white','csf','bone','skin','background'};
tiss_th  = [2,3,1]; % only used in case of noise estimation

% mask the maps if outputspace is MNI
maskmni = 0;   % use mask_ICV.nii to set outside brain values to zero
mask    = spm_data_read(fullfile(spm('Dir'),'tpm','mask_ICV.nii'));

% template from shoot/dartel (create template) 
Ntemp     = nifti(job.tmplfile);
%usetiss   = zeros(1,ntiss);  % typically grey, white, csf, etc. 
%usetiss(1:Ntemp(1).dat.dim(4)-1)=1;
DIM = Ntemp(1).dat.dim(1:3);

if ~segmalg
    % affine registration of template to MNI (given by TPM.nii)
    if job.space==2 % MNI space
        [p,f,e]=fileparts(job.tmplfile{1});
        f=fullfile(p,'templ2mni_affine.mat');
        if ~exist(f)
           templ2mni_affine = Mmni/spm_klaff(Ntemp(end),Ntpm);
           save(f,'templ2mni_affine');
        else
           load(fullfile(p,'templ2mni_affine.mat'));
        end
        roptions=struct('which',[1 0],'interp',7,'wrap',[0 0 0],'mask',0,'prefix',''); 
    end
end

% which subjects to process
subjects = job.subjects;
% spm batch call of this function contains number of subject(s) to be
% processed, e.g. all subjects=1:37; or parallel usage subjects=17;

% remove existing estimations from data folder
% not applied during code development to see all intermediate files
% files_delete{1} = {'^avg_','^BiasFiel','^wj','^j_','rc1','rc2','rc3','^v_rc1avg','^y_','^wc','c1','c2','c3','c1avg','c2avg','c3avg','.*seg8.mat','^r','^str',};
% files_delete{2} = {'^num_p','^snum_p','^mavg','^avg','^c1','^c2','^c3','.*seg8.mat','^mwc1avg','^mwc2avg','^mwc3avg','rc1avg','rc2avg','^smwc1avg','^smwc2avg','^smwc3avg','^p1_','^p2_','^p3','^sp3','^sp1','^sp2','^v_rc1avg','^ww','^iy','^y_avg','^wlabels'};
% if length(subjects)==nsub
%     for i = 1:numel(files_delete{1})
%         j = cellstr(spm_select('FPList',mpmdir{longit_use_mod},files_delete{1}{i}));
%         for k = 1:numel(j)
%             spm_unlink(j{k});
%         end
%     end
%     if longit_use_T1
%         for i = 1:numel(files_delete{1})
%             j = cellstr(spm_select('FPList',mpmdir{coreg_use_mod},files_delete{1}{i}));
%             for k = 1:numel(j)
%                 spm_unlink(j{k});
%             end
%         end
%     end
% end

% removes qform warnings during changing of hdrs
warning off;

% initialize spm
spm('defaults', 'PET');

% subjects preprocessing loop
for sub = subjects
    disp(['... Processing ' sprintf('sub_%04d',rid(sub))]);
    ind = find(rid(sub)==PP(:,1));
    % exclude cross sectional data ?
    if excl_cross
        if length(ind)==1
            disp 'only one scan found, subject ignored';
            break
        end
    end  
    subdir = fullfile(outdir,[sprintf('sub_%04d',rid(sub))]);
    if ~exist(subdir,'dir'),mkdir(subdir);end
 
            % prepare structural images for longitudinal registration
            [scans_pre, prec] = spm_traj_prep_reg(job.imagesStructural(ind,:),catjobss,intnorm,skullstrip,estnoise,noisesd,subdir);
             
            % specify options for groupwise model
            tim    = PP(ind,3);
            if usemedian % how to center
                tim = tim - median(tim); % spm default 
            else
                tim = tim - mean(tim); % using mean makes often more sense e.g. for rather assymmetric timing of follow ups
            end
            wparam = kron(wparam0,1./(abs(tim)+1/365));
            sparam = round(3*abs(tim)+2);
            Nscans = nifti(scans_pre); 
            
            % =================================================================
            % groupwise registration of all scans of a subject
            output = {'rigid','wdef'};
            if pushavg,output = {output{:},'wavg'};end  
            if rigid_within,output = {output{:},'wbia'};
            else, output = {output{:},'wjac'}; end
            Out{1} = spm_groupwise_ls(Nscans, output, prec, wparam, bparam, sparam);
            
            if rigid_within % apply differential bias field correction !!! not work yet
                Nbia=nifti(Out{1}.bia); 
                for j=1:numel(Nscans)
                    Nscans(j).dat(:,:,:)=Nscans(j).dat(:,:,:)./Nbia(j).dat(:,:,:);
                end
            end
            
            % ============================================================================
            % bring all scans of all maps into subject's midpoint image space 
            clear out2 Defjob;
            out2{1}.push.weight          = {''};
            out2{1}.push.savedir.saveusr = {subdir};
            out2{1}.push.fov.file        = Out{1}.avg;
            out2{1}.push.preserve        = 0;
            out2{1}.push.fwhm            = [0 0 0];
            Defjob.out                  = out2;
            for j=1:length(ind)
                clear comp;
                comp{1}.inv.comp{1}.def     = {Out{1}.def{j}};
                comp{1}.inv.space           = {scans_pre{j}};
                Defjob.comp                 = comp;
                fnames = {};
                for mo=1:4 % loop over other maps
                    if ~rigid_within    
                          % push mpm images to avg space     
                            eval(['fnames{1} = job.images' MPM_Lab{mo} '{ind(' num2str(j) ')};']);
                            Defjob.out{1}.push.fnames = fnames;
                            Outavgspace{j,mo} = spm_deformations(Defjob); 
                    else  % rigid body case is still missing, but might be important for short term changes !!!!!!!!!!!!
                            
                    end
                end
            end
            clear Vo 
            % ==================================================================================================
            % calculate an average image across timepoints for all maps (for checks and normalization with ants) 
            if length(ind)>1    
                flags.dmtx    = 1;
                flags.mask    = 0;
                flags.interp  = -5;
                flags.dtype   = 16;
                for mo=1:4
                    for j=1:length(ind)
                        Vi{j} = Outavgspace{j,mo}.warped{1};
                    end
                    Vostr  = fullfile(subdir,['avg_' MPM_Lab{mo} sprintf('_sub_%04d',rid(sub)) '.nii']);
                    Vo{mo} = spm_imcalc(Vi,Vostr,'sum(X)/size(X,1)',flags);
                    clear Vi;
                end
            else % in this case there is no longitudinal follow ups. just rename/copy files
                clear fnames;
                for mo=1:4
                    eval(['fnames{1} = job.images' MPM_Lab{mo} '{ind(1)};']);
                    copyfile(fnames{1}(1:end),fullfile(subdir,['avg_' MPM_Lab{mo} sprintf('_sub_%04d',rid(sub)) '.nii']));
                    Vo{mo}.fname = fullfile(subdir,['avg_' MPM_Lab{mo} sprintf('_sub_%04d',rid(sub)) '.nii']);
                end
            end    
            
        if ~segmalg % spm12 segmentation & geodesic shooting (not available yet)
           
        else % CAT12 toolbox segmentation and normalization 
            
           % get and write global volumes by jacobian integration in midpoint space 
           [p,f,e]=fileparts(Out{1}.avg);
           Nwj = {}; 
           for j=1:length(ind)  % jacobians
               Nwj={Nwj{:},Out{1}.jac{j}};
           end
           Nwj=nifti(Nwj);
          
           % =====================================================================================
           % cat12 segmentation and calculation of deformation, and cf between-subjects modulation
           catjob.data{1} = Out{1}.avg;
          
           Ntiss = {};
           for t=1:3 % tissue segments from midpoint
               Ntiss = {Ntiss{:},fullfile(p,'mri',['p' num2str(t) f e])};
               spm_unlink(Ntiss{t}); % this is for new segments for each penalization prepro parameter run
           end
           try
               cat_run(catjob);
           end
           if ~exist(Ntiss{1},'file') % affine did not work ?
                 try
                       catjob.extopts.APP=0;
                       cat_run(catjob);
                 end
           end
           if ~exist(Ntiss{1},'file')
                 try
                       catjob.extopts.APP=1;
                       cat_run(catjob);
                 end
           end
           if ~exist(Ntiss{1},'file')
                 try
                       catjob.extopts.APP=2;
                       cat_run(catjob);
                 end
           end
           Ntiss = nifti(Ntiss);
            [~,vx]=spm_get_bbox(Ntiss(1));
            clear G;
            for j=1:length(ind)
                for t=1:3 
                    vol=Nwj(j).dat(:,:,:).*Ntiss(t).dat(:,:,:);
                    G(j,t)=sum(vol(:))*prod(vx); 
                end
            end          
            G(:,4)=sum(G,2);
            save(fullfile(subdir,'global_volumes.mat'),'G');  % might be useful as covariates here and there         
            clear G Nwj Ntiss;
           
            % =================================================================
            % push mpms, midpoint image to template space 
            % note!! wp12 need to modulated (lin & nonlin for the tw smoothing)
            % mpms and avg don't     
            clear Defjob;
            fnames = {};[p,f,e]=fileparts(Out{1}.avg); 
            for mo=1:4
                for j=1:length(ind) 
                    fnames={fnames{:},Outavgspace{j,mo}.warped{1}};
                end
                fnames={fnames{:},Vo{mo}.fname}; % for pushing of avg image
            end        
            Defjob.comp{1}.def = {fullfile(p,'mri',['y_' f e])}; 
            Defjob.out{1}.pull.fnames = fnames; 
            Defjob.out{1}.pull.savedir.saveusr{1} = subdir;
            Defjob.out{1}.pull.interp = 5;
            Defjob.out{1}.pull.mask   = 0;
            Defjob.out{1}.pull.fwhm   = [0 0 0];
            spm_deformations(Defjob);
 
            % =================================================================
            % tissue weighted and gaussian smoothing of MPMs 
            % note: there is a diffusion filter alternative in the toolbox
            % folder
            if fwhm
                
                % warped modulated gray, white segments, ie. denominator in
                % the tisssue weighted smoothing g*w in Draganski, 2011
                [p,f,e]=spm_fileparts(Out{1}.avg);
                clear tc;
                for t=1:3
                    if nonlinonly
                         tc{t}=fullfile(p,'mri',['m0wp',num2str(t),f,e]);
                    else
                         tc{t}=fullfile(p,'mri',['mwp',num2str(t),f,e]);
                    end
                end    
                
                % template tissue classes 
                m_tc{1} = [Ntemp(end).dat.fname ',1'];
                m_tc{2} = [Ntemp(end).dat.fname ',2'];
                
                for smo=fwhm
                    for mo=1:4  % loop over mpms (except the T1)
                       for j=1:length(ind) % loop over timepoints
                            % get the warped mpm 
                            [p,f,e] = spm_fileparts(Outavgspace{j,mo}.warped{1});                        
                            ff = fullfile(p,['w',f,e]);  
                            for t=1:3 % loop over output tissue classes
                               % calculate numerator ws(phi) i.e. multiply warped segments and mpm
                                P{t} = spm_imcalc(strvcat(char(tc{t}),char(ff)),insert_pref(ff,['num_p' num2str(t) '_']),'(i1.*i2)');
                                P{t} = P{t}.fname;
                                % smooth numerator i.e. g*ws(phi)
                                N{t} = insert_pref(P{t},'s');spm_smooth(P{t},N{t},smo);
                                % smooth denominator i.e. g*w
                                M{t} = insert_pref(tc{t},'sden_');spm_smooth(tc{t},M{t},smo);
                                % divide 
                                Q{t} = spm_imcalc(strvcat(N{t},M{t},M{t}),insert_pref(P{t},['tws' num2str(smo)]),'(i1./i2).*(i3>0.0001)');
                            end
                        end
                    end
                    clear P N M Q;  
                end
            end   
            clear Outavgspace avgout tissout Out;  
         
            % delete files not needed
            % note: don't remove intermediate processing files during
            % code development
%             for i = 1:numel(files_delete{2})
%                     j = cellstr(spm_select('FPList',subdir,files_delete{2}{i}));
%                     for k = 1:numel(j)
%                         spm_unlink(j{k});
%                     end
%             end
            
        end % end of cat12 segmentation, maps normalization, and tissue weighted smoothing
    
end % end of subject loop 

% clean up data folder
% if length(subjects)==nsub
%     for i = 1:numel(files_delete{1})
%         j = cellstr(spm_select('FPList',mpmdir{longit_use_mod},files_delete{1}{i}));
%         for k = 1:numel(j)
%             spm_unlink(j{k});
%         end
%     end
%     
%     if longit_use_T1
%         for i = 1:numel(files_delete{1})
%             j = cellstr(spm_select('FPList',mpmdir{coreg_use_mod},files_delete{1}{i}));
%             for k = 1:numel(j)
%                 spm_unlink(j{k});
%             end
%         end
%     end
%     
% end

fprintf('%-40s: %30s\n','Longitudinal VBQ preprocessing completed',spm('time'))  

end

function [Xc,Xm]=gz_longit(X,pp,mode)
[~,ncol]=size(X);
% defining rids
rid=pp(1,1);
for scan=2:size(pp,1)
    if rid(end)~=pp(scan,1)
        rid=[rid;pp(scan,1)];
    end
end

% some useful operations on data in pp format
Xc=[];Xm=[];
for col=1:ncol
    mvec=[];
    vec=X(:,col);
    for sub=1:length(rid)
      ind=find(pp(:,1)==rid(sub));
      switch mode
          case 'mean'
              mvec=[mvec;mean(vec(ind))*ones(length(ind),1)];
          case 'min'
              mvec=[mvec;min(vec(ind))*ones(length(ind),1)];
          case 'max'
              mvec=[mvec;max(vec(ind))*ones(length(ind),1)];  
              
      end
    end
    vec=vec-mvec;
    Xc=[Xc;vec];Xm=[Xm;mvec];
end
end

function fout=insert_pref(f,p)
fout=strcat(spm_str_manip(f,'h'),filesep,p,spm_str_manip(f,'t'));
end
