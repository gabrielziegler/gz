function [SPM] = spm_traj_pre_vbm(job)
% SPM preprocessing for longitudinal Voxel-based Morphometry (VBM) 
fprintf('%-40s: %30s\n','Longitudinal VBM preprocessing started',spm('time'))  

% SPM preprocessing steps (please download and install from
%     https://www.fil.ion.ucl.ac.uk/spm/software/download/)
% (1) prepares structural images (MT, T1-MPRAGE or synthesized T1) 
%     for longitudnal registration in terms of skull stripping (if required), 
%     gloabl intensity normalization (if appropriate for certain projects)
% (2) longitudinal registration of the structural scans using symmetric 
%     diffeomorphic registration (Ashburner & Ridgway, 2013)
% (3) segmentation of midpoint average image in gray/white/csf tissue classes 
%     and normalization into template space using 
%     Computational Anatomy Toolbox for SPM (please download and install from 
%     http://www.neuro.uni-jena.de/cat/index.html); 
%     study-wise template (Ashburner, 2007, 2011) from DARTEL or
%     geodesic shooting needs to be provided or default cat template
%     can be used
% (4) modulation within-subject and cf. between-subject (if chosen)
%     accounts for volume changes and differences induced by registration 
%     within- or between subjects
% (5) Gaussian smoothing to account for registration inaccuracies?

% the general strategy used here is normalizing jacobians determinant images to template
% space (or MNI) and applying within- and between subject modulation on
% warped tissue segments of the obtained midpoint avgerage image from Ashburner & Ridgway
% 2013 longitudinal registration
%
% in this nspn application, segmention is done using cat12 toolbox (please install) only
% and normalization was performed using geodesic shooting to self-generated 
% study-wise template 

% estnoise
% 0 - use a given fixed noise sd
% 1 - estimate noise using JA's default but cf after skull/csf stripping
% 2 - use std of noise in white matter 

% skullstrip
% 0 - no skull stripping before longitudinal registration (JA's default)
% 1 - skullstrip using cat (still includes csf)
% 2 - skull and csf strip before registration (only includes gm/wm)

% warning: code not fully tested for all options and different application scenarios
% for any questions and requests please contact gabriel.ziegler@dzne.de

% last update: 12/07/2018
% the output image resolution is set to 1mm 

% which output tissue segments
outtiss   = [1,1,1,0,0,0];
indout    = find(outtiss);

% parameters for symmetric diffeomorphic registration
usemedian      = 0;                  % use median (JA's default) instead of mean (might be more suitable for asymmetric sampling)
excl_cross     = 0;                  % ignore cross-sectional scans (for observational studies) or push them through same pipeline
intnorm        = job.bias;           % uses cats global intensity normalized/bias corr./denoised images for longitudinal registration
                                     % this makes sense in case of undesired scanner changes over time 
skullstrip     = job.skullstrip;     % 0 - no skullstripping; 1 - removes only skull and background; 2 - background, skull and csf masked out 
estnoise       = job.estimatenoise;  % allows to obtain noise estimate for longit registration within white matter only
                                     % 0 - use given sd; 1 - use whole image to estimate; 2 - white matter estimate
noisesd        = job.noisesd;        % used as noise sd or additional factor when estnoise is false                                     
wparam0        = job.wparam;         % regularization parameters (might prefer to use less than JA defaults designed for dementia etc.)
rigid_within   = any(~isfinite(wparam0));
if rigid_within
    wparam0    = Inf*wparam0;        % rigid warps
end
bparam         = [0 0 1e8];          % use 1e6 default for MPRAGE; 1e8 for MPMs since less non-intensity corr. required
ord            = [3 3 3 0 0 0];

% use log transformed jacobian determinants
logjac      = 0;    % (not yet) 
pushavg     = 1;    % this pushes average to the output space 
pushjac     = 1;    % brings between-subject jacobian to output space   % does this work ??
writebias   = 0;    % (not yet) differential bias (careful with qMRI data time effects)

% parameters for shoot/dartel normalization  
% to simplify checks of normalization accuracy                           
modulate  = job.modulate;   % include between-subject variability in outputs

% parameters for segmentation 
segmalg = job.segmalg;
ngaus   = [2,3,2,3,4,2]; % numbers of Gaussians in each tissue class during segmentation

switch segmalg
    case 0 % spm12 segmentation (not implemented yet) 
     
    case 1 % cat12 segmentation 
      
        % parallel off
        catjobss.nproc = 0;
        
        % bias
        catjobss.opts  = struct('biasstr',0.25,'affreg','mni','samp',3);
       
        % segmentation
        catjobss.extopts.segmentation = struct('APP',1070,'NCstr',-Inf,'LASstr',0.5,'gcutstr',0.5,'cleanupstr',0.5,'WMHCstr',0.5,'WMHC',3,'restypes',struct('best',[0.5 0.3]));
                
        % registration
        catjobss.extopts.registration.darteltpm   = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_1_IXI555_MNI152.nii')};
        catjobss.extopts.registration.shootingtpm = {job.tmplfile{1}};
        catjobss.extopts.registration.regstr      = 0.5; % optimized default
        
        % resolution
        catjobss.extopts.vox = 1;
        
        % surface
        catjobss.surface = struct('pbtres',0.5,'scale_cortex',0.7,'add_parahipp',0.1,'close_parahipp',0);
        
        % admin  
        catjobss.admin = struct('ignoreErrors',1,'verb',1,'print',1);
        nonlinonly       = 0;   % (old stuff) 0 - affine and nonlin modulation or 1 - nonlinear modulation only 
      
        % output
        catjobss.output  = struct(...
            'GM',struct('native',0,'warped',0,'mod',0,'dartel',0),...
            'WM',struct('native',0,'warped',0,'mod',0,'dartel',0),...
            'CSF',struct('native',0,'warped',0,'mod',0,'dartel',0),...
            'atlas',struct('native',0,'dartel',0),...
            'label',struct('native',1,'warped',0,'dartel',0),...
            'bias',struct('native',intnorm,'warped',0,'dartel',0),...
            'las',struct('native',0,'warped',0,'dartel',0),...
            'jacobian',struct('warped',0),'warps',[0 0]); 
        % note: the global intensity corrected image (created in /mri/m*.nii) might be used
        catjobss.output.surface  = 0; 
        atlases.neuromorphometrics = 0;atlases.lpba40 = 0;atlases.cobra = 0;atlases.hammers = 0;atlases.ibsr = 0;atlases.aal = 0;atlases.mori = 0;atlases.anatomy = 0;
        catjobss.output.ROImenu.atlases = atlases;
       
        % for midpoint segmentation and normalization
        catjob = catjobss;
        catjob.output  = struct(...
            'GM',struct('native',1,'warped',outtiss(1)*(1-modulate),'mod',(nonlinonly+1)*outtiss(1)*modulate,'dartel',2),...
            'WM',struct('native',1,'warped',outtiss(2)*(1-modulate),'mod',(nonlinonly+1)*outtiss(2)*modulate,'dartel',2),...
            'CSF',struct('native',1,'warped',outtiss(3)*(1-modulate),'mod',(nonlinonly+1)*outtiss(3)*modulate,'dartel',2),...
            'atlas',struct('native',0,'dartel',0),...
            'label',struct('native',1,'warped',0,'dartel',0),...
            'bias',struct('native',0,'warped',pushavg,'dartel',0),...
            'las',struct('native',0,'warped',0,'dartel',0),...
            'jacobian',struct('warped',pushjac),'warps',[1 0]);% the average is useful in template space
        catjob.output.surface = 0;
        catjob.output.ROImenu.atlases = atlases;clear atlases;     
end

% size of Gaussian smoothing kernel applied 
fwhm = job.fwhm;

% main output directory
d   = spm_file(job.outdir{1},'cpath');
if ~exist(d,'dir')
    sts = mkdir(d);
    if ~sts, error('Error creating output directory "%s".',d); end
end
outdir = job.outdir{1};

% nnumber of scans and datadir
nscan         = size(job.imagesStructural,1);
[datadir,~,~] = spm_fileparts(char(job.imagesStructural(1)));
cd(datadir);

% load person period file and get within subject timing
% note: this contains information about the all input files
% since the code can be run on one folder with all files from a study
% here the selection of the scans of a single subject is made
% files should be renamed containing an ID number 
% spm batch call of this function contains number of subject(s) to be
% processed, e.g. subject=1:10;
% convention: PP=[id,age at scan]; this can be used during analysis as well
pp       = load(job.ppfile{1},'PP');
try   PP = pp.PP;catch PP = pp.pp;end
clearvars pp;

% use ages at observations to calculate timing of groupwise registration
PP(:,3) = gz_longit(PP(:,2),PP,'min');
rid     = unique(PP(:,1));
nsub    = length(rid);

% choose tissue probability map used during cat segmentation 
tpm     = fullfile(spm('Dir'),'tpm','TPM.nii');
if segmalg
    catjob.opts.tpm = {tpm};
    catjobss.opts.tpm = {tpm};
end
Ntpm    = nifti(tpm);ntiss=size(Ntpm.dat,4);
Mmni    = spm_get_space(tpm);
for t=1:ntiss
     tissue(t).tpm    = {[tpm ',' num2str(t)]};
     tissue(t).ngaus  = ngaus(t);
     tissue(t).warped = [0 0];
end

% mask the MNI outputs
maskmni = 0;   % (old code) use mask_ICV.nii to set outside brain values to zero
mask    = spm_data_read(fullfile(spm('Dir'),'tpm','mask_ICV.nii'));

% template from shoot/dartel (create template) 
Ntemp     = nifti(job.tmplfile);
DIM       = Ntemp(1).dat.dim(1:3);
usetiss   = zeros(1,ntiss);           % typically grey, white, csf, etc. 
usetiss(1:Ntemp(1).dat.dim(4)-1)=1;   % only used for shooting template
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

% remove existing estimations from data folder 
% note: no deletion during code development to see all intermediate images 
% files_delete{1} = {'^avg_','^wj','^j_','rc1','rc2','rc3','^v_rc1avg','^y_','^wc','c1','c2','c3','c1avg','c2avg','c3avg','.*seg8.mat'};
% files_delete{2} = {'^wj','^wc1avg'};
% if length(subjects)==nsub 
%     for i = 1:numel(files_delete{1})
%         j = cellstr(spm_select('FPList',datadir,files_delete{1}{i}));
%         for k = 1:numel(j)
%             spm_unlink(j{k});
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
    ind    = find(rid(sub)==PP(:,1));
    % exclude purely cross-sectional data ? 
    if excl_cross
        if length(ind)==1
            disp 'only one scan found, subject ignored';
            continue
        end
    end
    subdir = fullfile(outdir,[sprintf('sub_%04d',rid(sub))]);
    if ~exist(subdir,'dir'),mkdir(subdir);end
        
        % prepare structural images for longitudinal registration
        [scans_pre, prec] = spm_traj_prep_reg(job.imagesStructural(ind,:),catjobss,intnorm,skullstrip,estnoise,noisesd,subdir);
        
        % specify options for groupwise model
        tim = PP(ind,3);
        if usemedian % how to center time variable
            tim = tim - median(tim);
        else
            tim = tim - mean(tim); % using mean makes often more sense e.g. for rather assymmetric timing of follow ups
        end
        wparam = kron(wparam0,1./(abs(tim)+1/365));
        sparam = round(3*abs(tim)+2);
        Nscans = nifti(scans_pre); 
        
        % groupwise registration of all scans of a subject
        % =================================================================
        output = {'rigid'};
        if pushavg,output = {output{:},'wavg'};end  
        if rigid_within||writebias,output = {output{:},'wbia','bia','wjac'};
           else, output = {output{:},'wjac'}; 
        end
        Out{1} = spm_groupwise_ls(Nscans, output, prec, wparam, bparam, sparam);
                    
        if rigid_within % (not yet) apply differential bias field correction
                Nbia=nifti(Out{1}.bia); 
                for j=1:numel(Nscans)
                    Nscans(j).dat(:,:,:)=Nscans(j).dat(:,:,:)./Nbia(j).dat(:,:,:);
                end
        end
            
        if ~segmalg % use spm12 segmentation & geodesic shooting (not implemented yet)
      
        else % use CAT12 toolbox segmentation and normalization to template in MNI space 
           
           % get and write global volumes by jacobian integration in midpoint space 
           [p,f,e]=fileparts(Out{1}.avg);
           Nwj = {}; 
           for j=1:length(ind)  % jacobians
               Nwj={Nwj{:},Out{1}.jac{j}};
           end
           Nwj=nifti(Nwj);
          
           % cat12 segmentation and calculation of deformations, and cf between-subjects modulation
           % =====================================================================================
           catjob.data{1} = Out{1}.avg;
           Ntiss = {};
           for t=1:3 % tissue segments from midpoint
               Ntiss = {Ntiss{:},fullfile(p,'mri',['p' num2str(t) f e])};
               spm_unlink(Ntiss{t}); % this is for new segments for each penalization prepro parameter run
           end
           try 
               cat_run(catjob);
           end      
           if ~exist(Ntiss{1},'file') % did it work ?
                 try
                       catjob.extopts.APP=0;
                       cat_run(catjob);
                 end
           end
           if ~exist(Ntiss{1},'file') % still not yet ?
                 try
                       catjob.extopts.APP=1;
                       cat_run(catjob);
                 end
           end
           if ~exist(Ntiss{1},'file') % still nothing ?
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
                  vol    = Nwj(j).dat(:,:,:).*Ntiss(t).dat(:,:,:);
                  G(j,t) = sum(vol(:))*prod(vx); 
              end
           end          
           G(:,4) = sum(G,2);
           save(fullfile(subdir,'global_volumes.mat'),'G');           
           clear G Nwj Ntiss;
           
           % normalize within-subject jacobian det images to template space
           % ==============================================================
           % note: because tissue classes are modulated (affine+nonlin or nonlinonly),
           % there is no need for modulation of within subject-jacobians during apply deformations step 
           clear Defjob
           fnames = {}; 
           for j=1:length(ind)  % jacobians
                 fnames={fnames{:},Out{1}.jac{j}};
           end

           [p,f,e]=fileparts(Out{1}.avg);
           Defjob.comp{1}.def = {fullfile(p,'mri',['y_' f e])}; 
           Defjob.out{1}.pull.fnames = fnames; 
           Defjob.out{1}.pull.savedir.saveusr{1} = p;
           Defjob.out{1}.pull.interp = 5;
           Defjob.out{1}.pull.mask   = 0;
           Defjob.out{1}.pull.fwhm   = [0 0 0];
           spm_deformations(Defjob);
                      
           % get within-subject jacobian det images
           % ==============================================================
           Nwj={};
           for j=1:numel(scans_pre)
              [p,f,e]=fileparts(scans_pre{j}); 
              Nwj = {Nwj{:},fullfile(p,['wj_' f e])};
           end
           Nwj=nifti(Nwj);
           
           % get tissue segments of midpoint image
           % ==============================================================
           Ntiss={};[p,f,e]=fileparts(Out{1}.avg); 
           for t=indout
              if modulate
                  if nonlinonly    
                       Ntiss={Ntiss{:},fullfile(p,'mri',['m0wp' num2str(t) f e])}; 
                  else
                       Ntiss={Ntiss{:},fullfile(p,'mri',['mwp' num2str(t) f e])};
                  end                 
              else
                 Ntiss={Ntiss{:},fullfile(p,'mri',['wp'  num2str(t) f e])};
              end
           end    
           Ntiss=nifti(Ntiss);
           
           % ==============================================================
           % within-subject modulation in (mni) template space
           for j=1:numel(Nwj)            % loop over scans of subject
               for t=indout              % loop over tissue classes
                        [~,f,e]        = fileparts(Out{1}.jac{j});
                        outfile        = ['wp' num2str(t) f(2:end) e];
                        No             = nifti;
                        No.dat         = file_array(fullfile(subdir,outfile),size(Nwj(j).dat),'float32',0,1,0);
                        No.mat         = Nwj(j).mat;
                        No.mat0        = No.mat;
                        No.mat_intent  = 'Aligned';
                        No.mat0_intent = No.mat_intent;
                        No.descrip     = sprintf('realigned tissue segments');
                        create(No);
                        
                        if logjac
                            error('not yet implemented'); 
                        else  
                            No.dat(:,:,:) = Ntiss(t).dat(:,:,:).*Nwj(j).dat(:,:,:);
                        end
                                            
                        % =================================================
                        % Gaussian smoothing
           
                        if fwhm
                            for smo=fwhm
                                spm_smooth(No.dat.fname,fullfile(subdir,['s' num2str(smo) outfile]),smo);
                            end
                        end
                        clear No;
                end
           end
            
           % move midpoint & jacobian det images to subdir (only for checks)
           % ==============================================================
           
           [p,f,e]=fileparts(Out{1}.avg); 
          
           if pushavg
                movefile(fullfile(p,'mri',['wm' f e]),fullfile(subdir,['w' f e]));
           end
%             if pushjac
%                 movefile(fullfile(p,['jac_wp1' f e]),fullfile(subdir,['jac_rp1' f e]));
%            end
            
        end % end of segmentation and normalization
            
end % end of subject loop

% clean up data folder to save disk space
% not applied during code development to check all intermediate files
% if length(subjects)==nsub
% for i = 1:numel(files_delete{1})
%     j = cellstr(spm_select('FPList',datadir,files_delete{1}{i}));
%     for k = 1:numel(j)
%         spm_unlink(j{k});
%     end
% end
% end

fprintf('%-40s: %30s\n','Longitudinal VBM preprocesing completed',spm('time'))  

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
