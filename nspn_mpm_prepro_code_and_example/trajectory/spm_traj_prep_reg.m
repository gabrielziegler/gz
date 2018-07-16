function [scans_pre,prec] = spm_traj_prep_reg(scans_raw,catjobss,intnorm,skullstrip,estnoise,noisesd,outdir)
% this prepares images for longitudinal registration using spm_groupwise_ls

nts       = size(scans_raw,1);
scans_pre = scans_raw;
TPM_Lab = {'gray','white','csf','bone','skin','background'};
tiss_th  = [2,3,1]; % only used in case of noise estimation

% initial estimate of scanner noise using the whole scan (JA default) before int. norm
prec = ones(nts,1)*noisesd^(-2);
if estnoise
    prec = spm_noise_estimate(scans_raw).^(-2);
end        
% =================================================================
% intensity normalization, skullstripping, prepare noise estimate
if intnorm || skullstrip || estnoise % segmentation of each scan if necessary
                 
                 if intnorm
                    disp('... ... intensity normalization applied before longitudinal registration');
                 end
                 if skullstrip==1
                    disp('... ... skullstripping applied before longitudinal registration');
                 end 
                 if skullstrip==2
                    disp('... ... skull & csf stripping applied before longitudinal registration');
                 end
                 if estnoise==2
                    noise = zeros(nts,4); % for all tissue classes
                    disp('... ... segmentation of each scan for white matter noise estimate');
                 end
                 clear Segout;
                 for i=1:nts % loop over timepoints
                        clear matlabbatch noisemask; 
                        [pp,ff,e]        = fileparts(scans_raw{i});
                        % input images for registration and labels for
                        % masking (if available)
                        Vi               = {scans_raw{i}; ...
                                            fullfile(pp,'mri',['p0',ff,e])};
                        catjobss.data  = cellstr(scans_raw{i});
                        if ~exist(Vi{2},'file') % note that this is to safe time when re-running (careful with code changes !!!) 
                           try
                               cat_run(catjobss); % using 1070 cat 12 default
                           end
                        end
                        if ~exist(Vi{2},'file') % if affine registration does not work above, alternative modes are tried
                           try
                               catjobss.extopts.APP=1;
                               cat_run(catjobss);
                           end
                        end
                        if ~exist(Vi{2},'file')
                           try
                               catjobss.extopts.APP=2;
                               cat_run(catjobss);
                           end
                        end
                        if ~exist(Vi{2},'file')
                           try
                               catjobss.extopts.APP=0;
                               cat_run(catjobss);
                           end
                        end
                        if intnorm % intensity normalized images from mri subfolder
                            % define output scans
                            scans_pre{i} = fullfile(pp,'mri',['m' ff  e]);
                            
                            % Vi for skullstripping
                            Vi{1}= scans_pre{i};

                            % udpate the noise estimate
                            if estnoise
                                prec(i,1) = spm_noise_estimate(Vi{1})^(-2);
                            end
                        end
                        if skullstrip % creates a temporay mask file 
                                if ~exist(fullfile(pp,'msk'),'dir')
                                    mkdir(fullfile(pp,'msk'));
                                end
                                if ~exist(fullfile(pp,'str'),'dir')
                                    mkdir(fullfile(pp,'str'));
                                end
                                Vm           = fullfile(pp,'msk',[ff,'_msk' e]);
                                flags.dmtx   = 0;
                                flags.mask   = 0;
                                flags.interp = -5;
                                flags.dtype  = 16;
                                % get the smoothed brain mask 
                                spm_imcalc(Vi,Vm,'(i2>0)',flags);
                                spm_smooth(Vm,Vm,[2,2,2]);

                                % warning: the following is way to hacky and specific for one recent project the code was developed!!!!
                                if intnorm
                                    Vo = fullfile(pp,'str',[ff '_instr' e]);
                                else
                                    Vo = fullfile(pp,'str',[ff '_str' e]);
                                end
                                % deal with MT case
                                if ~isempty(findstr('_MT_',Vi{1})) % here we use the intensity normalzed image in case
                                    % truncating the MT in classes slightly (Vi
                                    % contains (cf normalized) MT and p0) and write it into Vo)
                                    % spm_imcalc(Vi,Vo,'i1.*(i2==0).*(i1<1)+i1.*(i1>0).*(i1<1).*(i2>0).*(i2<1.5)+i1.*(i2>=1.5)',flags); 
                                    switch skullstrip
                                        case 2 % removes csf
                                            spm_imcalc(Vi,Vo,'i1.*(i2-0.5).*(i2<1.5).*(i2>0) + i1.*(i2>=1.5).*(i1<4)',flags); 
                                        case 1 % leaves csf in 
                                            spm_imcalc(Vi,Vo,'i1.*(i2<1.5).*(i2>0)           + i1.*(i2>=1.5).*(i1<4)',flags); 
                                    end   
                                    Vi2{2} = Vm;
                                    Vi2{1} = Vo; 
                                    % apply mask to truncated image
                                    spm_imcalc(Vi2,Vo,'i1.*i2.*(i1>0)',flags); 
                                else % T1syn or all others
                                    switch skullstrip
                                        case 2 % removes csf
                                            spm_imcalc(Vi,Vo,'i1.*(i2-0.5).*(i2<1.5).*(i2>0) + i1.*(i2>=1.5)',flags); 
                                        case 1 % leaves csf in 
                                            spm_imcalc(Vi,Vo,'i1.*(i2<1.5).*(i2>0)           + i1.*(i2>=1.5)',flags); 
                                    end   
                                    Vi2{2} = Vm;
                                    Vi2{1} = Vo; 
                                    % apply soft mask to truncated image
                                    spm_imcalc(Vi2,Vo,'i1.*i2',flags); 
                                end
                                scans_pre{i}=Vo;
                                Vi{1}=Vo; % Vi{2} is the p0                 
                                clear tmpout Vo Vi2
                        
                            end % end of skullstrip
                            
                            % brainmask and tissue specific noise estimation
                            % ==============================================
                            if estnoise==2
                                [~,f,e]=fileparts(scans_pre{i});
                                clear noisemask
                                for tc=1:3
                                    flags=struct('dmtx',0,'mask',0,'interp',1,'dtype',4);
                                    tmpout=spm_imcalc({Vi{2}},fullfile(outdir,['brain_mask_' TPM_Lab{tc} '_' num2str(i) '.nii']),['(i1>' num2str(tiss_th(tc)-0.5) ').*(i1<' num2str(tiss_th(tc)+0.5) ')'],flags); 
                                    noisemask{tc} = spm_data_read(tmpout.fname);
                                end
                                noisemask{4}=noisemask{1}+noisemask{2}+noisemask{3};
                                % warning: no rician mixture (JA default) used
                                Nii=nifti(scans_pre{i});
                                for tc=1:4
                                    f           = Nii.dat(:,:,:);
                                    noise(i,tc) = nanstd(f(noisemask{tc}>0.5));
                                end
                            end
                      
                      end  % end of loop over timepoints
                        
end  % end of intnorm skullstrip estnoise   
                 
% which noise estimate for longitudinal registration
% ========================================================
switch estnoise
                     
    case 0 % no noise estimatino use given sd 
        for i=1:nts
             fprintf('... ... scan %i fixed noise sd used is %.3f\n',i, sqrt(1/prec(i,:)) ) ;
        end
                            
    case 1 % updating the noise precision when skullstripping or just display, save
                             
        if skullstrip
             for i=1:nts % estimates of scanner noise using the whole scan
                  prec(i,1) = spm_noise_estimate(scans_pre{i})^(-2);
             end
             if sum(isnan(prec))>0
                 error('there was nan in noise estimate, please change estimatenoise or skullstrip parameters');
             end
             for i=1:nts
                  fprintf('... ... stripped scan %i noise sd used is %.2f\n',i, sqrt(1/prec(i,:)) ) ;
             end
        else
            for i=1:nts
                  fprintf('... ... scan %i noise sd used is %.2f\n',i, sqrt(1/prec(i,:)) ) ;
            end
        end
        
        
    case 2 % use tissue specific noise estimates
                            
        for tc=1:3
              for i=1:nts
                 disp(sprintf('... ... scan %i noise sd in %s is %.3f',i,TPM_Lab{tc},noise(i,tc)));
              end
        end
        for kk=1:size(noise,1)
              noise(kk,isnan(noise(kk,:)))=min(noise(kk,:)); % change that
        end
        prec = (noisesd*noise(:,2)).^(-2); % use noise estimate within white matter segment 
                            % note in case of tissue specific estimates the given value for sd are
                            % used as an additional factor for the estimates
                            % (set to 1 otherwise)
        disp('... ... using white matter noise sd');
        disp(sprintf('... ... scaled by %.4f',noisesd));
                 
end

try
   save(fullfile(outdir,'noise.mat'),'noise','prec');   
catch
   save(fullfile(outdir,'noise.mat'),'prec');
end

end