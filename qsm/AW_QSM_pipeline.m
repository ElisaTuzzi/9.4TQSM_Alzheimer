% Reconstructs magnitude and phase images as well as QSM from acquisition weighted (AW) high resolution .rms data acquired at the 9.4T 
% It requires Linux systems with FSL, some Matlab functions and Nifti Toolbox 
%%  

clear all; clc
addpath('./reco_9T') % folder with AW reco functions
addpath('./qsm/')  % AW-QSM pipeline
datadir=['./data/'] % raw data folder
maindir=['/path/'] % path with all subjects folders 
b0dir='./B0_files/' % folder with MTE dicom files for B0 reading

cc = tic;  
res = [0.1322 0.1322 0.61]; % you need to know the exact resolution in mm
subj={'subj1' 'subj2' 'subj3'};
for s=1:size(subj,2)
    fname=[b0dir,subj{s},'.IMA']
    hdr=dicominfo(fname);
    hdr.ProtocolName % check 
    %**********************
    % David Balla, Max Planck Institute for Biological Cybernetics %
    % get a dicom header (hdr) from a scan with exactly the same slice orientation
    B0 = cross(hdr.ImageOrientationPatient(1:3),hdr.ImageOrientationPatient(4:6))
    subjdir=[maindir,subj{s},'/']
    cd(subjdir)
    %use reco routine in "reco_9T" folder for reading the weighted acq data
    data = ReadSiToolData('/path/filename.rms');  % uses reco function "in reco_9T" folder. %Rolph Pohmann, Max Planck Institute for Biological Cybernetics%
    disp('Done ReadSiToolData\n');
    toc(cc);
    %%
    % Rolf Pohmann %
    [nx ny nz nc] = size(squeeze(data));
    d = shiftdim(single(squeeze(data)),3); %single instead of double for faster reco
    d = d(:,:,:,4:end-2); %crop in a sensibile way - outer slices for 3D acq (nz should be even), other dimensions to avoid overflooding the memory
    mtx = size(d); mtx = mtx(2:end); %adapt the matrix size 
    recon = adaptiveCombine3(d,[6 6 3]); % Inati 2014 Matlab/mex-function for phase sensitive combination of coil images
    Mag=abs(recon); % magnitude image
    Ph=angle(recon); % phase image
    nii = make_nii(Mag,res,[0 0 0],16); %Nifti Toolbox function - Matlab File Exchange
    save_nii(nii,[subjdir,'Mag.nii']); %Nifti Toolbox function
    nii = make_nii(Ph,res,[0 0 0],16); %Nifti Toolbox function - Matlab File Exchange
    save_nii(nii,[subjdir, 'Ph.nii']); %Nifti Toolbox function
    disp('Done adaptive combination of coil images\n');
    toc(cc);
    %%
    phname=[subjdir,'Ph.nii']
    ph=spm_read_vols(spm_vol(phname));
    mtx=size(ph)
    
    %***********************
    %%% David Balla %%%
    uw = unwrapLaplacian(ph,mtx,res); %MEDI Toolbox function (Cornell)
    %uw = unwrapLaplacian_floor(ph,mtx,res); % edited elisa
    %uw = LaplacianPhaseUnwrap(angle(recon),res,round(mtx.*[0.4 0.4 1.1])); %STISuite function (Duke)
    nii = make_nii(uw,res,[0 0 0],16); %Nifti Toolbox function
    save_nii(nii,[subjdir,'uw.nii']); %Nifti Toolbox function
    disp('Done Phase Unwrapping\n');
    toc(cc);
    %%
    
    %FSL functions
    !bet ./Mag.nii ./Mag -m -n -f 0.05 -Z 
    !gunzip *.gz
    %!fslchfiletype NIFTI /path/Mag_mask
    ref=spm_vol([maindir,subj{s},'/Mag.nii']);
    vol=spm_vol([maindir,subj{s},'/Mag_mask.nii']);
    volnew=vol;
    volnew.mat=ref.mat;
    %volnew.fname=[maindir,subj{k},Mag_',thr{j},'_mask_NEW.nii']
    mtrx=spm_read_vols(vol);
    spm_write_vol(volnew,mtrx);
    nii = load_untouch_nii([maindir,subj{s},'/Mag_mask.nii']); %Nifti Toolbox function
    mask = nii.img;
    mmask = single(imfill(mask,8,'holes'));
    nii = mmask;
    nii = make_nii(mmask,res,[0 0 0],16); %Nifti Toolbox function
    save_nii(nii,[maindir,subj{s},'/mmask.nii']); %Nifti Toolbox function
    disp('Done Masking\n');
    toc(cc);
    %%
    
    %Make a good mask - it is very important to avoid external signal sources
    %for a good reco
    maskfile=[subjdir,'mmask.nii']
    mask=spm_read_vols(spm_vol(maskfile));
    uw = uw - median(uw(mask==1)); %  demeaning   
    
    %Background field removal for a local phase map
    %input argument 4 is the kernel size and can be adapted. Optimal is around
    %6-8 times the inplane resolution
    [ph,mask_ero] = resharp(uw,mask,res,res(3)*2,1e-3); %Sun, H. and Wilman, A. H. (2013) MRM
    nii = make_nii(ph,res,[0 0 0],16);
    save_nii(nii,[subjdir,'ph_RESH.nii']);
    mask_ero = single(mask_ero);
    disp('Done Backgroud field removal\n');
    toc(cc);
    %%

    X = QSM_iLSQR(ph,mask_ero,'H',B0,'voxelsize',res,'padsize',round(mtx.*[0 0 0.4]),'niter',15,'TE',16.5,'B0',9.3875); %STISuite function (Duke)
    nii = make_nii(X,res,[0 0 0],16);
    save_nii(nii,[subjdir,'SUSC.nii']);
    %Xnew=permute(flipud(X),[2,1,3]);
    %nii = make_nii(Xnew,res,[0 0 0],16);
    %save_nii(nii,[subjdir,'SUSCnew.nii']);
    %vol=spm_vol([subjdir,'SUSCnew.nii']);    
    %volnew.fname=[subjdir,'qsm_test.nii']
    %v=spm_read_vols(vol);
    %p=fliplr(permute(v,[2,1,3]));
    %ref=spm_vol([subjdir,'AWI_flr_RES.nii']);
    %rf=spm_read_vols(ref);
    %refnew=ref;
    %volnew.mat=ref.mat;
    %refnew.fname=[subjdir1,'qsm_testDil20.nii']
    %spm_write_vol(refnew,p);
    disp('Done QSM\n');
    toc(cc);
end
return
