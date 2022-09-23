% Calculates median QSM values inside the Harvard ROIs and the number of voxels surviving after thresholding median QSM values inside the same ROIs, in the acquisition weighted space  
% written by Gisela Hagberg and edited by Elisa Tuzzi

clear all, close all, clc
maindir='./path/'; % folder with all subjects folders
qdir='./data/';  % folder with qsm and other input files (segmentations and rois) of each subject 
outdir=[maindir, 'LOAD']; % folder in which "AWIQSM_Harv.mat" file with paramaneters required for the plaque load calculation is saved  
sidx=[1 4 5 2 3]; % subject index
subj={'AD01' 'HC01' 'HC02' 'AD02' 'AD03'}; % example

subjj=subj(sidx); 
thr=[0.01:0.0009:0.2]; % qsm thresholds
gmthr=.9;  % GM threshold
load('HarvROis.mat')   % Harvard ROIs 


for s=1:size(subjj,2)
    subjdir=[maindir,num2str(subjj{s}),qdir];   % folder with qsm and other input files (segmentations and rois) of each subject
    sodir=[maindir,num2str(subjj{s}),'/LOAD/']; % output folder of each subject
    cd (subjdir)
    seg{1}=dir('c1*.nii');% % GM resliced to QSM space
    seg{2}=dir('c2*.nii');% % WM  resliced to QSM space 
    seg{3}=dir('c3*.nii');% % CSF resliced to QSM space
    AAL=dir('rwnpCNT1NEW_Reor_2_Eco3_2_AWI.nii') % AAL roi in AW space
    HOctx=dir('cort*.nii'); % cortical regions in AW-QSM space
    HOsc=dir('subcort*.nii'); % subcortical regions in AW-QSm space
    qsmfile=dir('qsm.nii');
    qsmname=qsmfile.name;
    
    Qmap=spm_read_vols(spm_vol([subjdir, qsmname])); % ppm
    Qmap(isnan(Qmap))=0;
    Qmap=1*Qmap*1000; % ppb
    
    
    mask=zeros(size(Qmap));
    for t=1:3
        tis{t}=spm_read_vols(spm_vol([subjdir, seg{t}.name]));
        mask=mask+tis{t};
    end
    
    GM=tis{1};
    nvox(s)=sum(GM(:)>0.90); 
    tiv3(s)=sum(mask(:)>0.8); 
    mask=mask>0.8;
    
    %get rid of cerebellum and brain stem
    %LBLs=spm_read_vols(spm_vol([roidir, SUIT.name]));
    %LBLs=1-(LBLs>0.1);
    
    %now gets rid of Ncaud put and thalamus
    LBL=spm_read_vols(spm_vol([subjdir, AAL.name]));   % figure, imagesc(LBL(:,:,8))
    subC=double((LBL>7000 & LBL<7105));
    
    spm_smooth(subC,subC,[5 5 5]);
    for t=1:3
        A=((tis{t}.*mask.*(subC<0.05).*(abs(Qmap)>0.00001))>0.9); %  qsm inside the GM  
        QSMval(s,t)=median(Qmap(A>0));
        QSMgmst(s,t)=std(Qmap(A>0));
        nvTis(s,t)=sum(A(:)>0);
    end
    
    A=(mask)>0.5;
    tiv3noCB(s)=sum(A(:));
    
    
    % --------------------
    % Harv rois:
    LBL2=spm_read_vols(spm_vol([subjdir, HOsc.name]));% figure, imagesc(LBL2(:,:,8))% sub
    LBL=spm_read_vols(spm_vol([subjdir, HOctx.name]));% figure, imagesc(LBL(:,:,8))% cort
   
    RH=LBL2>11;
    RH=100*RH;
    
    % Cortical rois:
    LBL1=LBL+RH;
    val(1:2:96)=[1:48];
    val(2:2:96)=[1:48]+100;
    
    gm=tis{1};
    
    
    for k=1:96
        A=(round(LBL)==val(k)).*mask.*(tis{1}>0.9); % qsm in the Harv cortical regions
        A=A.*(abs(Qmap)>0.00001);
        QS(s,k)=median(Qmap(A>0));
        QSst(s,k)=std(Qmap(A>0));
        
        % n voxels surviving in QSM cortical areas:
        nvoxROI(s,k)=sum(A(A>0));
    end
   %%
   
    %  Subcortical rois:
    k=97;
    val=[1 12 2 13 3 14 4 15 5 16 6 17 7 18 9 19 10 20 11 21 8]; % L R interleaved
    for nn=1:21
        A=(round(LBL2)==val(nn)).*mask;    
        A=A.*(abs(Qmap)>0.00001);
        QS(s,k)=median(Qmap(A>0));
        QSst(s,k)=std(Qmap(A>0));
        
        % n voxels surviving in QSM subcortical regions:
        nvoxROI(s,k)=sum(A(A>0)); 
        k=k+1;
    end
    

    toc
    fprintf(num2str(s));
    cd (sodir)
    save AWIQSM_Harv  Q* nv* tiv3* 

end
    cd (outdir)
    save AWIQSM_Harv  Q* nv* tiv3* 
    return
 