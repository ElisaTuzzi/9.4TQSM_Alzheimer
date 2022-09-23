% Initialisation (main directories and variables)
% Computes the plaque load fration from selected ROis and from the cortex of AD patients and
% controls in the AW space

clear all, close all, clc
maindir='./path/'; % directory with all subjects folders 
indir='./LOAD/'; % folder with .mat file (AWIQSM_Harv.mat) required for plaque load calculation in hte AW space 
outdir='./LOAD/'  % same folder as indir 
subj={'AD01' 'HC01' 'HC02' 'AD02' 'AD03'};
adidx=[1 4 5]; % AD indexes
hcidx=[2 3]; % HC indexes
adsubj=subj(adidx); 
hcsubj=subj(hcidx);
thr=[0.01:0.00015:0.05]; % QSM thresholds in ppm
gmthr=.9;  % GM threshold
gms={'09'}
%%

% input files
load ([indir, 'AWIQSM_Harv.mat']) % n voxels with median QSM > threshold, inside the Harvard ROIs in AW space 
load('./PL_lod/HarvROis.mat')   % Harvard ROIs 

% define ROIs 
AWIidx=find(all(nvoxROI~=0)); % find the indexes of the ROis where all subjcts have n voxels~=0; 24 indexes identified
AWIroivox=nvoxROI(:,AWIidx); % n voxels in the identified indexes 
DK(AWIidx);  % Harvard ROIs where all subjects have n voxels~=0

% Cortical:
cidx=AWIidx(1:16); % 16 cortical ROIs identified
cort=DK(cidx) % 

%%******
% Harvard Cortical ROIs indexes:
val(1:2:96)=[1:48];
%val(2:2:96)=[1:48];
val(2:2:96)=[1:48]+100;
cortval=val(cidx); 
rois={'GM' 'Cort'}
%%


%************************
%plaque load calculation:
                                    
% ad:      
adfracvenroi=zeros(4,length(adsubj));
for j=1:length(adsubj)
    addir=[maindir,num2str(adsubj{j}),indir]; % input data file folder
    odir=[maindir,num2str(adsubj{j}),'/ANALYSIS/LOAD/'];
    cd (addir)
    j
    
    %qsm:
    adqfile=dir('qsm.nii');
    adqname=adqfile.name;
    adqvol=spm_vol(adqname);
    adq=spm_read_vols(adqvol);
    adq(isnan(adq))=0;
    mtx=size(adq);    

    % set gm & rois:
    adgmfile=dir('c1*.nii'); % GM in AWI space
    adcortfile=dir('rwnp1*.nii');  % cort in awi space
    %adsubfile=dir('rwnp2*.nii'); % sub in awi space
    adrois={adgmfile adcortfile}; % gm, cort, sub
   
    % mask:
    seg{1}=dir('c1*.nii'); % GM  registered to QSM space
    seg{2}=dir('c2*.nii'); % WM  registered to QSM space 
    seg{3}=dir('c3*.nii'); % CSF registered to QSM space
    
    mask=zeros(size(adq));
    for t=1:3
        tis{t}=spm_read_vols(spm_vol([addir, seg{t}.name])); 
        tis{t}=tis{t};
        mask=mask+tis{t};
    end
    mask=mask>0.8;
    gm=tis{1};

    % define veins
    ven=adq>0.05; % cut over 50 ppb
    
    
    
    %*********************************
    % plaque load from GM: 
    figure(j)
    for rs=1:size(adrois,2)
        adroiname=adrois{rs}.name;
        adroivol=spm_vol(adroiname);
        adroi=spm_read_vols(adroivol); % gm
        rs
        
        if rs==1   %% GM
            adroi=(1-ven).*((adroi.*mask.*(abs(adq)>0.00001))>0.9); % thresholded qsm. No vein contribution 
            adroi=adroi>0;  
            n=find(adroi==1); size(n);
            npixMASK=size(n,1) % n pixels inside the masked qsm
            
            % vein frac:
            adnv=size(find(ven==1),1) % n pixels with susceptibility > 50 ppb
            adfracven=adnv./size(find(adroi==1),1) % n fraction of vein pixels inside the masked qsm
            adfracvenaroi(rs, j)=adfracven
            
            adroimasked=adq.*adroi;  
            adroivolnew=adroivol; 
            adroivolnew.fname=[addir, num2str(rois{rs}),'09_admasked.nii'];
            spm_write_vol(adroivolnew,adroimasked); 
            adroiMasked=adroimasked.*(adroimasked>0); 
            adroivolnew.fname=[addir, num2str(rois{rs}), '09_adMasked.nii'];
            spm_write_vol(adroivolnew,adroiMasked); 
            adroivolnew.fname=[addir, num2str(rois{rs}),'_09.nii'];
            spm_write_vol(adroivolnew,adroi); 
            
                        
            %pl frac from gm:
            adf=zeros(size(thr));
            for t=1:length(thr)
                th=adroiMasked>thr(t); 
                n=find(th==1);
                npix=size(n,1);
                f=npix/npixMASK;
                adf(1,t)=f;
            end
           % figure(j)
            fad_mtx_gm(j,:)=adf;  
            filename = [odir,num2str(rois{rs}),'_adf.mat']
            p = fad_mtx_gm;
            save(filename, 'p');
            cd (addir)
           
            
        %******************** 
        % Cort ROIs:
        else    %% Cort
            for rr=1:size(cortval,2)  % indexes of cortical ROIs common to all subjects in AW space
                rr
                cval=(round(adroi)==cortval(rr)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(adq)>0.00001);
                commROI=(cval);
                commROI=commROI>0; 
                n=find(commROI==1);
                npixMASK=size(n,1)
                
                %masked qsm:
                adroimasked=adq.*commROI;
                adroivolnew=adroivol;
                adroivolnew.fname=[addir, num2str(cort{rr}), num2str(gms{1}),'_adqmasked.nii'];
                spm_write_vol(adroivolnew,adroimasked); 
                adroiMasked=adroimasked.*(adroimasked>0);
                adroivolnew=adroivol;
                adroivolnew.fname=[addir, num2str(cort{rr}), num2str(gms{1}),'_adroiMasked.nii'];
                spm_write_vol(adroivolnew,adroiMasked); 
                adroivolnew.fname=[addir, num2str(cort{rr}), num2str(gms{1}),'.nii'];
                spm_write_vol(adroivolnew,commROI); 
                
                
                %plaque fraction from common cortical rois:
                adf=zeros(size(thr));
                for t=1:length(thr)
                    th=adroiMasked>thr(t); 
                    n=find(th==1);
                    npix=size(n,1);
                    f=npix/npixMASK;
                    adf(1,t)=f;
                end
                comm_fad_cort(rr,:)=adf;  
            end
            close all
            filename = [odir,'CortcommROIs_adf.mat']
            p2 = comm_fad_cort;
            save(filename, 'p2');    
            cd (addir)            
            comm_fad_mtx_cort(j,:,:)=comm_fad_cort;
            
            
            %*********************
            % Cortex:
            adroi=(round(adroi)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(adq)>0.00001);
            adroi=adroi>0;  
            n=find(adroi==1);
            npixMASK=size(n,1)
            
            % vein fraction:
            adnv=size(find(ven==1),1)
            adfracven=adnv./size(find(adroi==1),1)
            adfracvenaroi(rs, j)=adfracven
            
            %masked qsm:
            adroimasked=adq.*adroi;    
            adroiMasked=adroimasked.*(adroimasked>0);
            adroivolnew=adroivol;
            adroivolnew.fname=[addir, num2str(rois{rs}), num2str(gms{1}),'_adroiMasked.nii'];
            spm_write_vol(adroivolnew,adroiMasked); 
            adroivolnew.fname=[addir, num2str(rois{rs}), num2str(gms{1}),'.nii'];
            spm_write_vol(adroivolnew,adroi); 
            
                        
            %plaque fraction from the cortex:
            adf=zeros(size(thr));
            for t=1:length(thr)
                th=adroiMasked>thr(t); 
                n=find(th==1);
                npix=size(n,1);
                f=npix/npixMASK;
                adf(1,t)=f;
            end
            fad_mtx_cort(j,:)=adf; 
            filename = [odir,num2str(rois{rs}),'_adf.mat']
            p3 = fad_mtx_cort;
            save(filename, 'p3');
        end
        cd(addir)
    end
end
close all
cd(outdir)
save ADPL_AWsp fad_mtx_gm comm_fad_mtx_cort fad_mtx_cort thr adfracvenaroi

%%


% HC:          
hcfracvenroi=zeros(4,length(hcsubj));
for k=1:length(hcsubj)
    hcdir=[maindir,num2str(hcsubj{k}),indir];
    odir=[maindir,num2str(hcsubj{k}),'/ANALYSIS/LOAD/'];
    cd (hcdir)
    k
    
    %qsm:
    hcqfile=dir('qsm_flr_Reor.nii');
    hcqname=hcqfile.name;
    hcqvol=spm_vol(hcqname);
    hcq=spm_read_vols(hcqvol);
    hcq(isnan(hcq))=0;
    
    % vein fraction:
    ven=hcq>0.05; % vein contribution; cuts over 50 ppb 
    mtx=size(hcq);  
    
    % set gm & rois:
    hcgmfile=dir('c1*.nii'); % GM in AWI space
    hccortfile=dir('rwnp1*.nii');  % cort in awi space
    %hcsubfile=dir('rwnp2*.nii'); % sub in awi space
    hcrois={hcgmfile hccortfile}; % gm, aal, cort, sub
    
    % mask:
    seg{1}=dir('c1*.nii');% % GM  registered to QSM space
    seg{2}=dir('c2*.nii');% % WM  registered to QSM space 
    seg{3}=dir('c3*.nii');% % CSF registered to QSM space
    
    mask=zeros(size(hcq));
    for t=1:3
        tis{t}=spm_read_vols(spm_vol([hcdir, seg{t}.name])); 
        mask=mask+tis{t};
    end
    mask=mask>0.8;
    gm=tis{1};
    
    
    %*********************************
    % QSM from GM: (th same as above)
    for rs=1:size(hcrois,2)
        hcroiname=hcrois{rs}.name;
        hcroivol=spm_vol(hcroiname);
        hcroi=spm_read_vols(hcroivol); % gm
        
        rs

        %thresholding:
        if rs==1   %% GM
            hcroi=(1-ven).*((hcroi.*mask.*(abs(hcq)>0.00001))>0.9);  %always the same voxels  
            hcroi=hcroi>0;  
            n=find(hcroi==1);
            npixMASK=size(n,1)
   
            % vein fraction:
            hcnv=size(find(ven==1),1)
            hcfracven=hcnv./size(find(hcroi==1),1)
            hcfracvenaroi(rs, k)=hcfracven
            
            % masked qsm
            hcroimasked=hcq.*hcroi;  
            hcroivolnew=hcroivol; 
            hcroivolnew.fname=[hcdir, num2str(rois{rs}),'09_hcmasked.nii'];
            spm_write_vol(hcroivolnew,hcroimasked); 
            hcroiMasked=hcroimasked.*(hcroimasked>0); 
            hcroivolnew.fname=[hcdir, num2str(rois{rs}),'09_hcMasked.nii'];
            spm_write_vol(hcroivolnew,hcroiMasked); 
            hcroivolnew.fname=[hcdir, num2str(rois{rs}),'_09.nii'];
            spm_write_vol(hcroivolnew,hcroi); 
            
            
            %plaque fraction from gm:
            hcf=zeros(size(thr));
            for t=1:length(thr)
                th=hcroiMasked>thr(t); 
                n=find(th==1);
                npix=size(n,1);
                f=npix/npixMASK;
                hcf(1,t)=f;
            end
            fhc_mtx_gm(k,:)=hcf; 
            filename = [odir,num2str(rois{rs}),'_hcf.mat']
            p = fhc_mtx_gm;
            save(filename, 'p');
            cd (hcdir)
            %%
            
            
           %*******************
           % Cort ROIs:
        else   %% Cort
            for rr=1:size(cortval,2)
                rr
                cval=(round(hcroi)==cortval(rr)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(hcq)>0.00001);
                commROI=(cval); 
                commROI=commROI>0;
                n=find(commROI==1);
                npixMASK=size(n,1)

                %masked qsm:
                hcroimasked=hcq.*commROI;
                hcroivolnew=hcroivol;
                hcroivolnew.fname=[hcdir, num2str(cort{rr}), num2str(gms{1}),'_adqmasked.nii'];
                spm_write_vol(hcroivolnew,hcroimasked); 
                hcroiMasked=hcroimasked.*(hcroimasked>0);
                hcroivolnew=hcroivol;
                hcroivolnew.fname=[hcdir, num2str(cort{rr}), num2str(gms{1}),'_adqmasked.nii'];
                spm_write_vol(hcroivolnew,hcroiMasked); 
                hcroivolnew.fname=[hcdir, num2str(cort{rr}), num2str(gms{1}),'.nii'];
                spm_write_vol(hcroivolnew,commROI); 
                
                
                %plaque fraction from common rois:
                hcf=zeros(size(thr));
                for t=1:length(thr)
                    th=hcroiMasked>thr(t); 
                    n=find(th==1);
                    npix=size(n,1);
                    f=npix/npixMASK;
                    hcf(1,t)=f;
                end
                comm_fhc_cort(rr,:)=hcf;  
            end
            close all
            filename = [odir,'CortcommROIs_hcf.mat']
            p2 = comm_fhc_cort;
            save(filename, 'p2');    
            cd (hcdir)            
            comm_fhc_mtx_cort(k,:,:)=comm_fhc_cort;
            %%
            
            %*************
            % Cortex:
            hcroi=(round(hcroi)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(hcq)>0.00001);
            hcroi=hcroi>0;  
            n=find(hcroi==1);
            npixMASK=size(n,1)
            
            % vein frac:
            hcnv=size(find(ven==1),1)
            hcfracven=hcnv./size(find(hcroi==1),1)
            hcfracvenaroi(rs, k)=hcfracven

            %masked qsm:
            hcroimasked=hcq.*hcroi;  
            hcroiMasked=hcroimasked.*(hcroimasked>0);
            hcroivolnew=hcroivol;
            hcroivolnew.fname=[hcdir, num2str(rois{rs}), num2str(gms{1}),'_hcroiMasked.nii'];
            spm_write_vol(hcroivolnew,hcroiMasked); 
            hcroivolnew.fname=[hcdir, num2str(rois{rs}), num2str(gms{1}),'.nii'];
            spm_write_vol(hcroivolnew,hcroi); 
            mask=mask>0.8;
            
            %plaque fraction from the cortex:
            hcf=zeros(size(thr));
            for t=1:length(thr)
                th=hcroiMasked>thr(t); 
                n=find(th==1);
                npix=size(n,1);
                f=npix/npixMASK;
                hcf(1,t)=f;
            end
            fhc_mtx_cort(k,:)=hcf; 
            filename = [odir,num2str(rois{rs}),'_hcf.mat']
            p3 = fhc_mtx_cort;
            save(filename, 'p3');
        end
        cd(hcdir)
    end
end
%close all
cd(outdir)
save HCPL_AWsp fhc_mtx_gm  comm_fhc_mtx_cort fhc_mtx_cort thr hcfracvenaroi

return