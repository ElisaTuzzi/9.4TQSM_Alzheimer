% Initialisation (main directories and variables)
% Computes the plaque load fration from selected ROis and from the cortex of AD patients and
% controls in the MTE space

clear all, close all, clc
maindir='./path/'; % directory with all subjects folders 
indir='./LOAD/';  % folder with .mat file (MTEQSM_Harv.mat) required for plaque load calculation in hte MTE space  
outdir='./LOAD/'  % same folder as indir  
subj={'AD01' 'HC01' 'HC02' 'AD02' 'AD03'};
adidx=[1 4 5]; % AD indexes
hcidx=[2 3]; % HC indexes
adsubj=subj(adidx); 
hcsubj=subj(hcidx);
thr=[0.01:0.00015:0.05]; % QSM thresholds in ppm
gmthr=.9;  % GM threshold
gms={'09'}
cat={'CoR', 'NoC'}; % NAV correction; no NAV correction
sus={'PARAM' 'DIAM'}; 
tk={'12', '3'}; 
ks={'12', '16'}; 

load('./PL_load/HarvROis.mat')   % Harvard ROIs 


%**************************
% plaque load calculation:

%ad:      
adfracvenroi=zeros(4,length(adsubj));
for j=1:length(adsubj)
    addir=[maindir,num2str(adsubj{j}),indir];
    odir=[maindir,num2str(adsubj{j}),'/ANALYSIS/LOAD/'];
    cd (addir)
    j
    
    %qsm:
    adqfile=dir('qsmTke-*_FLIP.nii')
    qfname=adqfile(1,1).name;
    str1=[findstr(qfname,'_')];
    str2=[findstr(qfname,'-')];
    strfile=[qfname(str1(2)+19:str2(5)+3)]
    for ct=1:size(cat,2)
        ct
        for ss=1 % param
            ss
            for tt=1:size(tk,2)
                tt
                for kk=1:size(ks,2)
                    kk
                    
                    %QSM:
                    adqname=['qsmTke-',num2str(tk{tt}), '_ks', num2str(ks{kk}),'_', num2str(cat{ct}),'s', num2str(adsubj{j}), strfile, num2str(sus{ss}),'_FLIP.nii']
                    adqvol=spm_vol(adqname);
                    adq=spm_read_vols(adqvol);
                    adq(isnan(adq))=0;
                    mtx=size(adq);    

                    % set gm & rois:
                    adgmfile=dir('*2_ECO3_c1*.nii');% GM in MTE space
                    adcortfile=dir('*2_ECO3_*cort*Reor.nii'); %  cort in MTE space
                    adrois={adgmfile adcortfile}; % gm, cort

                    % mask:
                    seg{1}=dir('*2_ECO3_c1*.nii');% % GM in QSM space
                    seg{2}=dir('*2_ECO3_c2*.nii');% % WM  in QSM space 
                    seg{3}=dir('*2_ECO3_c3*.nii');% % CSF in QSM space
                    
                    mask=zeros(size(adq));
                    for t=1:3
                        tis{t}=spm_read_vols(spm_vol([addir, seg{t}.name])); 
                        tis{t}=tis{t};
                        mask=mask+tis{t};
                    end
                    mask=mask>0.8;
                    gm=tis{1};
                    
                    % vein definition
                    ven=adq>0.05; % cut over 50 ppb

                   %***************************************
                    %Harvard common cortical ROIs definition:
                    
                    % comm Harv:
                    nvroi=[indir,'/MTEQSM_Harv_', num2str(cat{ct}),'_Tk', num2str(tk{tt}),'_ks', num2str(ks{kk}), '_', num2str(sus{ss}),'.mat'] 

                    % index definition:
                    MTEidx=find(all(nvoxROI~=0)); % find indexes of Harvard ROis in which all subj have voxels
                    MTEroivox=nvoxROI(:,MTEidx); 
                    DK(MTEidx);

                    % cort:
                    cidx=MTEidx(1:78) % cortical indexes of Harvard ROIs identified
                    cort=DK(cidx)  
                    
                    % indexes of bilateral cortical ROIs:
                    val(1:2:96)=[1:48]; 
                    val(2:2:96)=[1:48]+100;
                    cortval=val(cidx); 
                 
                    
                    %*********************************
                    % QSM from GM: 
                    figure(j)
                    for rs=1:size(adrois,2)
                        adroiname=adrois{rs}.name;
                        adroivol=spm_vol(adroiname);
                        adroi=spm_read_vols(adroivol); % gm
                        rs
                        
                        %thresholding:
                        if rs==1   %% GM
                            adroi=(1-ven).*((adroi.*mask.*(abs(adq)>eps))>0.9);   
                            n=find(adroi==1); size(n);
                            npixMASK=size(n,1)

                            % vein fraction:
                            adnv=size(find(ven==1),1)
                            adfracven=adnv./size(find(adroi==1),1)
                            adfracvenaroi(rs, j)=adfracven

                            % masked qsm
                            adroimasked=adq.*adroi; 
                            adroivolnew=adroivol; 
                            adroivolnew.fname=[addir, num2str(rois{rs}),'09_admasked.nii'];
                            spm_write_vol(adroivolnew,adroimasked); 
                            adroiMasked=adroimasked.*(adroimasked>0);  
                            adroivolnew.fname=[addir, num2str(rois{rs}), '09_adMasked.nii'];
                            spm_write_vol(adroivolnew,adroiMasked); 
                            adroivolnew.fname=[addir, num2str(rois{rs}),'_09.nii'];
                            spm_write_vol(adroivolnew,adroi); 
                            %***

                            %plaque fraction from gm:
                            adf=zeros(size(thr));
                            for t=1:length(thr)
                                th=adroiMasked>thr(t); 
                                n=find(th==1);
                                npix=size(n,1);
                                f=npix/npixMASK;
                                adf(1,t)=f;
                            end
                            fad_mtx_gm(j,:)=adf;  
                            filename = [odir,num2str(rois{rs}),'_adf.mat']
                            p = fad_mtx_gm;
                            save(filename, 'p');
                        

                        %********************************
                        % plaque load from Cortical ROIs:
                        
                        else   %% Cort
                            for rr=1:size(cortval,2)
                                rr
                                cval=(round(adroi)==cortval(rr)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(adq)>0.00001);
                                commROI=(cval); 
                                commROI=commROI>0; %figure, for g=1:96, imagesc(commROI(:,:,g)),colorbar, g=gca; title(num2str(cort{(rr)})), pause(0.01), end
                                n=find(commROI==1);
                                npixMASK=size(n,1)

                                % masked qsm
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
                            
                                % plaque fracion calculation:
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


                           %**********************
                            % Cortex:
                            adroi=(round(adroi)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(adq)>0.00001);
                            adroi=adroi>0;  
                            n=find(adroi==1);
                            npixMASK=size(n,1)
                         
                            
                            % vein frac:
                            adnv=size(find(ven==1),1)
                            adfracven=adnv./size(find(adroi==1),1)
                            adfracvenaroi(rs, j)=adfracven

                            % masked qsm:
                            adroimasked=adq.*adroi;    
                            adroiMasked=adroimasked.*(adroimasked>0);
                            adroivolnew=adroivol;
                            adroivolnew.fname=[addir, num2str(rois{rs}), num2str(gms{1}),'_adroiMasked.nii'];
                            spm_write_vol(adroivolnew,adroiMasked); 
                            adroivolnew.fname=[addir, num2str(rois{rs}), num2str(gms{1}),'.nii'];
                            spm_write_vol(adroivolnew,adroi); 


                            %pl frac from cort rois:
                            adf=zeros(size(thr));
                            for t=1:length(thr)
                                th=adroiMasked>thr(t); 
                                n=find(th==1);
                                npix=size(n,1);
                                f=npix/npixMASK;
                                adf(1,t)=f;
                            end
                            fad_mtx_cort(j,:)=adf; %figure,
                            filename = [odir,num2str(rois{rs}),'_adf.mat']
                            p3 = fad_mtx_cort;
                            save(filename, 'p3');
                        end
                        cd (addir)
                    end
                    %cd(loaddir), 
                    flname=[outdir,'/ADPL_MTEsp_Tk',num2str(tk{tt}),'_', num2str(cat{ct}), '_ks', num2str(ks{kk}), '_', num2str(sus{ss}),'.mat']
                    save(flname,'fad_mtx_gm','comm_fad_mtx_cort','fad_mtx_cort','thr','adfracvenaroi')
                    cd (addir) 
                end
            end
        end
    end
end
close all
        

% HC:          
hcfracvenroi=zeros(4,length(hcsubj));
for k=1:length(hcsubj)
    hcdir=[maindir,num2str(hcsubj{k}),indir];
    odir=[maindir,num2str(hcsubj{k}),'/ANALYSIS/LOAD/'];
    cd (hcdir)
    k
    
    %qsm:
    hcqfile=dir('qsmTke-*_FLIP.nii')
    qfname=hcqfile(1,1).name
    str1=[findstr(qfname,'_')];
    str2=[findstr(qfname,'-')];
    strfile=[qfname(str1(2)+19:str2(5)+3)]
    for ct=1:size(cat,2)
        ct
        for ss=1
            ss
            for tt=1:size(tk,2)
                tt
                for kk=1:size(ks,2)
                    kk
                    
                    %QSM:
                    hcqname=['qsmTke-',num2str(tk{tt}), '_ks', num2str(ks{kk}),'_',num2str(cat{ct}),'s', num2str(hcsubj{k}), strfile, num2str(sus{ss}),'_FLIP.nii']
                    hcqvol=spm_vol(hcqname);
                    hcq=spm_read_vols(hcqvol);
                    hcq(isnan(hcq))=0;
                    mtx=size(hcq);    

                    % set gm & rois:
                    hcgmfile=dir('*2_ECO3_c1*.nii');% GM in MTE space
                    hccortfile=dir('*2_ECO3_*cort*Reor.nii'); %  cort in MTE space
                    hcrois={hcgmfile hccortfile}; % gm, cort

                    % mask:
                    seg{1}=dir('*2_ECO3_c1*.nii');% % GM in QSM space
                    seg{2}=dir('*2_ECO3_c2*.nii');% % WM  in QSM space 
                    seg{3}=dir('*2_ECO3_c3*.nii');% % CSF in QSM space
                    
                    mask=zeros(size(hcq));
                    for t=1:3
                        tis{t}=spm_read_vols(spm_vol([hcdir, seg{t}.name])); 
                        tis{t}=tis{t};
                        mask=mask+tis{t};
                    end
                    mask=mask>0.8;

                    ven=hcq>0.05; % cut over 50 ppb
                    gm=tis{1};


                   %***************************************
                    %Harvard common cortical ROIs definition:
                    
                    % comm Harv:
                    nvroi=[indir,'/MTEQSM_Harv_', num2str(cat{ct}),'_Tk', num2str(tk{tt}),'_ks', num2str(ks{kk}), '_', num2str(sus{ss}),'.mat'] 

                    % index definition:
                    MTEidx=find(all(nvoxROI~=0)); % find indexes of Harvard ROis in which all subj have voxels
                    MTEroivox=nvoxROI(:,MTEidx); 
                    DK(MTEidx);

                    % cort:
                    cidx=MTEidx(1:78) % cortical indexes of Harvard ROIs identified
                    cort=DK(cidx)  
                    
                    % indexes of bilateral cortical ROIs:
                    val(1:2:96)=[1:48]; 
                    val(2:2:96)=[1:48]+100;
                    cortval=val(cidx); 
                 
                    
                    %*********************************
                    % QSM from GM: (th same as above)
                    
                    figure(j)
                    for rs=1:size(hcrois,2)
                        hcroiname=hcrois{rs}.name;
                        hcroivol=spm_vol(hcroiname);
                        hcroi=spm_read_vols(hcroivol); % gm
                        rs
                        
                        %thresholding:
                        if rs==1   %% GM
                            hcroi=(1-ven).*((hcroi.*mask.*(abs(hcq)>eps))>0.9);   
                            n=find(hcroi==1); size(n);
                            npixMASK=size(n,1)

                            % vein fraction:
                            hcnv=size(find(ven==1),1)
                            hcfracven=hcnv./size(find(hcroi==1),1)
                            hcfracvenaroi(rs, j)=hcfracven

                            % masked qsm
                            hcroimasked=hcq.*hcroi; 
                            hcroivolnew=hcroivol; 
                            hcroivolnew.fname=[hcdir, num2str(rois{rs}),'09_hcmasked.nii'];
                            spm_write_vol(hcroivolnew,hcroimasked); 
                            hcroiMasked=hcroimasked.*(hcroimasked>0);  
                            hcroivolnew.fname=[hcdir, num2str(rois{rs}), '09_hcMasked.nii'];
                            spm_write_vol(hcroivolnew,hcroiMasked); 
                            hcroivolnew.fname=[hcdir, num2str(rois{rs}),'_09.nii'];
                            spm_write_vol(hcroivolnew,hcroi); 
                            %***

                            %plaque fraction from gm:
                            hcf=zeros(size(thr));
                            for t=1:length(thr)
                                th=hcroiMasked>thr(t); 
                                n=find(th==1);
                                npix=size(n,1);
                                f=npix/npixMASK;
                                hcf(1,t)=f;
                            end
                            fhc_mtx_gm(j,:)=hcf;  
                            filename = [odir,num2str(rois{rs}),'_hcf.mat']
                            p = fhc_mtx_gm;
                            save(filename, 'p');
                        

                        %********************************
                        % plaque load from Cortical ROIs:
                        
                        else   %% Cortical ROIs
                            for rr=1:size(cortval,2)
                                rr
                                cval=(round(hcroi)==cortval(rr)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(hcq)>0.00001);
                                commROI=(cval); 
                                commROI=commROI>0; figure, for g=1:96, imagesc(commROI(:,:,g)),colorbar, g=gca; title(num2str(cort{(rr)})), pause(0.01), end
                                n=find(commROI==1);
                                npixMASK=size(n,1)

                                % masked qsm
                                hcroimasked=hcq.*commROI;
                                hcroivolnew=hcroivol;
                                hcroivolnew.fname=[hcdir, num2str(cort{rr}), num2str(gms{1}),'_hcqmasked.nii'];
                                spm_write_vol(hcroivolnew,hcroimasked); 
                                hcroiMasked=hcroimasked.*(hcroimasked>0);
                                hcroivolnew=hcroivol;
                                hcroivolnew.fname=[hcdir, num2str(cort{rr}), num2str(gms{1}),'_hcroiMasked.nii'];
                                spm_write_vol(hcroivolnew,hcroiMasked); 
                                hcroivolnew.fname=[hcdir, num2str(cort{rr}), num2str(gms{1}),'.nii'];
                                spm_write_vol(hcroivolnew,commROI); 
                            
                                % plaque fracion calculation:
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


                           %**********************
                            % from the Cortex:
                            hcroi=(round(hcroi)).*mask.*(tis{1}>gmthr).*(1-ven).*(abs(hcq)>0.00001);
                            hcroi=hcroi>0;  
                            n=find(hcroi==1);
                            npixMASK=size(n,1)
                         
                            
                            % vein frac:
                            hcnv=size(find(ven==1),1)
                            hcfracven=hcnv./size(find(hcroi==1),1)
                            hcfracvenaroi(rs, j)=hcfracven

                            % masked qsm:
                            hcroimasked=hcq.*hcroi;    
                            hcroiMasked=hcroimasked.*(hcroimasked>0);
                            hcroivolnew=hcroivol;
                            hcroivolnew.fname=[hcdir, num2str(rois{rs}), num2str(gms{1}),'_hcroiMasked.nii'];
                            spm_write_vol(hcroivolnew,hcroiMasked); 
                            hcroivolnew.fname=[hcdir, num2str(rois{rs}), num2str(gms{1}),'.nii'];
                            spm_write_vol(hcroivolnew,hcroi); 


                            %pl frac from the cortex:
                            hcf=zeros(size(thr));
                            for t=1:length(thr)
                                th=hcroiMasked>thr(t); 
                                n=find(th==1);
                                npix=size(n,1);
                                f=npix/npixMASK;
                                hcf(1,t)=f;
                            end
                            fhc_mtx_cort(j,:)=hcf; %figure,
                            filename = [odir,num2str(rois{rs}),'_hcf.mat']
                            p3 = fhc_mtx_cort;
                            save(filename, 'p3');
                        end
                        cd(hcdir)
                    end
                    %cd(outdir)
                    flname=[outdir,'/HCPL_MTEsp_Tk',num2str(tk{tt}),'_', num2str(cat{ct}), '_ks', num2str(ks{kk}),'_', num2str(sus{ss}),'.mat']
                    save(flname,'fhc_mtx_gm','comm_fhc_mtx_cort','fhc_mtx_cort','thr','hcfracvenaroi')
                    cd (hcdir) 
                    end
                end
            end
        end
    end
close all
return

