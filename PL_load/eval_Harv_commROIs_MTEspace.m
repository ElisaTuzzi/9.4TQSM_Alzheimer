% Calculates median QSM values inside the Harvard ROIs and the number of voxels surviving after thresholding median QSM values inside the same ROIs, in the MTE space  
% written by Gisela Hagberg and edited by Elisa Tuzzi

close all, clear all, clc
maindir='./path/'; % directory with all subjects folders 
indir='./data/'; % folder with all input files (qsm, rois, segmentations)
outdir=[maindir, 'LOAD']; % folder in which "MTEQSM_Harv.mat" file with paramaneters required for the plaque load calculation is saved 

subj={'AD01' 'HC01' 'HC02' 'AD02' 'AD03'};
sidx=[1 4 5 2 3]; % subject index

cat={'CoR', 'NoC'}; % NAV correction; no NAV correction
sus={'PARAM'}; 
tk={'12', '3'}; % Tikhonov param
ks={'12', '16'};  % kernel size

subjj=subj(sidx);
gmthr=.9; % GM threshold
load('HarvROis.mat')   % Harvard ROIs 


for s=1:size(subjj,2)
    subjdir=[maindir,num2str(subjj{s}),indir];   
    %sodir=[maindir,num2str(subjj{s}),'/LOAD/'];
    cd (subjdir)
    seg{1}=dir('*2_ECO3_c1*.nii');% % GM in MTE-QSM space
    seg{2}=dir('*2_ECO3_c2*.nii');% % WM  in MTE-QSM space 
    seg{3}=dir('*2_ECO3_c3*.nii');% % CSF in MTE-QSM space
    AAL=dir('*2_ECO3_*AAL_*Reor.nii');  % AAL roi in MTe space
    %SUIT=dir('rwnpCNTs55SUIT*.nii'); % 
    %JHU=dir('rwnpCNTs55JHU*.nii'); % 
    HOctx=dir('*2_ECO3_*cort*Reor.nii'); % cortical regions in MTE-QSM space
    HOsc=dir('*2_ECO3_*sub*Reor.nii') ;% subcortical regions in MTE-QSM space

    qsmfile=dir('qsmTke-*_FLIP.nii'); % qsm 
    qfname=qsmfile(1,1).name;
    str1=[findstr(qfname,'_')];
    str2=[findstr(qfname,'-')];
    strfile=[qfname(str1(2)+19:str2(5)+3)];
    for ct=1:size(cat,2)
        for tt=1:size(tk,2)
            for ss=1:size(sus,2)
                for kk=1:size(ks,2)
                qsmname=['qsmTke-',num2str(tk{tt}), '_ks', num2str(ks{kk}),'_',num2str(cat{ct}),'s', num2str(subjj{s}), strfile, num2str(sus{ss}),'_FLIP.nii']

                Qmap=spm_read_vols(spm_vol([subjdir, qsmname])); % ppm
                Qmap(isnan(Qmap))=0;
                Qmap=1*Qmap*1000; % ppb

                mask=zeros(size(Qmap));
                for t=1:3
                    tis{t}=spm_read_vols(spm_vol([subjdir, seg{t}.name]));
                    mask=mask+tis{t};
                end
                GM=tis{1};
                nvox(s)=sum(GM(:)>0.98); 
                tiv3(s)=sum(mask(:)>0.5);
                mask=mask>0.8;

                %get rid of cerebellum and brain stem
                %LBLs=spm_read_vols(spm_vol([roidir, SUIT.name]));
                %LBLs=1-(LBLs>0.1);

                %now get rid of Ncaud put and thalamus
                LBL=spm_read_vols(spm_vol([subjdir, AAL.name]));   % figure, imagesc(LBL(:,:,8))
                subC=double((LBL>7000 & LBL<7105));
                spm_smooth(subC,subC,[5 5 5]);
                for t=1:3
                    A=((tis{t}.*mask.*(subC<0.05).*(abs(Qmap)>0.00001))>0.9);
                    QSMval(s,t)=median(Qmap(A>0));
                    QSMgmst(s,t)=std(Qmap(A>0));
                    nvTis(s,t)=sum(A(:)>0);
                end

                A=(mask)>0.5;
                tiv3noCB(s)=sum(A(:));
                %%
                %GlPall=spm_read_vols(spm_vol([addir, GlP.name]));
                %Put=spm_read_vols(spm_vol([addir, PUT.name]));

                % --------------------
                % Harv rois:
                LBL2=spm_read_vols(spm_vol([subjdir, HOsc.name]));% figure, imagesc(LBL2(:,:,8))% sub
                LBL=spm_read_vols(spm_vol([subjdir, HOctx.name]));% figure, imagesc(LBL(:,:,8))% cort


                RH=LBL2>11;
                RH=100*RH;

                % Cort:
                LBL=LBL+RH;
                val(1:2:96)=[1:48];
                val(2:2:96)=[1:48]+100;

                gm=tis{1};
                for k=1:96
                    A=(round(LBL)==val(k)).*mask.*(tis{1}>0.9);  
                    A=A.*(abs(Qmap)>eps); 
                    QS(s,k)=median(Qmap(A>0));
                    QSst(s,k)=std(Qmap(A>0));
                    
                    % n voxels surviving in QSM cort areas:
                    nvoxROI(s,k)=sum(A(A>0));
                end

                %  Subcort:
                k=97;
                val=[1 12 2 13 3 14 4 15 5 16 6 17 7 18 9 19 10 20 11 21 8];%L R interleaved
                for nn=1:21
                    A=(round(LBL2)==val(nn)).*mask;     
                    A=A.*(abs(Qmap)>eps); 

                    QS(s,k)=median(Qmap(A>0));
                    QSst(s,k)=std(Qmap(A>0));
                    
                    % n voxels surviving in QSM subcort areas:
                    nvoxROI(s,k)=sum(A(A>0));
                    k=k+1;
                end

                toc
                fprintf(num2str(s));
                flname=[outdir,'/MTEQSM_Harv_', num2str(cat{ct}),'_Tk', num2str(tk{tt}),'_ks', num2str(ks{kk}), '_', num2str(sus{ss}),'.mat']
                save(flname,'QS','QSst','QSMval','QSMgmst','nvoxROI','nvTis','nvox','tiv3noCB','tiv3', 'Qmap' )
                close all
            end
        end
    end
end

end

return


