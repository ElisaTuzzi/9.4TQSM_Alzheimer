% returns the ROIs with a significant difference between AD and controls for each method
% plots the mean plaque load curves of the AD and HC groups, as well as the curves of each single subject in the 16 cortical common ROIs

close all, clear all, clc
leg={'ADmean' 'HCmean' 'AD05' 'AD09' 'AD10' 'AD12' 'AD14' 'AD15' 'AD17' 'HC06' 'HC08' 'HC11' 'HC16' 'HC18' 'HC19' 'HC20'}

%**********************
%MTE NAV cor Tk-3 ks1.2:

load(['./LOAD/ADPL_MTEsp_Tk3_CoR_ks12_PARAM.mat'])  % .mat file with the plaque load from the AD patients in MTE space
load(['./LOAD/HCPL_MTEsp_Tk3_CoR_ks12_PARAM.mat'])  % .mat file with the plaque load from controls in MTE space
load('./LOAD/MTEQSM_Harv.mat') % .mat file for ROIs definition. Contains the number of voxels with median QSM>thr of each ROI in MTE space  
load('./Pl_load/HarvROis.mat')   % Harvard ROIs 

MTEidx=find(all(nvoxROI~=0));  % finds the ROIs where the number of voxels with median QSM>thr differs from zero in all subjects 
MTEroivox=nvoxROI(:,MTEidx) ;% 
DK(MTEidx);

% cort:
% index definition:
cidx=MTEidx(1:78); % cortical indexes of Harvard ROIs identified
cort=DK(cidx); 
cmnidx=[1,5:7,13,15:18,24:28,32,33]; % indexes of ROIs common to MTE and AW space 
DK(cmnidx);


% indexes of bilateral cortical ROIs:
val(1:2:96)=[1:48]; 
val(2:2:96)=[1:48]+100;
cortval=val(cidx); 
pp=find(cortval<100); % 1x39
unilatCort=cort(pp);
cmnROIs=unilatCort(cmnidx);

%**********************
rname={'Front Pole', 'Insular Cort', 'Sup Front Gyr', 'Mid Front Gyr' ,'Inf Temp Gyr temp-occ part', 'Inf Front Gyr Pars Oper', 'Precentral Gyr' ,'Temp Fus Cort, post div' ,'Sup Temp Gyr, ant div', 'Sup Temp Gyr, post div' ,'Mid Temp Gyr, post div' ,'Mid Temp Gyr, temp-occ div' ,'Postcentr Gyr', 'Sup Par lob' ,'Supram Gyr ant div' ,'Supram Gyr post div', 'Angular Gyrus' ,'Later Occ Cort sup div' ,'Later Occ Cort, inf div','Calc', 'Front Medial Cort' ,'SMAf' ,'Subcall Cort', 'Paracing Gyr', 'Cing Gyr ant div' ,'Cing Gyr post div', 'Precun Cort' ,'Cuneal Cort', 'Front Orb Cort' ,'Para Hipp Gyr, post div' ,'Lingual Gyrus', 'Front Operc Cort' ,'Central Operc Cort', 'Pariet Operc Cort' ,'Planum Polare', 'Heschls Temp' ,'Planum Temporale' ,'SupraCalc Cort' ,'Occ Pole'} ;

% ttest:
admean=squeeze(mean(comm_fad_mtx_cort(:,1:2:end,:),1));
hcmean=squeeze(mean(comm_fhc_mtx_cort(:,1:2:end,:),1));
comm_fad_mtx_cort_new=comm_fad_mtx_cort(:,1:2:end,:);
comm_fhc_mtx_cort_new=comm_fhc_mtx_cort(:,1:2:end,:);

TTst=zeros(16,2); %    plots plaque load in the 16 ROIs common to AW and MTE space
for m=1:size(cmnidx,2)
    figure(m), plot(thr.*1000,admean(cmnidx(m),:), 'r-','LineWidth', 3.5), hold on, plot(thr.*1000, hcmean(cmnidx(m),:),'b-','LineWidth', 3.5),title(num2str(cmnROIs{m})), 
    hold on, plot(thr.*1000, squeeze(comm_fad_mtx_cort_new(:,m,:)),  '-.','LineWidth', 1.5)
    hold on, plot(thr.*1000, squeeze(comm_fhc_mtx_cort_new(:,m,:)),  ':','LineWidth', 1),
    [h, p]=ttest2(admean(cmnidx(m),:), hcmean(cmnidx(m),:)); 
    %text(30, 0.22, 'H='), text(32, 0.22, num2str(h(m),:)), text(36, 0.22, 'p='), text(38, 0.22, num2str(p(m),:))
    %hold on, plot([thr(1).*1000 thr(end).*1000],[0.15 0.15],'k','LineWidth',2), hold on, 
    legend(leg), lgd = legend;
    lgd.NumColumns = 3;
    TT=[h, p];
    TTst(m,:)=TT;
end
pidx=find(TTst(:,2)<0.05)';
sigTTst=TTst(pidx,:);
h=sigTTst(:,1);
p=sigTTst(:,2);
Param_Cor_Tk3_ks12_p05=rname(pidx)';
a=p;
T = table(Param_Cor_Tk3_ks12_p05, a)
%%


%************************
%MTE PARAM Cor Tk12 ks16:

close all, clear all,
load(['./LOAD/ADPL_MTEsp_Tk12_CoR_ks16_PARAM.mat'])  % plaque load of controls in MTE space
load(['./LOAD/HCPL_MTEsp_Tk12_CoR_ks16_PARAM.mat'])  % plaque load of AD patients in MTE space
load('./LOAD/MTEQSM_Harv.mat') % n voxels of each ROI in MTE space  
load('./PL_load/HarvROis.mat')   % Harvard ROIs 

MTEidx=find(all(nvoxROI~=0));  % finds the ROIs where the number of voxels with median QSM>thr differs from zero in all subjects 
MTEroivox=nvoxROI(:,MTEidx) ;% 
DK(MTEidx);

% cort:
% index definition:
cidx=MTEidx(1:78); % cortical indexes of Harvard ROIs identified
cort=DK(cidx); 
cmnidx=[1,5:7,13,15:18,24:28,32,33]; % indexes of ROIs common to MTE and AW space 
DK(cmnidx);


% indexes of bilateral cortical ROIs:
val(1:2:96)=[1:48]; 
val(2:2:96)=[1:48]+100;
cortval=val(cidx); 
pp=find(cortval<100); % 1x39
unilatCort=cort(pp);
cmnROIs=unilatCort(cmnidx);

%**********************
rname={'Front Pole', 'Insular Cort', 'Sup Front Gyr', 'Mid Front Gyr' ,'Inf Temp Gyr temp-occ part', 'Inf Front Gyr Pars Oper', 'Precentral Gyr' ,'Temp Fus Cort, post div' ,'Sup Temp Gyr, ant div', 'Sup Temp Gyr, post div' ,'Mid Temp Gyr, post div' ,'Mid Temp Gyr, temp-occ div' ,'Postcentr Gyr', 'Sup Par lob' ,'Supram Gyr ant div' ,'Supram Gyr post div', 'Angular Gyrus' ,'Later Occ Cort sup div' ,'Later Occ Cort, inf div','Calc', 'Front Medial Cort' ,'SMAf' ,'Subcall Cort', 'Paracing Gyr', 'Cing Gyr ant div' ,'Cing Gyr post div', 'Precun Cort' ,'Cuneal Cort', 'Front Orb Cort' ,'Para Hipp Gyr, post div' ,'Lingual Gyrus', 'Front Operc Cort' ,'Central Operc Cort', 'Pariet Operc Cort' ,'Planum Polare', 'Heschls Temp' ,'Planum Temporale' ,'SupraCalc Cort' ,'Occ Pole'} ;

% ttest:
admean=squeeze(mean(comm_fad_mtx_cort(:,1:2:end,:),1));
hcmean=squeeze(mean(comm_fhc_mtx_cort(:,1:2:end,:),1));
comm_fad_mtx_cort_new=comm_fad_mtx_cort(:,1:2:end,:);
comm_fhc_mtx_cort_new=comm_fhc_mtx_cort(:,1:2:end,:);

%ttest
TTst=zeros(16,2); %    plots plaque load in the 16 ROIs common to AW and MTE space
for m=1:size(cmnidx,2)
    figure(m), plot(thr.*1000,admean(cmnidx(m),:), 'r-','LineWidth', 3.5), hold on, plot(thr.*1000, hcmean(cmnidx(m),:),'b-','LineWidth', 3.5),title(num2str(cmnROIs{m})), 
    hold on, plot(thr.*1000, squeeze(comm_fad_mtx_cort_new(:,m,:)),  '-.','LineWidth', 1.5)
    hold on, plot(thr.*1000, squeeze(comm_fhc_mtx_cort_new(:,m,:)),  ':','LineWidth', 1),
    [h, p]=ttest2(admean(cmnidx(m),:), hcmean(cmnidx(m),:)); 
    %text(30, 0.22, 'H='), text(32, 0.22, num2str(h(m),:)), text(36, 0.22, 'p='), text(38, 0.22, num2str(p(m),:))
    %hold on, plot([thr(1).*1000 thr(end).*1000],[0.15 0.15],'k','LineWidth',2), hold on, 
    legend(leg), lgd = legend;
    lgd.NumColumns = 3;
    TT=[h, p];
    TTst(m,:)=TT;
end
pidx=find(TTst(:,2)<0.05)';
sigTTst=TTst(pidx,:);
h=sigTTst(:,1);
p=sigTTst(:,2);
Param_Cor_Tk3_ks12_p05=rname(pidx)';
a=p;
T = table(Param_Cor_Tk12_ks16_p05, a)

return


