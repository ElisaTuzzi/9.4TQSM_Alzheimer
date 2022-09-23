% returns the ROIs with a significant difference between AD and controls 
% plots the mean plaque load curves of AD and HC groups, as well as the curves of each single subject, from the cortical ROIs common to all subjects in the Acquisition Weighted space

clear all, clc
leg={'ADmean' 'HCmean' 'AD05' 'AD09' 'AD10' 'AD12' 'AD14' 'AD15' 'AD17' 'HC06' 'HC08' 'HC11' 'HC16' 'HC18' 'HC19' 'HC20'}
load('./LOAD/HCPL_AWsp.mat')  % .mat file with the plaque load from controls in AW space
load('./LOAD/ADPL_AWsp.mat')  % .mat file with the plaque load from AD patients in AW space
load ('./PL_load/AWIQSM_Harv.mat') % .mat file for ROIs definition. Contains the number of voxels with median QSM>thr of each ROI in AW space  

idx=find(all(nvoxROI~=0)); % finds the ROIs where the number of voxels with median QSM>thr differs from zero in all subjects 
AWIroivox=nvoxROI(:,idx); 
DK(idx);

%figure, imagesc(nvoxROI), colorbar, g=gca; set(g,'clim',[0 1e4]), ylabel('Subjs'), xlabel('ROIs'), title('#voxHarvROI all subjs') % image nvoxROIs all subjs
%AWIidx=find(all(nvoxROI~=0)) % find ROis in which all subj have voxels  

% cort:
cidx=idx(1:16);
cort=DK(cidx)
%%******

% values definition:
% Cort val:
val(1:2:96)=[1:48];
val(2:2:96)=[1:48]+100;
cortval=val(cidx); %

% Cortical ROIs defined in AW space
rname={'Frontal Pole' 'Inferior Temporal Gyrus temporo-occipital part' 'Inferior Frontal Gyrus pars opercularis' 'Precentral Gyrus' 'Postcentral Gyrus' 'Supramarginal Gurys anterior division' 'Supramarginal Gyrus posterior division' 'Angular Gyrus' 'Lateral Occipital Cortex superior division' 'Paracingulate Gyrus' 'Cingulate Gyrus anterior division' 'Cingulate Gyrus posterior division' 'Precuneous Cortex' 'Cuneal Cortex' 'Frontal Operculum Cortex' 'Central Opercular Cortex'};  
close all

% ttest:
admean=squeeze(mean(comm_fad_mtx_cort,1));
hcmean=squeeze(mean(comm_fhc_mtx_cort,1));
TTst=zeros(size(admean,1),2);
for m=1:size(c,2) %[1     3     5     7     9    10    11    13    15    17    19    21    23    25    27    28]
    figure(m), plot(thr.*1000,admean(m,:), 'r-','LineWidth', 3.5), hold on, plot(thr.*1000, hcmean(m,:),'b-','LineWidth', 3.5),title(num2str(rname{m})), 
    hold on, plot(thr.*1000, squeeze(comm_fad_mtx_cort(:,m,:)),  '-.','LineWidth', 1.5)
    hold on, plot(thr.*1000, squeeze(comm_fhc_mtx_cort(:,m,:)),  ':','LineWidth', 1),
    [h, p]=ttest2(admean(m,:), hcmean(m,:)); 
    text(30, 0.22, 'H='), text(32, 0.22, num2str(h)), text(36, 0.22, 'p='), text(38, 0.22, num2str(p)),
    %hold on, plot([thr(1).*1000 thr(end).*1000],[0.15 0.15],'k','LineWidth',2), 
    hold on, legend(leg), lgd = legend;
    lgd.NumColumns = 3;
    TT=[h, p];
    TTst(m,:)=TT;
end
%pidx=find(TTst(:,2)<0.0016)'; % after Bonferroni correction
pidx=find(TTst(:,2)<0.05)';
sigTTst=TTst(pidx,:);
h=sigTTst(:,1);
p=sigTTst(:,2);
Param=rname(pidx)';
a=p;
T = table(Param, a)

return  

  
