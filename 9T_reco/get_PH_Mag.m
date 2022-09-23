% Reconstracts mag and phase images from high resolution acquisition
% weighted images acquired at 9.4T
% Uses .rms files obtained with SiTools.
% Please change the path of your data in the elisa.bat file to use SiTools

%% Rolf Pohmann, Max Planck Institute for Biological Cybernetics %%
 

function [mag,ph,mtx]=get_PH_mag(filename)  % input file has to be .rms

ima=ReadSiToolData(filename);%data from all the 30 coils;  
[nx ny nz nc] = size(squeeze(ima));
d = shiftdim(single(squeeze(ima)),3); %single instead of double for faster reco
%d = single(squeeze(ima));
d = d(:,:,:,4:end-2); %crop in a sensibile way - outer slices for 3D acq (nz should be even), other dimensions to avoid overflooding the memory
mtx = size(d); mtx = mtx(2:end); %adapt the matrix size
%   Adaptive recon based on Walsh et al.
%   Walsh DO, Gmitro AF, Marcellin MW.
%   Adaptive reconstruction of phased array MR imagery. 
%   Magn Reson Med. 2000 May;43(5):682-90.
%
%    and
%
%   Mark Griswold, David Walsh, Robin Heidemann, Axel Haase, Peter Jakob. 
%   The Use of an Adaptive Reconstruction for Array Coil Sensitivity Mapping and Intensity Normalization, 
%   Proceedings of the Tenth  Scientific Meeting of the International Society for Magnetic Resonance in Medicine pg 2410 (2002)
%
%   better iterative phase correction based on 
%   Inati et al., ISMRM 2014, #4407

recon = adaptiveCombine(d,[6 6 3]); %Inati 2014; Matlab/mex-function for phase sensitive combination of coil images

%recon = adaptiveCombine(d);
Mag=abs(recon);
Ph=angle(recon);

end
