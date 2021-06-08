function stereoAcuity(file)
% This function analyze the results of eRDS7 psi-marg-grid algorithm to visualize the
% data and the threshold
%
% data structure is (in psi.history)
%   1       current psi trial
%   2       current disparity shown in "
%   3       psi.correct or not
%   4-6     current estimates for thres / slope / neg_slope using sum method
%   7      final 75% threshold estimate using sum parameters
%   8       trial # (different from psi.trial)
% 9-10      grid range on threshold

close all;
[eRDSpath,~]=fileparts(fileparts(mfilename('fullpath'))); %path to erds folder
[~,filename,ext] = fileparts(file); 
if isempty(ext); ext='.mat'; end
load(fullfile(eRDSpath,'dataFiles',[filename,ext]),'psi1','psi2','expe');
expe.filename = filename;
addpath(fullfile(eRDSpath,'eRDS_functions'));
expe.eRDSpath = eRDSpath; dispi('Data file: ',expe.filename);
dispi('Duration: ',round(expe.duration,1),' min');
disp('-------------------------------------------------------');
dispi('    composite threshold');
disp('-------------------------------------------------------');
psi=psi1;        psi.history = [psi1.history; psi2.history]; psi.history = sortrows(psi.history,8);
figure('Color', 'w','Units','normalized','Position',[0 0 0.9 0.9]);
psi = recomputeFromData4(psi, 1);
psi.final_threshold=round(min(psi.maxAllowerThreshold,psi.threshold),1);
dispi('Final threshold: ',psi.final_threshold,' arcsec');
saveas(gcf,fullfile(expe.eRDSpath,'figures', [expe.filename,'.fig']));
saveas(gcf,fullfile(expe.eRDSpath,'figures', [expe.filename,'.png']));
end















