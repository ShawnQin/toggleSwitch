% ================================================================
% this porgram summary the simulated results and export them into 
% a file that R can laod and redraw figures
% last revised on 11/10/2016
% ================================================================

clear
clc
DataFolder = '/home/shan/Documents/MATLAB/criticlity/toggleSwitch';

%get all the mat file in this folder
AllFile = dir(fullfile(DataFolder,filesep,'*.mat'));
% only keep the file name and ensure all vectors are columns
AllFile = {AllFile.name}.'; 
% get indices of files matching regular expression
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(AllFile,str,'once'));
pattern1 = '_tau-100_se[0-9\._]+N20';

%for the amplitude data
AllFile_amplitude = AllFile(FIND(pattern1));

% all the correlation time scale
pattern2 = 'se0.1_N20';
AllFile_timescale = AllFile(FIND(pattern2));

% summary of extrinsic noise with different amplitude of noise, stored in a
% data structure
SummaryAmpli = struct('ampl',[],'a1',[],'var1',[],'var1_std',[],'cv1',[],'cv1_std',[],...
    'Fano1',[],'Fano1_std',[],'lagAuto1',[],'lagAuto1_std',[],...
    'var2',[],'var2_std',[],'cv2',[],'cv2_std',[],'Fano2',[],'Fano2_std',[],...
    'lagAuto2',[],'lagAuto2_std',[],'corr',[],'corr_std',[]);

SummaryTimeScale = struct('timeScale',[],'a1',[],'var1',[],'var1_std',[],...
    'cv1',[],'cv1_std',[],'lagAuto1',[],'lagAuto1_std',[],'var2',[],'var2_std',[],...
    'cv2',[],'cv2_std',[],'lagAuto2',[],'lagAuto2_std',[],'corr',[],'corr_std',[],...
    'Fano1',[],'Fano1_std',[],'Fano2',[],'Fano2_std',[]);

%sort the file into a order
SORT = @(str,cellArray)  regexp(cellArray,str,'match');
amplitude = nan(length(AllFile_amplitude),1);
timescale = nan(length(AllFile_timescale),1);
s1 = '(?<= *se)[\d.]+';
% s2 = '[\d.]+(?= _se*)';
s2 = '(?<= *tau-)[\d.]+';
for i0 = 1:length(AllFile_amplitude)
    amplitude(i0) = str2num(char(SORT(s1,char(AllFile_amplitude(i0)))));
end
[ap,amporder] = sort(amplitude);

for i0 = 1:length(AllFile_timescale)
    timescale(i0) = str2num(char(SORT(s2,char(AllFile_timescale(i0)))));
end
[ts,timeorder] = sort(timescale);

SummaryAmpli.ampl = ap;
SummaryTimeScale.timeScale = ts;
for i0 = 1:length(AllFile_amplitude)
    tempData = load(char(fullfile(DataFolder,filesep,AllFile_amplitude(amporder(i0)))));
    SummaryAmpli.a1 = [SummaryAmpli.a1,tempData.a1list'];
    SummaryAmpli.var1 = [SummaryAmpli.var1,nanmean(tempData.variance1,2)];
    SummaryAmpli.var1_std = [SummaryAmpli.var1_std,nanstd(tempData.variance1,0,2)];
    SummaryAmpli.var2 = [SummaryAmpli.var2,nanmean(tempData.variance2,2)];
    SummaryAmpli.var2_std = [SummaryAmpli.var2_std,nanstd(tempData.variance2,0,2)];
    
    SummaryAmpli.cv1 = [SummaryAmpli.cv1,nanmean(sqrt(tempData.variance1)./tempData.meanVal1,2)];
    SummaryAmpli.cv1_std = [SummaryAmpli.cv1_std,nanstd(sqrt(tempData.variance1)./tempData.meanVal1,0,2)];
    SummaryAmpli.cv2 = [SummaryAmpli.cv2,nanmean(sqrt(tempData.variance2)./tempData.meanVal2,2)];
    SummaryAmpli.cv2_std = [SummaryAmpli.cv2_std,nanstd(sqrt(tempData.variance2)./tempData.meanVal2,0,2)];
    
    SummaryAmpli.Fano1 = [SummaryAmpli.Fano1,nanmean(tempData.variance1./tempData.meanVal1,2)];
    SummaryAmpli.Fano1_std = [SummaryAmpli.Fano1_std,nanstd(tempData.variance1./tempData.meanVal1,0,2)];
    SummaryAmpli.Fano2 = [SummaryAmpli.Fano2,nanmean(tempData.variance2./tempData.meanVal2,2)];
    SummaryAmpli.Fano2_std = [SummaryAmpli.Fano2_std,nanstd(tempData.variance2./tempData.meanVal2,0,2)];
    
    SummaryAmpli.lagAuto1 = [SummaryAmpli.lagAuto1,nanmean(tempData.lagAuto1,2)];
    SummaryAmpli.lagAuto1_std = [SummaryAmpli.lagAuto1_std,nanstd(tempData.lagAuto1,0,2)];
    SummaryAmpli.lagAuto2 = [SummaryAmpli.lagAuto2,nanmean(tempData.lagAuto2,2)];
    SummaryAmpli.lagAuto2_std = [SummaryAmpli.lagAuto2_std,nanstd(tempData.lagAuto2,0,2)];
    
    SummaryAmpli.corr = [SummaryAmpli.corr,nanmean(tempData.corrCLE,2)];
    SummaryAmpli.corr_std = [SummaryAmpli.corr_std,nanstd(tempData.corrCLE,0,2)];
    
end

for j0 = 1:length(AllFile_timescale)
    tempData = load(char(fullfile(DataFolder,filesep,AllFile_timescale(timeorder(j0)))));
    SummaryTimeScale.a1 = [SummaryTimeScale.a1,tempData.a1list'];
    SummaryTimeScale.var1 = [SummaryTimeScale.var1,nanmean(tempData.variance1,2)];
    SummaryTimeScale.var1_std = [SummaryTimeScale.var1_std,nanstd(tempData.variance1,0,2)];
    SummaryTimeScale.var2 = [SummaryTimeScale.var2,nanmean(tempData.variance2,2)];
    SummaryTimeScale.var2_std = [SummaryTimeScale.var2_std,nanstd(tempData.variance2,0,2)];
    
    SummaryTimeScale.cv1 = [SummaryTimeScale.cv1,nanmean(sqrt(tempData.variance1)./tempData.meanVal1,2)];
    SummaryTimeScale.cv1_std = [SummaryTimeScale.cv1_std,nanstd(sqrt(tempData.variance1)./tempData.meanVal1,0,2)];
    SummaryTimeScale.cv2 = [SummaryTimeScale.cv2,nanmean(sqrt(tempData.variance2)./tempData.meanVal2,2)];
    SummaryTimeScale.cv2_std = [SummaryTimeScale.cv2_std,nanstd(sqrt(tempData.variance2)./tempData.meanVal2,0,2)];
    
    SummaryTimeScale.Fano1 = [SummaryTimeScale.Fano1,nanmean(tempData.variance1./tempData.meanVal1,2)];
    SummaryTimeScale.Fano1_std = [SummaryTimeScale.Fano1_std,nanstd(tempData.variance1./tempData.meanVal1,0,2)];
    SummaryTimeScale.Fano2 = [SummaryTimeScale.Fano2,nanmean(tempData.variance2./tempData.meanVal2,2)];
    SummaryTimeScale.Fano2_std = [SummaryTimeScale.Fano2_std,nanstd(tempData.variance2./tempData.meanVal2,0,2)];
    
    SummaryTimeScale.lagAuto1 = [SummaryTimeScale.lagAuto1,nanmean(tempData.lagAuto1,2)];
    SummaryTimeScale.lagAuto1_std = [SummaryTimeScale.lagAuto1_std,nanstd(tempData.lagAuto1,0,2)];
    SummaryTimeScale.lagAuto2 = [SummaryTimeScale.lagAuto2,nanmean(tempData.lagAuto2,2)];
    SummaryTimeScale.lagAuto2_std = [SummaryTimeScale.lagAuto2_std,nanstd(tempData.lagAuto2,0,2)];
    
    
    SummaryTimeScale.corr = [SummaryTimeScale.corr,nanmean(tempData.corrCLE,2)];
    SummaryTimeScale.corr_std = [SummaryTimeScale.corr_std,nanstd(tempData.corrCLE,0,2)];
end

%save the summary data
saveFile1 = 'ExtriFig1SimuAmp.mat';
saveFile2 = 'ExtriFig1SimuTime.mat';
save(saveFile1,'-struct','SummaryAmpli')
save(saveFile2,'-struct','SummaryTimeScale')