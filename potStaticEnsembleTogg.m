% =========================================================================
% this grogram summary data of toggle switch model with static ensemble
% data
% and plot the EWS
% =========================================================================

clear
clc
DataFolder = '/home/shan/Documents/MATLAB/criticlity/toggleSwitch/toggCLEEnsemble';

%get all the mat file in this folder
AllFile = dir(fullfile(DataFolder,filesep,'*.mat'));
% only keep the file name and ensure all vectors are columns
AllFile = {AllFile.name}.'; 
% get indices of files matching regular expression
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(AllFile,str,'once'));
pattern1 = '_se[0-9\._]+N20_2[68]';

%for the amplitude data
AllFile_amplitude = AllFile(FIND(pattern1));

%preallocate data strut to store final data
SummaryAmpli = struct('ampl',[],'a1',[],'var1',[],'var1_std',[],'cv1',[],'cv1_std',[],...
    'Fano1',[],'Fano1_std',[],'var2',[],'var2_std',[],'cv2',[],'cv2_std',[],'Fano2',[],'Fano2_std',[],...
    'corr',[],'corr_std',[]);

%sort the file into a order
SORT = @(str,cellArray)  regexp(cellArray,str,'match');
amplitude = nan(length(AllFile_amplitude),1);

s1 = '(?<= *se)[\d.]+';
[ap,amporder] = sort(amplitude);

for i0 = 1:length(AllFile_amplitude)
    amplitude(i0) = str2num(char(SORT(s1,char(AllFile_amplitude(i0)))));
end


%outliner index
outInx = [7,19]; % the index of outliners when Signame = 0.25

SummaryAmpli.ampl = sort(amplitude);
for i0 = 1:length(AllFile_amplitude)
    tempData = load(char(fullfile(DataFolder,filesep,AllFile_amplitude(amporder(i0)))));
    %get ride of outliners
    if(i0 == 3)
        tempData.meanVal(outInx(1),1,outInx(2)) = nan;
        tempData.meanVal(outInx(1),2,outInx(2)) = nan;
        tempData.variance(outInx(1),1,outInx(2)) = nan;
        tempData.variance(outInx(1),2,outInx(2)) = nan;
        tempData.corrCLE(outInx(1),1,outInx(2)) = nan;
        tempData.corrCLE(outInx(1),2,outInx(2)) = nan;
    end
%     SummaryAmpli.a1 = [SummaryAmpli.a1,tempData.a1list'];
%     SummaryAmpli.var1 = [SummaryAmpli.var1,nanmean(squeeze(tempData.variance(:,1,:)),2)];
%     SummaryAmpli.var1_std = [SummaryAmpli.var1_std,nanstd(squeeze(tempData.variance(:,1,:)),0,2)];
%     SummaryAmpli.var2 = [SummaryAmpli.var2,nanmean(squeeze(tempData.variance(:,2,:)),2)];
%     SummaryAmpli.var2_std = [SummaryAmpli.var2_std,nanstd(squeeze(tempData.variance(:,2,:)),0,2)];
%     
%     SummaryAmpli.cv1 = [SummaryAmpli.cv1,nanmean(sqrt(squeeze(tempData.variance(:,1,:)))./squeeze(tempData.meanVal(:,1,:)),2)];
%     SummaryAmpli.cv1_std = [SummaryAmpli.cv1_std,nanstd(sqrt(squeeze(tempData.variance(:,1,:)))./squeeze(tempData.meanVal(:,1,:)),0,2)];
%     SummaryAmpli.cv2 = [SummaryAmpli.cv2,nanmean(sqrt(squeeze(tempData.variance(:,2,:)))./squeeze(tempData.meanVal(:,2,:)),2)];
%     SummaryAmpli.cv2_std = [SummaryAmpli.cv2_std,nanstd(sqrt(squeeze(tempData.variance(:,2,:)))./squeeze(tempData.meanVal(:,2,:)),0,2)];
%     
%     SummaryAmpli.Fano1 = [SummaryAmpli.Fano1,nanmean(squeeze(tempData.variance(:,1,:))./squeeze(tempData.meanVal(:,1,:)),2)];
%     SummaryAmpli.Fano1_std = [SummaryAmpli.Fano1_std,nanstd(squeeze(tempData.variance(:,1,:))./squeeze(tempData.meanVal(:,1,:)),0,2)];
%     SummaryAmpli.Fano2 = [SummaryAmpli.Fano2,nanmean(squeeze(tempData.variance(:,2,:))./squeeze(tempData.meanVal(:,2,:)),2)];
%     SummaryAmpli.Fano2_std = [SummaryAmpli.Fano2_std,nanstd(squeeze(tempData.variance(:,2,:))./squeeze(tempData.meanVal(:,2,:)),0,2)];
    
    
%     SummaryAmpli.corr = [SummaryAmpli.corr,nanmean(tempData.corrCLE,2)];
%     SummaryAmpli.corr_std = [SummaryAmpli.corr_std,nanstd(tempData.corrCLE,0,2)];
    
    SummaryAmpli.a1{i0} = tempData.a1list';
    SummaryAmpli.var1{i0} = nanmean(squeeze(tempData.variance(:,1,:)),2);
    SummaryAmpli.var1_std{i0} = nanstd(squeeze(tempData.variance(:,1,:)),0,2);
    SummaryAmpli.var2{i0} = nanmean(squeeze(tempData.variance(:,2,:)),2);
    SummaryAmpli.var2_std{i0} = nanstd(squeeze(tempData.variance(:,2,:)),0,2);
    
    SummaryAmpli.cv1{i0} = nanmean(sqrt(squeeze(tempData.variance(:,1,:)))./squeeze(tempData.meanVal(:,1,:)),2);
    SummaryAmpli.cv1_std{i0} = nanstd(sqrt(squeeze(tempData.variance(:,1,:)))./squeeze(tempData.meanVal(:,1,:)),0,2);
    SummaryAmpli.cv2{i0} = nanmean(sqrt(squeeze(tempData.variance(:,2,:)))./squeeze(tempData.meanVal(:,2,:)),2);
    SummaryAmpli.cv2_std{i0} = nanstd(sqrt(squeeze(tempData.variance(:,2,:)))./squeeze(tempData.meanVal(:,2,:)),0,2);
    
    SummaryAmpli.Fano1{i0} = nanmean(squeeze(tempData.variance(:,1,:))./squeeze(tempData.meanVal(:,1,:)),2);
    SummaryAmpli.Fano1_std{i0} = nanstd(squeeze(tempData.variance(:,1,:))./squeeze(tempData.meanVal(:,1,:)),0,2);
    SummaryAmpli.Fano2{i0} = nanmean(squeeze(tempData.variance(:,2,:))./squeeze(tempData.meanVal(:,2,:)),2);
    SummaryAmpli.Fano2_std{i0} = nanstd(squeeze(tempData.variance(:,2,:))./squeeze(tempData.meanVal(:,2,:)),0,2);
    
    SummaryAmpli.corr{i0} = nanmean(tempData.corrCLE(:,:,1),2);
    SummaryAmpli.corr_std{i0} = nanstd(tempData.corrCLE(:,:,1),0,2);
end

%% plot the above data
%system size of this data set is 20
OMEGA = 20;


%load LNA data, make sure the loaded data of LNA is the corresponding one
Folder1 = '/home/shan/Documents/MATLAB/criticlity/toggleSwitch/figure and data';
filenameAmpLNA = [Folder1, filesep,'toggSwiLNATheoryAmpl_N20.mat'];
AmpLNA = load(filenameAmpLNA);

% index of selected different amplitude of EN
inx_se = [2,3,4];

%% set the font size, symbol size and line width
errorBarSize = 1.5;
errorMarkerSize  = 12;
LineWidth = 3;
labelFontSize = 36;
axisFontSize = 30;
ticketWidth = 1.5;

%size of figure
Xtick = 4:2:16;
figureSize = [0 0 8 6];
axis_var = [3,16,1e2,5e3];
axis_cv = [3,16,0.05,0.4];
axis_Fano = [3,16,1,25];
axis_corr = [3,16,-1,0.5];

colorSet = [3, 110, 184;224, 135, 51; 202, 39, 44; 0, 0, 0]/256;

%% plot the variance of B
figure(1)
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

hold on
for i0 = 1:3;
    errorbar(SummaryAmpli.a1{i0}, SummaryAmpli.var2{i0},SummaryAmpli.var2_std{i0},...
    'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'Color',colorSet(i0,:),'MarkerFaceColor','auto');
    plot(AmpLNA.a1,AmpLNA.varB(:,inx_se(i0)),'LineWidth',LineWidth,'Color',colorSet(i0,:))
end
plot(AmpLNA.a1,AmpLNA.var0B,'-','Color',colorSet(4,:),'LineWidth',LineWidth)
hold off
axis([3,16,100,1e4])
set(gca,'Xtick',Xtick,'YScale','log','LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
xlabel('a1','FontSize',labelFontSize,'FontName','Helvetica')
ylabel('variacne of B','FontSize',labelFontSize,'FontName','Helvetica')
fileNamePdf = [DataFolder,filesep,'GlpVar.pdf'];
fileNameFig = [DataFolder,filesep,'GlpVar.fig'];
print('-dpdf',fileNamePdf)
saveas(gcf,fileNameFig)
%% plot cv
figure(2)
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
cvB = nan(size(AmpLNA.varB,1),length(inx_se));
for i0 = 1:length(inx_se)
    cvB(:,i0) = sqrt(AmpLNA.varB(:,inx_se(i0)))./AmpLNA.meanB/OMEGA;
    plot(AmpLNA.a1,cvB(:,i0),'LineWidth',LineWidth,'Color',colorSet(i0,:))
    errorbar(SummaryAmpli.a1{i0}, SummaryAmpli.cv2{i0},SummaryAmpli.cv2_std{i0},...
    'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'MarkerFaceColor',...
    'auto','Color',colorSet(i0,:));
end
cv0B = sqrt(AmpLNA.var0B)./AmpLNA.meanB/OMEGA;
plot(AmpLNA.a1,cv0B,'-','LineWidth',LineWidth,'Color',colorSet(4,:))
hold off
axis(axis_cv)
set(gca,'Xtick',Xtick,'Ytick',0:0.1:0.4,'LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
xlabel('a1','FontSize',labelFontSize,'FontName','Helvetica')
ylabel('coefficient of variation','FontSize',labelFontSize,'FontName','Helvetica')
fileNameCvPdf = [DataFolder,filesep,'Glpcv.pdf'];
fileNameFanoFig = [DataFolder,filesep,'Glpcv.fig'];
print('-dpdf',fileNameCvPdf)
saveas(gcf,fileNameFanoFig)


%% plot Fano factor
figure(3)
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
FanoB = nan(size(AmpLNA.varB,1),length(inx_se));
for i0 = 1:length(inx_se)
    FanoB(:,i0) = AmpLNA.varB(:,inx_se(i0))./AmpLNA.meanB/OMEGA;
    plot(AmpLNA.a1,FanoB(:,i0),'LineWidth',LineWidth,'Color',colorSet(i0,:))
    errorbar(SummaryAmpli.a1{i0}, SummaryAmpli.Fano2{i0},SummaryAmpli.Fano2_std{i0},...
        'o','MarkerSize',errorMarkerSize,'Color',colorSet(i0,:),'LineWidth',errorBarSize,...
    'MarkerFaceColor','auto');
end
Fano0B = AmpLNA.var0B./AmpLNA.meanB/OMEGA;
plot(AmpLNA.a1,Fano0B,'k-','LineWidth',LineWidth,'Color',colorSet(4,:))
hold off
axis([3,16,0.5,25])
set(gca,'Xtick',Xtick,'YScale','log','LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
xlabel('a1','FontSize',labelFontSize,'FontName','Helvetica')
ylabel('Fano factor','FontSize',labelFontSize,'FontName','Helvetica')
fileNameFanoPdf = [DataFolder,filesep,'GlpFano.pdf'];
fileNameFanoFig = [DataFolder,filesep,'GlpFano.fig'];
print('-dpdf',fileNameFanoPdf)
saveas(gcf,fileNameFanoFig)

%% plot correlation coefficient
figure(4)
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
for i0 = 1:length(inx_se);
    errorbar(SummaryAmpli.a1{i0}, SummaryAmpli.corr{i0},SummaryAmpli.corr_std{i0},...
    'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'Color',colorSet(i0,:),...
    'MarkerFaceColor','auto');
    plot(AmpLNA.a1,AmpLNA.corr(:,inx_se(i0)),'LineWidth',LineWidth,'Color',colorSet(i0,:))
end
plot(AmpLNA.a1,AmpLNA.corr0,'k-','LineWidth',LineWidth,'Color',colorSet(4,:))
hold off
axis(axis_corr)
set(gca,'Xtick',Xtick,'Ytick',-1:0.2:0.4,'LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
xlabel('a1','FontSize',labelFontSize,'FontName','Helvetica')
ylabel('correlation coeff','FontSize',labelFontSize,'FontName','Helvetica')
fileNameCurrPdf = [DataFolder,filesep,'Glpcorr.pdf'];
fileNameCurrFig = [DataFolder,filesep,'Glpcorr.fig'];
print('-dpdf',fileNameCurrPdf)
saveas(gcf,fileNameCurrFig)


%% also plot the figure for the other variable, here is gene A
%% plot the variance of A
figure(5)
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

hold on
for i0 = 1:3;
    errorbar(SummaryAmpli.a1{i0}, SummaryAmpli.var1{i0},SummaryAmpli.var1_std{i0},...
    'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'Color',colorSet(i0,:),'MarkerFaceColor','auto');
    plot(AmpLNA.a1,AmpLNA.varA(:,inx_se(i0)),'LineWidth',LineWidth,'Color',colorSet(i0,:))
end
plot(AmpLNA.a1,AmpLNA.var0A,'-','Color',colorSet(4,:),'LineWidth',LineWidth)
hold off
axis([3,16,100,1e4])
set(gca,'Xtick',Xtick,'YScale','log','LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
xlabel('a1','FontSize',labelFontSize,'FontName','Helvetica')
ylabel('variacne of B','FontSize',labelFontSize,'FontName','Helvetica')
fileNamePdf = [DataFolder,filesep,'GlpVar_A.pdf'];
fileNameFig = [DataFolder,filesep,'GlpVar_A.fig'];
print('-dpdf',fileNamePdf)
saveas(gcf,fileNameFig)



