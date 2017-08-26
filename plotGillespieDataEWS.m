% =========================================================================
% this grogram summary data of toggle switch modeled by Gillespie alogrithm
% and plot the EWS
% =========================================================================

clear
clc
DataFolder = '/home/shan/Documents/MATLAB/criticlity/toggleSwitch/toggGlpDataN20';

%get all the mat file in this folder
AllFile = dir(fullfile(DataFolder,filesep,'*.mat'));
% only keep the file name and ensure all vectors are columns
AllFile = {AllFile.name}.'; 
% get indices of files matching regular expression
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(AllFile,str,'once'));
pattern1 = '_tau-100_se[0-9\._]+N20';

%for the amplitude data
AllFile_amplitude = AllFile(FIND(pattern1));

%preallocate data strut to store final data
SummaryAmpli = struct('ampl',[],'a1',[],'var1',[],'var1_std',[],'cv1',[],'cv1_std',[],...
    'Fano1',[],'Fano1_std',[],'lagAuto1',[],'lagAuto1_std',[],...
    'var2',[],'var2_std',[],'cv2',[],'cv2_std',[],'Fano2',[],'Fano2_std',[],...
    'lagAuto2',[],'lagAuto2_std',[],'corr',[],'corr_std',[]);

%sort the file into a order
SORT = @(str,cellArray)  regexp(cellArray,str,'match');
amplitude = nan(length(AllFile_amplitude),1);

s1 = '(?<= *se)[\d.]+';
[ap,amporder] = sort(amplitude);

for i0 = 1:length(AllFile_amplitude)
    amplitude(i0) = str2num(char(SORT(s1,char(AllFile_amplitude(i0)))));
end


SummaryAmpli.ampl = ap;
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

figureSize = [0 0 8 6];

colorSet = [3, 110, 184;224, 135, 51; 202, 39, 44; 0, 0, 0]/256;

%% plot the variance of B
figure(1)
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

hold on
for i0 = 1:3
    errorbar(SummaryAmpli.a1(:,i0), SummaryAmpli.var2(:,i0),SummaryAmpli.var2_std(:,i0),...
    'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'Color',colorSet(i0,:),'MarkerFaceColor','auto');
    plot(AmpLNA.a1,AmpLNA.varB(:,inx_se(i0)),'LineWidth',LineWidth,'Color',colorSet(i0,:))
end
plot(AmpLNA.a1,AmpLNA.var0B,'-','Color',colorSet(4,:),'LineWidth',LineWidth)
hold off
axis([3,16,100,8000])
set(gca,'Xtick',3:2:16,'Ytick',[100,1000,5000],'YScale','log','LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
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
    errorbar(SummaryAmpli.a1(:,i0), SummaryAmpli.cv2(:,i0),SummaryAmpli.cv2_std(:,i0),...
    'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'MarkerFaceColor',...
    'auto','Color',colorSet(i0,:));
end
cv0B = sqrt(AmpLNA.var0B)./AmpLNA.meanB/OMEGA;
% cvB = sqrt(AmpLNA.varB(:,inx_se))./AmpLNA.meanB(:,inx_se)
% errorbar(SummaryAmpli.a1, SummaryAmpli.cv2,SummaryAmpli.cv2_std,...
%     'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'MarkerFaceColor',...
%     'auto','Color',colorSet(i0));
% plot(AmpLNA.a1,cvB,'LineWidth',LineWidth)
plot(AmpLNA.a1,cv0B,'-','LineWidth',LineWidth,'Color',colorSet(4,:))
hold off
axis([3,16,0,0.4])
set(gca,'Xtick',3:2:16,'Ytick',0:0.1:0.4,'LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
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
    errorbar(SummaryAmpli.a1(:,i0), SummaryAmpli.Fano2(:,i0),SummaryAmpli.Fano2_std(:,i0),...
        'o','MarkerSize',errorMarkerSize,'Color',colorSet(i0,:),'LineWidth',errorBarSize,...
    'MarkerFaceColor','auto');
end
Fano0B = AmpLNA.var0B./AmpLNA.meanB/OMEGA;
% cvB = sqrt(AmpLNA.varB(:,inx_se))./AmpLNA.meanB(:,inx_se)
% errorbar(SummaryAmpli.a1(:,inx_se), SummaryAmpli.Fano2(:,inx_se),SummaryAmpli.Fano2_std(:,inx_se),...
%     'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'MarkerFaceColor','auto');
% plot(AmpLNA.a1,FanoB,'LineWidth',LineWidth)
plot(AmpLNA.a1,Fano0B,'k-','LineWidth',LineWidth,'Color',colorSet(4,:))
hold off
axis([3,16,1,32])
set(gca,'Xtick',3:2:16,'Ytick',[1,10,20],'YScale','log','LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
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
for i0 = 1:length(inx_se)
    errorbar(SummaryAmpli.a1(:,i0), SummaryAmpli.corr(:,i0),SummaryAmpli.corr_std(:,i0),...
    'o','MarkerSize',errorMarkerSize,'LineWidth',errorBarSize,'Color',colorSet(i0,:),...
    'MarkerFaceColor','auto');
    plot(AmpLNA.a1,AmpLNA.corr(:,inx_se(i0)),'LineWidth',LineWidth,'Color',colorSet(i0,:))
end
plot(AmpLNA.a1,AmpLNA.corr0,'k-','LineWidth',LineWidth,'Color',colorSet(4,:))
hold off
axis([3,16,-1,0.4])
set(gca,'Xtick',3:2:16,'Ytick',-1:0.25:0.25,'LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
xlabel('a1','FontSize',labelFontSize,'FontName','Helvetica')
ylabel('correlation coeff','FontSize',labelFontSize,'FontName','Helvetica')
fileNameCurrPdf = [DataFolder,filesep,'Glpcorr.pdf'];
fileNameCurrFig = [DataFolder,filesep,'Glpcorr.fig'];
print('-dpdf',fileNameCurrPdf)
saveas(gcf,fileNameCurrFig)


%% plot auto correlation coefficient
figure(5)
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
for i0 = 1:length(inx_se)
    errorbar(SummaryAmpli.a1(:,i0), SummaryAmpli.lagAuto2(:,i0),SummaryAmpli.lagAuto2_std(:,i0),...
        'o','MarkerSize',errorMarkerSize,'Color',colorSet(i0,:),'LineWidth',errorBarSize);
    plot(AmpLNA.a1,AmpLNA.lagB(:,inx_se(i0)),'LineWidth',LineWidth,'Color',colorSet(i0,:))
end
plot(AmpLNA.a1,AmpLNA.lag0B,'k-','LineWidth',LineWidth,'Color',colorSet(4,:))
hold off
axis([3,16,0,1])
set(gca,'Xtick',3:2:16,'YLim',[0.3,1],'LineWidth',ticketWidth,'FontName','Helvetica','FontSize',axisFontSize)
xlabel('a1','FontSize',labelFontSize,'FontName','Helvetica')
ylabel('autocorrelation coeff','FontSize',labelFontSize,'FontName','Helvetica')
fileNameLagPdf = [DataFolder,filesep,'Glplag.pdf'];
fileNameLagFig = [DataFolder,filesep,'Glplag.fig'];
print('-dpdf',fileNameLagPdf)
saveas(gcf,fileNameLagFig)