% this program replot all the figures in manuscripts
% revised on 11/13/2016

clear
clc


%% load data for figure 1
Folder1 = '/home/shan/Documents/MATLAB/criticlity/toggleSwitch/figure and data';
filenameAmpSimu =[Folder1, filesep,'ExtriFig1SimuAmp.mat'];
filenameTimeSimu = [Folder1, filesep,'ExtriFig1SimuTime.mat'];
filenameAmpLNA = [Folder1, filesep,'toggSwiLNATheoryAmpl_N20.mat'];
filenameTimeLNA = [Folder1, filesep,'toggSwiLNATheoryTimeScale_N20_se0.1.mat'];

AmpSim = load(filenameAmpSimu);
TimeSimu  = load(filenameTimeSimu);
AmpLNA = load(filenameAmpLNA);
TimeLNA = load(filenameTimeLNA);

%system size of this data set is 20
OMEGA = 20;
se = AmpSim.ampl; % different strength of extrinsic noise
inx_se = 1:3;
time_scale = TimeSimu.timeScale;
inx_ts = 1:3;

%% plot the variance of B
figure(1)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');

hold on
h1 = errorbar(AmpSim.a1(:,inx_se), AmpSim.var2(:,inx_se),AmpSim.var2_std(:,inx_se),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(AmpLNA.a1,AmpLNA.varB(:,inx_se+1),'LineWidth',2)
plot(AmpLNA.a1,AmpLNA.var0B,'k-','LineWidth',2)
hold off
axis([3,16,100,8000])
set(gca,'Xtick',3:2:16,'Ytick',[100,1000,5000],'YScale','log','LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('variacne of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNamePdf = [Folder1,filesep,'fig1Var.pdf'];
fileNameFig = [Folder1,filesep,'fig1Var.fig'];
print('-dpdf',fileNamePdf)
saveas(gcf,fileNameFig)

%% plot cv
figure(2)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
cvB = nan(size(AmpLNA.varB,1),length(inx_se));
for i0 = 1:length(inx_se)
    cvB(:,i0) = sqrt(AmpLNA.varB(:,i0+1))./AmpLNA.meanB/OMEGA;
end
cv0B = sqrt(AmpLNA.var0B)./AmpLNA.meanB/OMEGA;
% cvB = sqrt(AmpLNA.varB(:,inx_se))./AmpLNA.meanB(:,inx_se)
h1 = errorbar(AmpSim.a1(:,inx_se), AmpSim.cv2(:,inx_se),AmpSim.cv2_std(:,inx_se),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(AmpLNA.a1,cvB,'LineWidth',2)
plot(AmpLNA.a1,cv0B,'k-','LineWidth',2)
hold off
axis([3,16,0,0.4])
set(gca,'Xtick',3:2:16,'Ytick',0:0.1:0.4,'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('correlation coefficient of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameCvPdf = [Folder1,filesep,'fig1cv.pdf'];
fileNameFanoFig = [Folder1,filesep,'fig1cv.fig'];
print('-dpdf',fileNameCvPdf)
saveas(gcf,fileNameFanoFig)


%% plot Fano factor
figure(3)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
FanoB = nan(size(AmpLNA.varB,1),length(inx_se));
for i0 = 1:length(inx_se)
    FanoB(:,i0) = AmpLNA.varB(:,i0+1)./AmpLNA.meanB/OMEGA;
end
Fano0B = AmpLNA.var0B./AmpLNA.meanB/OMEGA;
% cvB = sqrt(AmpLNA.varB(:,inx_se))./AmpLNA.meanB(:,inx_se)
errorbar(AmpSim.a1(:,inx_se), AmpSim.Fano2(:,inx_se),AmpSim.Fano2_std(:,inx_se),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(AmpLNA.a1,FanoB,'LineWidth',2)
plot(AmpLNA.a1,Fano0B,'k-','LineWidth',2)
hold off
axis([3,16,1,32])
set(gca,'Xtick',3:2:16,'Ytick',[1,10,20],'YScale','log','LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('correlation coefficient of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameFanoPdf = [Folder1,filesep,'fig1Fano.pdf'];
fileNameFanoFig = [Folder1,filesep,'fig1Fano.fig'];
print('-dpdf',fileNameFanoPdf)
saveas(gcf,fileNameFanoFig)

%% plot correlation coefficient
figure(4)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
errorbar(AmpSim.a1(:,inx_se), AmpSim.corr(:,inx_se),AmpSim.corr_std(:,inx_se),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(AmpLNA.a1,AmpLNA.corr(:,inx_se+1),'LineWidth',2)
plot(AmpLNA.a1,AmpLNA.corr0,'k-','LineWidth',2)
hold off
axis([3,16,-1,0.4])
set(gca,'Xtick',3:2:16,'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('Fano factor of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameCurrPdf = [Folder1,filesep,'fig1corr.pdf'];
fileNameCurrFig = [Folder1,filesep,'fig1corr.fig'];
print('-dpdf',fileNameCurrPdf)
saveas(gcf,fileNameCurrFig)


%% plot auto correlation coefficient
figure(5)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
errorbar(AmpSim.a1(:,inx_se), AmpSim.lagAuto2(:,inx_se),AmpSim.lagAuto2_std(:,inx_se),...
    'o','MarkerSize',12,'LineWidth',0.75);
plot(AmpLNA.a1,AmpLNA.lagB(:,inx_se+1),'LineWidth',2)
plot(AmpLNA.a1,AmpLNA.lag0B,'k-','LineWidth',2)
hold off
axis([3,16,0,1])
set(gca,'Xtick',3:2:16,'YLim',[0.3,1],'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('lag one autocorrelation of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameLagPdf = [Folder1,filesep,'fig1lag.pdf'];
fileNameLagFig = [Folder1,filesep,'fig1lag.fig'];
print('-dpdf',fileNameLagPdf)
saveas(gcf,fileNameLagFig)


%% for different time scale 
inx_se_time = [1,2,4];
%plot var with different time scale
figure(6)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
h1 = errorbar(TimeSimu.a1(:,inx_se_time), TimeSimu.var2(:,inx_se_time),TimeSimu.var2_std(:,inx_se_time),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(TimeLNA.a1,TimeLNA.varB(:,[1,3,5]),'LineWidth',2)
plot(TimeLNA.a1,TimeLNA.var0B,'k-','LineWidth',2)
hold off
axis([3,16,200,2000])
set(gca,'Xtick',3:2:16,'Ytick',[200,500,2000],'YScale','log','LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('variacne of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameVarTimePdf = [Folder1,filesep,'fig1VarTime.pdf'];
fileNameFigVarTime = [Folder1,filesep,'fig1VarTime.fig'];
print('-dpdf',fileNameVarTimePdf)
saveas(gcf,fileNameFigVarTime)

%% plot cv with different time scale
figure(7)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
LNA_inx = [1,3,5];  %select certain time scale
cvB = nan(size(TimeLNA.varB,1),length(LNA_inx));
for i0 = 1:length(LNA_inx)
    cvB(:,i0) = sqrt(TimeLNA.varB(:,LNA_inx(i0)))./TimeLNA.meanB/OMEGA;
end
cv0B = sqrt(TimeLNA.var0B)./TimeLNA.meanB/OMEGA;
% cvB = sqrt(AmpLNA.varB(:,inx_se))./AmpLNA.meanB(:,inx_se)
h1 = errorbar(TimeSimu.a1(:,inx_se_time), TimeSimu.cv2(:,inx_se_time),TimeSimu.cv2_std(:,inx_se_time),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(TimeLNA.a1,cvB,'LineWidth',2)
plot(TimeLNA.a1,cv0B,'k-','LineWidth',2)
hold off
axis([3,16,0,0.4])
set(gca,'Xtick',3:2:16,'Ytick',0:0.1:0.4,'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('correlation coefficient of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameTimeCvPdf = [Folder1,filesep,'fig1cvTime.pdf'];
fileNameTimeCvFig = [Folder1,filesep,'fig1cvTime.fig'];
print('-dpdf',fileNameTimeCvPdf)
saveas(gcf,fileNameTimeCvFig)

%% plot Fano factor time scale
figure(8)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
FanoB = nan(size(TimeLNA.varB,1),length(LNA_inx));
for i0 = 1:length(LNA_inx)
    FanoB(:,i0) = TimeLNA.varB(:,LNA_inx(i0))./TimeLNA.meanB/OMEGA;
end
Fano0B = TimeLNA.var0B./TimeLNA.meanB/OMEGA;
% cvB = sqrt(AmpLNA.varB(:,inx_se))./AmpLNA.meanB(:,inx_se)
h1 = errorbar(TimeSimu.a1(:,inx_se_time), TimeSimu.Fano2(:,inx_se_time),TimeSimu.Fano2_std(:,inx_se_time),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(TimeLNA.a1,FanoB,'LineWidth',2)
plot(TimeLNA.a1,Fano0B,'k-','LineWidth',2)
hold off
axis([3,16,1,15])
set(gca,'Xtick',3:2:16,'Ytick',[1,5,10],'YScale','log','LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('correlation coefficient of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameTimeFanoPdf = [Folder1,filesep,'fig1FanoTime.pdf'];
fileNameTimeFanoFig = [Folder1,filesep,'fig1FanoTime.fig'];
print('-dpdf',fileNameTimeFanoPdf)
saveas(gcf,fileNameTimeFanoFig)

%% plot corr time scale
figure(9)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
errorbar(TimeSimu.a1(:,inx_se_time), TimeSimu.corr(:,inx_se_time),TimeSimu.corr_std(:,inx_se_time),...
    'o','MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','auto');
plot(TimeLNA.a1,TimeLNA.corr(:,[1,3,5]),'LineWidth',2)
plot(TimeLNA.a1,TimeLNA.corr0,'k-','LineWidth',2)
hold off
axis([3,16,-1,0])
set(gca,'Xtick',3:2:16,'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('Fano factor of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameTimeCurrPdf = [Folder1,filesep,'fig1corrTime.pdf'];
fileNameTimeCurrFig = [Folder1,filesep,'fig1corrTime.fig'];
print('-dpdf',fileNameTimeCurrPdf)
saveas(gcf,fileNameTimeCurrFig)


%% plot auto correlation time scale
figure(10)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
errorbar(TimeSimu.a1(:,inx_se_time), TimeSimu.lagAuto2(:,inx_se_time),TimeSimu.lagAuto2_std(:,inx_se_time),...
    'o','MarkerSize',12,'LineWidth',1,'MarkerFaceColor','auto');
plot(TimeLNA.a1,TimeLNA.lagB(:,[1,3,5]),'LineWidth',2)
plot(TimeLNA.a1,TimeLNA.lag0B,'k-','LineWidth',2)
hold off
axis([3,16,0,1])
set(gca,'Xtick',3:2:16,'YLim',[0.3,1],'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('lag one autocorrelation of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameLagTimePdf = [Folder1,filesep,'fig1lagTime.pdf'];
fileNameLagTimeFig = [Folder1,filesep,'fig1lagTime.fig'];
print('-dpdf',fileNameLagTimePdf)
saveas(gcf,fileNameLagTimeFig)

%% ==============================================================================================
% for figure 2 mutual activation motif 


%% ================================================================================================
% for figure 3 negative feedback

%% ================================================================================================
% for figure 4 system size effect
systemSizeFile = [Folder1,filesep,'toggSwiLNATheorySysSize_tau0.1_se0.1.mat'];
sysSize = load(systemSizeFile);

% size-dependent Fanofactor
figure(41)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
plot(sysSize.a1,sysSize.var0B(:,1)./sysSize.meanB/sysSize.Size(1),'k--','Linewidth',3)
for i0 = 1:length(sysSize.Size)
    plot(sysSize.a1,sysSize.varB(:,i0)./sysSize.meanB/sysSize.Size(i0),'b-','Linewidth',3)
end
axis([3,16,1,20])
set(gca,'Xtick',3:2:16,'Ytick',[1,5,10,20],'YScale','log','LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('Fano factor','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameSizeFanoPdf = [Folder1,filesep,'fig1FanoSize.pdf'];
fileNameSizeFanoFig = [Folder1,filesep,'fig1FanoSize.fig'];
print('-dpdf',fileNameSizeFanoPdf)
saveas(gcf,fileNameSizeFanoFig)

% size-dependent correlation coefficient
figure(42)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
plot(sysSize.a1,sysSize.corr0(:,1),'k--','Linewidth',3)
for i0 = 1:length(sysSize.Size)
    plot(sysSize.a1,sysSize.corr(:,i0),'b-','Linewidth',3)
end
axis([3,16,-1,0.3])
set(gca,'Xtick',3:2:16,'Ytick',-1:0.2:0.2,'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('correlation coefficient','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameSizeCorrPdf = [Folder1,filesep,'fig1CorrSize.pdf'];
fileNameSizeCorrFig = [Folder1,filesep,'fig1CorrSize.fig'];
print('-dpdf',fileNameSizeCorrPdf)
saveas(gcf,fileNameSizeCorrFig)

% size-dependent lag one auto correlation
figure(43)
set(gcf,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
hold on
plot(sysSize.a1,sysSize.lag0B(:,1),'k--','Linewidth',3)
for i0 = 1:length(sysSize.Size)
    plot(sysSize.a1,sysSize.lagB(:,i0),'b-','Linewidth',3)
end
axis([3,16,0.2,1])
set(gca,'Xtick',3:2:16,'Ytick',0.2:0.2:1,'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('lag-one autocorrelation','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
fileNameSizelagPdf = [Folder1,filesep,'fig1lagSize.pdf'];
fileNameSizelagFig = [Folder1,filesep,'fig1lagSize.fig'];
print('-dpdf',fileNameSizelagPdf)
saveas(gcf,fileNameSizelagFig)


%% ================================================================================================
% detrending of time serials
% timeSerialFileIntr = [Folder1,filesep,'toggSwiTimeSerialDetrend_N50_se0_tauc10.mat'];
% timeSerialFileExt = [Folder1,filesep,'toggSwiTimeSerialDetrend_N50_se0.1_tauc10.mat'];
timeSerialFileIntr = [Folder1,filesep,'exampleTimeDetrIntri.mat'];
timeSerialFileExt = [Folder1,filesep,'exampleTimeDetrExtri.mat'];

timeDetrIntr = load(timeSerialFileIntr);
timeDetrExtr = load(timeSerialFileExt);

%randomly plot 10 for selection
% selectedInx = randi(500,10,1);
% figure(11)
% hold on
% for i0 = 1:10;
%     plot(timeDetrExtr.time{selectedInx(i0)},timeDetrExtr.Fano2{selectedInx(i0)})
%     plot(timeDetrIntr.time{selectedInx(i0)},timeDetrIntr.Fano2{selectedInx(i0)},'k--')
%     timeDetrExtr.kcv2(selectedInx(i0))
% end
% hold off
% 
% figure(12)
% hold on
% for i0 = 1:10;
%     plot(timeDetrExtr.time{selectedInx(i0)},timeDetrExtr.corr{selectedInx(i0)})
%     plot(timeDetrIntr.time{selectedInx(i0)},timeDetrIntr.corr{selectedInx(i0)},'k--')
% end
% hold off
% 
% figure(13)
% hold on
% for i0 = 1:10;
%     plot(timeDetrExtr.time{selectedInx(i0)},timeDetrExtr.lag2{selectedInx(i0)})
%     plot(timeDetrIntr.time{selectedInx(i0)},timeDetrIntr.lag2{selectedInx(i0)},'k--')
% end
% hold off

exampleInx = 1;
figure(51)
set(gcf,'Units','inches','Position',[0 0 8 4],...
'PaperPositionMode','auto');
hold on
plot(timeDetrExtr.time{exampleInx},timeDetrExtr.Fano2{exampleInx},'r-','LineWidth',2)
plot(timeDetrIntr.time{exampleInx},timeDetrIntr.Fano2{exampleInx},'k-','LineWidth',2)
set(gca,'xlim',[0,800],'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('time','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('Fano factor of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
exampleTimeFanoPdf = [Folder1,filesep,'fig5TimeFano.pdf'];
exampleTimeFanoFig = [Folder1,filesep,'fig5TimeFano.fig'];
print('-dpdf',exampleTimeFanoPdf)
saveas(gcf,exampleTimeFanoFig)
    
figure(52)
set(gcf,'Units','inches','Position',[0 0 8 4],...
'PaperPositionMode','auto');
hold on
plot(timeDetrExtr.time{exampleInx},timeDetrExtr.corr{exampleInx},'r-','LineWidth',2)
plot(timeDetrIntr.time{exampleInx},timeDetrIntr.corr{exampleInx},'k-','LineWidth',2)
set(gca,'xlim',[0,800],'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('time','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('correlation coefficient','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
exampleTimeCorrPdf = [Folder1,filesep,'fig5TimeCorr.pdf'];
exampleTimeCorrFig = [Folder1,filesep,'fig5TimeCorr.fig'];
print('-dpdf',exampleTimeCorrPdf)
saveas(gcf,exampleTimeCorrFig)

figure(53)
set(gcf,'Units','inches','Position',[0 0 8 4],...
'PaperPositionMode','auto');
hold on
plot(timeDetrExtr.time{exampleInx},timeDetrExtr.lag2{exampleInx},'r-','LineWidth',2)
plot(timeDetrIntr.time{exampleInx},timeDetrIntr.lag2{exampleInx},'k-','LineWidth',2)
set(gca,'xlim',[0,800],'LineWidth',1.5,'FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('time','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
exampleTimeLag2Pdf = [Folder1,filesep,'fig5TimeLag2.pdf'];
exampleTimeLag2Fig = [Folder1,filesep,'fig5TimeLag2.fig'];
print('-dpdf',exampleTimeLag2Pdf)
saveas(gcf,exampleTimeLag2Fig)