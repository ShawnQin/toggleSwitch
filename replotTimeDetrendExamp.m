% this program replot the last figure in manuscript, so that its size match
% the requirement of PRE and PRL with width 3.4 and height 2.8
% this also faciliate the polish process in AI
% last revised on 07/01/2017

clc
clear

%% load the original fiure, default in 
figFolder = '/Users/shan/Documents/MATLAB/criticlity/ExtrinsicNoiseDataFigure';
subFolder = 'figure4Detrend';
figResPath = fullfile(figFolder,subFolder,'togg-detrend-residule-example-01132017.fig');
expResFig = openfig(figResPath);
expResFig.WindowStyle = 'normal';
expResFig.DockControls = 'off';
expResFig.PaperPositionMode = 'auto';


%% set the parameter of graphics
figureSize = [0,0,3.4,2.5];
axisFontSize = 14;
tickWidth = 1;
labelFontSize = 14;
colorSet = [3,110,184;224,135,51;202,39,44;0,0,0]/256;

%% reset the figure
%set(expFig,'Units','inches','Position',figureSize,...
%    'PaperPositionMode', 'auto')
expResFig.Units = 'inches';
expResFig.Position = figureSize;
expResFig.CurrentAxes.XLim = [0,800];
expResFig.CurrentAxes.XTick = 0:200:800;
set(expResFig.CurrentAxes,'LineWidth',tickWidth,'FontName','Helvetica',...
    'FontSize',axisFontSize)
xlabel('Time','FontName','Helvetica','FontSize',labelFontSize)
ylabel('Residule','FontName','Helvetica','FontSize',labelFontSize)
%set(gca,'XLabel','Time','YLabel','Residule','FontSize',labelFontSize)
figNamefig = fullfile(figFolder,subFolder,'timeDetrExpfprPRE.fig');
figNamepdf = fullfile(figFolder,subFolder,'timeDetrExpfprPRE.pdf');
saveas(expResFig,figNamefig)
expResFig.PaperUnits = 'inches';
%expFig.PaperPosition = figureSize;
print('-painters','-dpdf',figNamepdf)

%% load the example IN trajectory
INTrajFile = fullfile(figFolder,subFolder,'togg_detrend_IN-example-01132017.fig');
INTrajFig = openfig(INTrajFile);
INTrajFig.WindowStyle = 'normal';
INTrajFig.DockControls = 'off';
INTrajFig.PaperPositionMode = 'auto';

INTrajFig.Units = 'inches';
INTrajFig.Position = figureSize;
INTrajFig.CurrentAxes.XTick = [];
INTrajFig.CurrentAxes.YTick = 0:250:1000;
set(INTrajFig.CurrentAxes,'LineWidth',tickWidth,'FontName','Helvetica',...
    'FontSize',axisFontSize)
%xlabel('Time','FontName','Helvetica','FontSize',labelFontSize)
ylabel('Residule','FontName','Helvetica','FontSize',labelFontSize)
%set(gca,'XLabel','Time','YLabel','Residule','FontSize',labelFontSize)
figNameTrajFig = fullfile(figFolder,subFolder,'trajINexpPRE.fig');
figNameTrajPdf = fullfile(figFolder,subFolder,'trajINexpPRE.pdf');
saveas(INTrajFig,figNameTrajFig)
INTrajFig.PaperUnits = 'inches';
%expFig.PaperPosition = figureSize;
%print(INTrajFig,figNameTrajPdf,'-painters','-dpdf')
print('-painters','-dpdf',figNameTrajPdf)


%% load example Fano factor
FanoExpFile = fullfile(figFolder,subFolder,'fig5TimeFano_example.fig');
FanoFig = openfig(FanoExpFile);
FanoFig.WindowStyle = 'normal';
FanoFig.DockControls = 'off';
FanoFig.PaperPositionMode = 'auto';

FanoFig.Units = 'inches';
FanoFig.Position = figureSize;
set(FanoFig.CurrentAxes,'LineWidth',tickWidth,'FontName','Helvetica',...
    'FontSize',axisFontSize,'FontWeight','normal')
xlabel('Time','FontName','Helvetica','FontSize',labelFontSize,'FontWeight','normal')
ylabel('Fano factor of B','FontName','Helvetica','FontSize',labelFontSize)

figNameFanofig = fullfile(figFolder,subFolder,'FanoExpPRE.fig');
figNameFanoPdf = fullfile(figFolder,subFolder,'FanoExpPRE.pdf');
saveas(FanoFig,figNameFanofig)
FanoFig.PaperUnits = 'inches';
print('-painters','-dpdf',figNameFanoPdf)