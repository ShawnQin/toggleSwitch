% this program replot the bifurcation diagram of the double positive feed
% back loop. The imported data are two .csv file in this folder
% default output are two bifurcation diagram

clear
clc

% set the folder to save figures
figFolder = '/Users/shan/Documents/MATLAB/criticlity/ExtrinsicNoiseDataFigure/mutualActiv';

%load the data
rawB = csvread('mutualActiveBifurB.csv');
rawA = csvread('mutualActiveBifurA.csv');
dataA = rawA; % set 0 to nan
dataA(dataA == 0) = nan;
dataB = rawB; % set 0 to nan
dataB(dataB == 0) = nan;


%% plot the bifurcation diagram of A
fh1 = figure(1);  % this is the bifurcation of A

% set the figure size to mathc previous figures
set(fh1,'Units','inches','Position',[0 0 5.7 4.3],...
'PaperPositionMode','auto');
plot(dataA(:,1),dataA(:,2:4),'LineWidth',3)
set(gca,'XLim',[0.1,0.7],'XTick',0.1:0.1:0.7,'YTick',0:0.2:1,...
    'LineWidth',1,'FontSize',20)
xlabel('K1','FontSize',24)
ylabel('Expression of A','FontSize',24)

saveFileAfig = fullfile(figFolder,'bifurA.fig');
saveFileAPdf = fullfile(figFolder,'bifurA.pdf');
saveas(fh1,saveFileAfig)
fh1.PaperUnits = 'inches';
print('-painters','-dpdf',saveFileAPdf)

%% plot the bifurcation diagram of B
fh2 = figure(2);  % this is the bifurcation of A

% set the figure size to mathc previous figures
set(fh2,'Units','inches','Position',[0 0 5.7 4.3],...
'PaperPositionMode','auto');
plot(dataB(:,1),dataB(:,2:4),'LineWidth',3)
set(gca,'XLim',[0.1,0.6],'XTick',0.1:0.1:0.6,'YTick',0:0.2:1,...
    'LineWidth',1,'FontSize',20)
xlabel('K1','FontSize',24)
ylabel('Expression of B','FontSize',24)

saveFileBfig = fullfile(figFolder,'bifurB.fig');
saveFileBPdf = fullfile(figFolder,'bifurB.pdf');
saveas(fh2,saveFileBfig)
fh2.PaperUnits = 'inches';
print('-painters','-dpdf',saveFileBPdf)
