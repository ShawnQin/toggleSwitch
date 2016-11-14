% this program replot all the figures in manuscripts
% revised on 11/13/2016

clear
clc


%load data for figure 1
Folder1 = '/Users/shan/Documents/MATLAB/criticlity/ExtrinsicNoiseDataFigure/figure1';
filenameAmpSimu =[Folder1, filesep,'ExtriFig1SimuAmp.mat'];
filenameTimeSimu = [Folder1, filesep,'ExtriFig1SimuTime.mat'];
filenameAmpLNA = [Folder1, filesep,'toggSwiLNATheoryAmpl_N20.mat'];
filenameTimeLNA = [Folder1, filesep,'toggSwiLNATheoryAmpl_N20.mat'];

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

%plot the variance of B

figure('Units','inches','Position',[0 0 8 7],...
'PaperPositionMode','auto');
figure(1)
hold on
h1 = errorbar(AmpSim.a1(:,inx_se), AmpSim.var2(:,inx_se),AmpSim.var2_std(:,inx_se),...
    'o','MarkerSize',15,'LineWidth',2);
plot(AmpLNA.a1,AmpLNA.varB(:,inx_se+1),'LineWidth',2)
plot(AmpLNA.a1,AmpLNA.var0B,'k-','LineWidth',2)
hold off
axis([3,16,100,8000])
set(gca,'Xtick',3:2:16,'YScale','log','FontName','Times','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Helvetica','FontWeight','Bold')
ylabel('variacne of B','FontSize',30,'FontName','Helvetica','FontWeight','Bold')

% plot cv
