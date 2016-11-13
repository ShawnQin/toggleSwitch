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
figure(1)
hold on
errorbar(AmpSim.a1(:,inx_se), AmpSim.var2(:,inx_se),AmpSim.var2_std(:,inx_se),'o','MarkerSize',15)
plot(AmpLNA.a1,AmpLNA.varB(:,inx_se+1),'LineWidth',3)
plot(AmpLNA.a1,AmpLNA.var0B,'k-','LineWidth',3)
hold off
set(gca,'YScale','log','FontSize',24,'FontWeight','Bold')
xlabel('a1','FontSize',30,'FontName','Myriad Pro','FontWeight','Bold')
ylabel('variacne of B','FontSize',30,'FontName','Myriad Pro','FontWeight','Bold')

% plot cv
