%this program use MatCont to do bifurcation analysis and plot nice diagram
% A insect outbreak model subject to extrinsic noise
% the fixed point with EN is defined as the extrema of quasistationary
% probability distribution through FPE
% last revised on 12/24/2016

close all
clear all
clc

global cds sys

%%%%% Set continuation pause environment variables %%%%%
%%
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

%%%%% Set system %%%%%
syshandle=@toggSwiMatcont;  %Specify system file


%%
%%%%% ODE Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Define a few intermediate functions for ODE integration %%%%%
SubFunHandles=feval(syshandle);  %Get function handles from system file
RHShandle=SubFunHandles{2};      %Get function handle for ODE

%Set ODE parameter
a1 = 2;
a2 = 10;
a0 = 0.03;
n1 = 2;
n2 = 2;

%Set ODE initial condition
xinit=[0.5;10];
tspan = [0 1e3];

%%%%%% Define an anynomous function to pass to integrator %%%%%
RHS_no_param=@(t,x)RHShandle(t,x,a1,a2,a0,n1,n2); 

%%%%% Set ODE integrator options %%%%%
options=odeset;
options=odeset(options,'RelTol',1e-5);
options=odeset(options,'maxstep',1e-1);

%%%%%% Integrate until a steady state is found. %%%%%
[tout xout]=ode15s(RHS_no_param,tspan,xinit,options);

figure()
plot(tout,xout,'-k','linewidth',2)
%axis([0 100 -1 1])

%%
%%%%% Continuation from equilibrium %%%%%%%%%%%%%%%%

%%%%%% Set initial condition as the endpoint of integration.  Use
%%%%%% to bootstrap the continuation.
%xinit=xout(length(xout));
xinit=xout(end,:)';

pvec=[a1,a2,a0,n1,n2]'; % Initialize parameter vector

ap=1; % Denote 'active' parameter for continuation, number indicate which parameter

%%%%% Initialize continuation %%%%%
[x0,v0]=init_EP_EP(syshandle, xinit, pvec, ap); %Initialize equilibrium

%%%%% Initialize Matcont options %%%%%
opt=contset;
opt=contset(opt,'MaxNumPoints',1000); %Set numeber of continuation steps
opt=contset(opt,'MaxStepsize',.05);  %Set max step size
opt=contset(opt,'Singularities',1);  %Monitor singularities
opt=contset(opt,'Eigenvalues',1);    %Output eigenvalues 
opt=contset(opt,'InitStepsize',0.01); %Set Initial stepsize

%%%%% Continuation %%%%%
%Equilibrium continuation, the output contain all the information about the
%fixed point: 
% x1   the equilibrium value and corresponding parameter
% v1   the tangent vectors along the curve, both x1 and v1 define the curve
% s1   an array of structure contain the information about the sigularities
% s1.index index of the sigularity point in x
[x1,v1,s1,h1,f1]=cont(@equilibrium,x0,v0,opt); 


figure()
cpl(x1,v1,s1,[3 2]); %plot the continuation from equilibrium forward
%%
%%%%% Continuation from equilibrium backward %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x0,v0]=init_EP_EP(syshandle, xinit, pvec, ap); %Initialize equilibrium
opt=contset(opt,'Backward',1);
[x2,v2,s2,h2,f2]=cont(@equilibrium,x0,v0,opt);

figure()
cpl(x2,v2,s2,[3 2]);
%%
%%%%% Plotting script.  x=continuation info  f=eigenvalues %%%%%%%%%
%%%%% Extract eigenvalues and match with curves for stability %%%%%%

figure()

xeqcurve=x1;

%%%%% This is the last eigenvalue.  
%%%%% These are ordered smallest to largest (real part).
%%%%% So the last one determines stability.
%now this is a 3-variables system,which has three eigenvaluse,
%the stability is determined by the largest one
minevaleq=max(real(f1),[],1);  

L=length(xeqcurve(1,:)); % total number of points

curveind=1;
lengthind=0;
maxlengthind=0;
evalstart=floor(heaviside(minevaleq(1))); %stability type of first point
datamateq=zeros(4,L);  %data matrix for equilibrium points, stable and unstable each two rows
eigDataeq = zeros(4,L); % this matrix store the inforamtion of the largest eigen value of each branch

for i=1:L
    evalind=floor(heaviside(minevaleq(i))); %eigen value, dertermines the stability
    if evalstart~=evalind
        curveind=curveind+1;
        i;
        evalstart=evalind; %reset the current stability type
        maxlengthind=max(lengthind,maxlengthind);   %length of the "branch"
        lengthind=0;
    end
    datamateq(1,i) = xeqcurve(3,i); % This is the parameter that is varied.
    datamateq(2,i) = xeqcurve(2,i); % This is the dependent axis of the bifurcation plot.  The one you wish to plot
    datamateq(3,i) = evalind;   % stability type,0 unstable,1 stable
    datamateq(4,i) = curveind;  %branch of curves type
    
    eigDataeq(1,i) = xeqcurve(3,i);  %control parameter
    eigDataeq(2,i) = minevaleq(i);   % largest eigenvalue
    eigDataeq(3,i) = evalind;        % stability index
    eigDataeq(4,i) = curveind;       % branches index
    lengthind=lengthind+1;
end

maxlengthind=max(maxlengthind,lengthind);

curveindeq=curveind;

for i=1:curveindeq   %cycle trough different "branches"
    index=find(datamateq(4,:)==i);
    eval(['curve' num2str(i) 'eq' '=datamateq(1:3,index);']); %different branches
    eval(['eigen' num2str(i) 'eq' '=eigDataeq(1:3,index);']); %different branches
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the parameter of the figure
figure()
figureSize = [0,0,8,6];
set(gcf,'Units','inches','Position',figureSize,'PaperPositionMode','auto')
axis([4,18,0,10])
LineWidth = 3;
LabelFontSize = 36;
LabelFontName = 'Helvetica';
AxisFontSize = 30;
ticketWidth = 1.5;

for i=1:curveindeq
    stability=eval(['curve' num2str(i) 'eq(3,1)']);
    if stability==0
        plotsty='-';
    else
        plotsty='--';
    end
    
    plotcolor='k';

    plotstr=strcat(plotcolor,plotsty);
    
    plot(eval(['curve' num2str(i) 'eq(1,:)']),eval(['curve' num2str(i) 'eq(2,:)']),...
        plotstr','Linewidth',LineWidth)
    hold on
end

%%%%% Final adjustments to make things look nicer %%%%%
xlabel('a1','Fontsize',LabelFontSize,'FontName',LabelFontName)
ylabel('level of B','Fontsize',LabelFontSize,'FontName',LabelFontName)
%title('saddle node bifurcation','fontsize',fontsizevar)
%set(gcf,'units','inches','pos',[0 0 3.1*scalefactor 2.61*scalefactor])
xlim([4,18])
ylim([0,10])
set(gca,'Xtick',4:2:18,'Ytick',0:2:10,'LineWidth',ticketWidth,'FontName',...
    LabelFontName,'FontSize',AxisFontSize)

%save the figure
cf = pwd; %current folder
fileNameBifurFig = [cf,filesep,'figure and data',filesep,'bifur_B','.fig'];
fileNameBifurPdf = [cf,filesep,'figure and data',filesep,'bifur_B','.pdf'];
print('-dpdf',fileNameBifurPdf)
saveas(gcf,fileNameBifurFig)
 
%%
%plot the change of eigenvalue
figure()
figureSize = [0,0,8,6];
set(gcf,'Units','inches','Position',figureSize,'PaperPositionMode','auto')
% axis([0.35,0.65,-2,1])
LineWidth = 3;
LabelFontSize = 36;
LabelFontName = 'Helvetica';
AxisFontSize = 30;
ticketWidth = 1.5;
 
for i=1:curveindeq   %only plot the stable branches
    stability=eval(['eigen' num2str(i) 'eq(3,1)']);
    if stability==0
        plotsty='-';
    else
        plotsty='--';
    end
    
    plotcolor='k';

    plotstr=strcat(plotcolor,plotsty);
    
    plot(eval(['eigen' num2str(i) 'eq(1,:)']),eval(['eigen' num2str(i) 'eq(2,:)']),...
        plotstr','Linewidth',LineWidth)
    hold on
end

xlabel('a1','Fontsize',LabelFontSize,'FontName',LabelFontName)
ylabel('\lambda_{max}','Fontsize',LabelFontSize,'FontName',LabelFontName)
xlim([4,18])
ylim([-1,0.6])
set(gca,'Xtick',4:2:18,'Ytick',-1:0.5:0.5,'LineWidth',ticketWidth,'FontName',...
    LabelFontName,'FontSize',AxisFontSize)
fileNameEigenFig = [cf,filesep,'figure and data',filesep,'eigen','.fig'];
fileNameEigenPdf = [cf,filesep,'figure and data',filesep,'eigen','.pdf'];
print('-dpdf',fileNameEigenFig)
saveas(gcf,fileNameEigenPdf)