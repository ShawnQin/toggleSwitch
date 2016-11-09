%**********************************************************************
% this grogram repeat the work done be Micheal Assaf in 2013, published on
% PRL. This is a very simple auto activation gene circuit
% The authors discussed the effect of extrinsic noise in detail 
% Here, I simluate this model using modiied Gillespie algorithm
% to show how extrinsic nose may affect the distribution, bistable region
% and even te auto correlation function 
% last revised on Dec. 20,2015
% ------------------------------------------------------------------
% the main function calculate the variance and lag one autocorrelation for
% EN strength. The distribution and autocorrelation of numerical methods
% the swithing time between high to low states and low to high states
%************************************************************************
function selfactivation()
 
%parameters that will be used in all the subroutines
% syms x a0 x0
% z= diff(-a0-(1-a0)*x^2/(x^2+x0^2)+x,x);
% % FPE_inxtrinsic()
a0 = 0.05;
OMIGA = 300;  %system size
x0list = [0.3:0.02:0.6,0.65,0.7];
x0 = 0.43;
tauc = 10;
NUM = 5e3;
N = 10;
runTime = 1e3;
sigmae = 0.1;
param = [a0,x0,OMIGA,tauc,sigmae];
%__________________________________________________________________________
%simulate the distribution
% slefActiDist(param,NUM,x0list);

%__________________________________________________________________________


fileName = '/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData/stochSteady_tau_1000_0.1.mat';
numericialSelfPromCLEExtrinsic(param,runTime,fileName,N);
%__________________________________________________________________________
%calcuate variance and autocorrelation use CLE 
x0list = 0.3:0.005:0.52;
CLE_noiseAmpli(param,x0list)

%__________________________________________________________________________
%calculate variance use WKB at adia
x0 = 0.3:0.005:0.52;
for i = 1:length(x0);    
    [var_in,var_ex] = var_WKB(a0,x0(i),300,0.1);
    VAR(i,:) =[var_in,var_ex];
end
plot(x0',VAR)

%__________________________________________________________________________

% numeric method to calculate variance based on distribution
% fileName = 'self-activ-distri-1e-2-1221.csv';
% N = 20;
% numeric_var_atuo(param,x0list,N,fileName);
% [x0,t_swi] = switchTime(param,x0list,fileName,N);

%__________________________________________________________________________

%numeric method of switching time
% N = 100;
% fileName = 'stochSteady_tau_10_0.1.mat';
fileName = '/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData/stochSteady_tau_1000_0.1.mat';
% [x0,t_swi]=switchTime(param,x0list,fileName,N);

%__________________________________________________________________________

%numeric result of swithing time
% plotSwiTime()

%__________________________________________________________________________

%WKB method with only intrinsic noise
% WKBswiRate(a0,OMIGA)

%__________________________________________________________________________

%bifurcation diagram
% bifurcation(a0,OMIGA)
%__________________________________________________________________________
%plot the distribution and calcualate the variacne from numeric simulation
plotDistri(x0list)
%__________________________________________________________________________

end
function dydt = odeSelfActi(t,y,a0,OMIGA,x0)
        dydt = a0*OMIGA + (1-a0)*OMIGA*y^2/(y^2+ (OMIGA*x0)^2) - y;
end
function [time,Y]=YGlpExtri(inV,OMIGA,a0,x0,totTime,sigmae,tauc)
%ths is a modified Gillespie algorithm
%here the sigame is not the same as in the paper

%first, produce all the extrinsic noise\
dt = 0.1;
Next = round(1.1*totTime/dt);
mu = exp(-dt/tauc);
ss = sqrt(2*sigmae^2/tauc);  %revised on Jan 4th,2016

Dext = ss*sqrt((1-mu^2)*tauc/2);
Xtemp = 0;
Xe = zeros(Next,1);
for i1 = 1:Next;
    Xtemp = Xtemp*mu+ Dext*randn;
    Xe(i1) = Xtemp;
end

Y(1,1) = inV;
Y(1,2) = 0;  % the second column store extrinsic noise
S = [1,-1];
deg = 1;
t = 0;
te = dt;
i = 1;
j = 0;
while(t < totTime)
    a = [OMIGA*(a0 + (1-a0)*Y(i,1)^2/(Y(i,1)^2+(OMIGA*x0)^2)); deg*Y(i,1)];
    [tau,inx] = min(-1./a.*log(rand(2,1)));
    if(t+tau<te)
        Y(i+1,1) = Y(i,1) + S(:,inx);
        Y(i+1,2) = deg-1;
        t = t + tau;  
        time(i+1)= t;
        i = i+1;
    else
        deg = 1+Xe(j+1);
        if(deg<0)
            deg = 0;%negtive degradation is impossible
        end
        t = te;
        te = te+dt;   %updata the extrinsic time point
        j = j+1;
    end
    
end
time = [0;time(1:end-1)'];
end
function [time,Y] = CLEselfActi(inV,param,totTime)

%parameters
a0 = param(1);
x0 = param(2);
OMIGA = param(3);
tauc = param(4);
s = param(5);  %revised on Dec 27,2015

sigmae = sqrt(2*s^2/tauc);

dt = 0.05;
N = round(totTime/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = [inV;0];
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);
for i1 = 1:N;
    determin = OMIGA*(a0 + (1-a0)*Xtemp(1)^2/(Xtemp(1)^2+(x0*OMIGA)^2))-(1+Xtemp(2))*Xtemp(1);
    stoch = [sqrt(abs(OMIGA*(a0 + (1-a0)*Xtemp(1)^2/(Xtemp(1)^2+(x0*OMIGA)^2)))),-sqrt(abs((1+Xtemp(2))*Xtemp(1)))];
    Xtemp(1) = Xtemp(1) + determin*dt + stoch*randn(2,1)*sqrt(dt);
    Xtemp(2) = Xtemp(2)*mu+ Dext*randn;
    if(Xtemp(2)<-1)
        Xtemp(2) = -1;
    end
    X(i1,:) = Xtemp;
end

time = (0:dt:N*dt)';
Y = [[inV;0]';X];
% plot(time,Y,'LineWidth',2)
% legend('gene1','gene2')
% xlabel('time','FontSize',24,'FontName','calibri')
% ylabel('gene expression','FontSize',24,'FontName','calibri')

end
function t_swi = CLEswitch_ext(inV,param,range)
%range define the vicinity of a state
a0 = param(1);
x0 = param(2);
OMIGA = param(3);
tauc = param(4);
s = param(5);

sigmae = sqrt(2*s^2/tauc);  %revised on Dec 27

dt = 0.1;
Xtemp = [inV;0];
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);
TMAX = 1e6; %maximun running time
    time = 0;
    while(time < TMAX)
        determin = OMIGA*(a0 + (1-a0)*Xtemp(1)^2/(Xtemp(1)^2+(x0*OMIGA)^2))-(1+Xtemp(2))*Xtemp(1);
        stoch = [sqrt(abs(OMIGA*(a0 + (1-a0)*Xtemp(1)^2/(Xtemp(1)^2+(x0*OMIGA)^2)))),-sqrt(abs((1+Xtemp(2))*Xtemp(1)))];
        Xtemp(1) = Xtemp(1) + determin*dt + stoch*randn(2,1)*sqrt(dt);
        Xtemp(2) = Xtemp(2)*mu+ Dext*randn;
        if(Xtemp(2)<-1)
            Xtemp(2) = -1;
        end
        time = time + dt;
        if(Xtemp(1)<range(2) && Xtemp(1)>range(1))
            t_swi = time;
            return
        end
    end
    t_swi = nan;   %if no swithing happen, return nan

end

function slefActiDist(param,NUM,x0list)
%this program find the numeric metastable states by locating the max mum and minimun
%of ditribution curve
%NUM is the total time of simulations

allIni = lhsdesign(NUM,1);
t_sim = 2000;
N_s = 1e3;    % the transient time which the system has not relax to steady state
Data  = zeros(NUM*N_s,length(x0list));
for j0 = 1:length(x0list);
    param(2) = x0list(j0);
    for i0 = 1:NUM;
        [time,Y]=CLEselfActi(allIni(i0),param,t_sim);
        Data((i0-1)*N_s+1:i0*N_s,j0) = Y(end-N_s+1:end,1);
    end

end
    save('selfDistri-ne-1228.mat','Data')
end
function [A,Ae] = jacMatrx(param,vec)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
syms y E a0 x0 OMIGA tauc 
J1 = jacobian(a0*OMIGA + (1-a0)*OMIGA*y^2/(y^2+ (OMIGA*x0)^2) -y,y);
A = double(subs(J1,{y,a0,x0,OMIGA},{vec,param(1),param(2),param(3)}));
% J2 = jacobian([a0*OMIGA+(1-a0)*OMIGA*y^2/(y^2+(OMIGA*x0)^2)-(1+E)*y;-1/tauc*E],[y E]);
J2 = jacobian([a0*OMIGA+(1+E)*(1-a0)*OMIGA*y^2/(y^2+(OMIGA*x0)^2)-y;-1/tauc*E],[y E]);
Ae = double(subs(J2,{y,E,a0,x0,OMIGA,tauc},{vec,0,param(1),param(2),param(3),param(4)}));

end
function CLE_noiseAmpli(param,x0list)
%this function comparing the correlation variance, and lag one
%auto correlation for CLE with extrinsic noise and without 
% sigmae define diffent extrinsic noise amplitude

a0 = param(1);
OMIGA = param(3);
x0 = param(2);
tauc = param(4);
se = [1e-3,1e-2,0.1,0.3];
sigmae = sqrt(se.^2*2/tauc);

% sigmae = [1e-3,0.1];
tspan = [0 1000];
y0 = 300;
[t,y] = ode15s(@odeSelfActi,tspan,y0,[],a0,OMIGA,x0);

variance_theory = zeros(length(x0list),1); %only intrinsic noise
variance_ext_cle_theory = zeros(length(x0list),length(sigmae));
lagOneLe_theory = zeros(length(x0list),1);
lagOneCLEext_theory = zeros(length(x0list),length(sigmae));

for j1 = 1:length(x0list);
    param(2)= x0list(j1);
    [t,y] = ode15s(@odeSelfActi,tspan,y0,[],a0,OMIGA,x0list(j1));
    steadyODE(j1) = y(end);
    for j2 = 1:length(sigmae)
        [A,Ae] = jacMatrx(param,steadyODE(j1));

        Dc = a0*OMIGA + (1-a0)*OMIGA*steadyODE(j1)^2/(steadyODE(j1)^2+ (OMIGA*x0list(j1))^2) + steadyODE(j1);
        Dcle = [a0*OMIGA + (1-a0)*OMIGA*steadyODE(j1)^2/(steadyODE(j1)^2+ (OMIGA*x0list(j1))^2) + steadyODE(j1),0;...
             0,sigmae(j2)^2];
 
        C1 = Dc/2/abs(A);
        variance_theory(j1) = C1;
        lagOneLe_theory(j1) = expm(A);
        
        C2 = lyap(Ae,Dcle);
%         corrcoef_ext_cle_theory(j1,j2) = C2(1,2)/(sqrt(C2(1,1)*C2(2,2)));
        variance_ext_cle_theory(j1,j2) = C2(1,1);
        G2 = C2*expm(Ae');
        lagOneCLEext_theory(j1,j2) = G2(1,1)/C2(1,1);

    end
end

%plot the theoretic results
colorSet = [ones(4,1)*60/256,[0;0.3;0.7;1],ones(4,1)*189/256];

%variance
figure(1)
plot(x0list,variance_theory,'r-','Linewidth',3)
hold on
for j4 = 1:length(sigmae);
    plot(x0list,variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--') 
end
hold off
legend('CLE','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('x0','FontSize',24,'FontWeight','Bold')
ylabel('variance','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

%auto correlation
figure(2)
plot(x0list,lagOneLe_theory,'r-','LineWidth',3)
hold on
for j5 = 1:length(sigmae);
    plot(x0list,lagOneCLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:),'LineStyle','--')    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('x0','FontSize',24,'FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
end

function [x0,t_swi]=switchTime(param,x0list,fileName,N)
%this fucntion calcuate the noise induced swith time 
%hishstates are the QSD mean of high states
%lowstates are the QSD mean of lower states

%compare three different situation: without EN, with EN but different amplitute
% threeState= csvread(fileName);
load(fileName); %modified on Dec 28,2015
inx_h = ~isnan(steady(:,2));
high = steady(inx_h,1);
low = steady(inx_h,3);
min_saddle = min(steady(inx_h,2));
x0 = x0list(inx_h);

%____________________________________________________________
%for intrinsic noise
% x0 = (0.47:0.01:0.51);
% for i = 1:length(x0);
%     y0 = steady(param(1),x0(i));
%     high(i) = param(3)*max(y0);
%     low(i) = param(3)*min(y0);
%     saddle(i) = param(3)*y0(2);
% end
% min_saddle = min(saddle);
%__________________________________________________________
t_swi = zeros(N,length(high));
MAXTRY = 2*N;   % maximum trying times
for i0 = 1:length(high); 
    param(2) = x0(i0);
   
    % firt from high to low state
    lowRef = [low(i0) - 20,min([low(i0) + 20,min_saddle])]; %the border of lower state
    n = 0;
    flag = 0;
    while(n<N && flag<MAXTRY)
        ts = CLEswitch_ext(high(i0),param,lowRef);
        if(~isnan(ts))
            n = n+1;
            t_swi(n,i0) = ts;
        end
            flag = flag + 1;
    end
end
saveFile = ['selfSwiTime_tau',num2str(param(4)),'_se',num2str(param(5)),'_',date,'.mat'];
save(saveFile,'t_swi','x0')
end
function plotDistri(x0list)

%first find all the minimum and maximum of distribution
%then calculate the variance at each metastate

%      load('/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData/selfDistri-1220.mat')
    load('selfDistri_tau1e3_ne-1228'); 
%     load('selfDistri-tau10_se0_1.mat');
    steady = nan(length(x0list),3);
    for i0 = 1:length(x0list);
        [f,xi]=ksdensity(Data(:,i0),'npoints',400,'support',[-4 400]);
        [maxTab, minTab] = peakdet(f,1e-4,xi);
        semilogy(xi,f)
        hold on
        plot(maxTab(:,1),maxTab(:,2),'ro')

        if (size(maxTab,1)==1&& maxTab(1,1) > 50)
            steady(i0,1) = maxTab(1,1);
        elseif(size(maxTab,1)==1&& maxTab(1,1) < 50)
            steady(i0,3) = maxTab(1,1);
        elseif(size(maxTab,1)==2)
            steady(i0,1) = max(maxTab(:,1));
            steady(i0,3) = min(maxTab(:,1));
        end
        if(~isempty(minTab))
            plot(minTab(:,1),minTab(:,2),'g*')
            steady(i0,2) = minTab(:,1);
        end
        
        %theoretical
%         area(xi,f)
%         hold on
    end
%     save('stochSteady_tau_10_0.1.mat','x0list','steady')    figure(2)
    plot(x0list',steady,'ro')
    
%     threeState= csvread('self-activ-distri-1e-4-1221.csv');
% %     threeState= csvread('self-activation-distribution-1120.csv');
%     plot(threeState(:,1),threeState(:,2:4),'o')
    inx_h = ~isnan(steady(:,2));
    high = steady( ~isnan(steady(:,1)),2);  %all the high state
    saddle = steady(:,2);
    saddle(isnan(saddle))=0;
    var_high = nan(length(x0list),1);
    for j0 = 1:length(high);
        var_high(j0) = var(Data(Data(:,j0)>saddle(j0),j0));
    end
    figure(3)
    plot(x0list',var_high,'r--o')
    xlabel('x0','FontSize',24,'FontWeight','Bold')
    ylabel('variance at high state','FontSize',24,'FontWeight','Bold')
    set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',2)
    legend('simu-\sigma_{E}=0.1,\tau_{c} = 10^{3}','simu-\sigma_{E}=0.1,\tau_{c} = 10','simu-\sigma_{E}=0','WKB-\sigma_{E}=0.1,\tau_{c} = 10^{3}','WKB-\sigma_{E}=0')
end
function numeric_var_atuo(param,x0list,N,fileName)
%this funtion calculate the varinace auto correlation numerically

%first,load the file which contain the points of parameters for bistable
%region. a csv file
threeState= csvread(fileName);
allHL = threeState(:,[2,4]);
inx_h = find(~isnan(threeState(:,3)));
variCLEext = zeros(size(threeState,1),N);
detr_variCLEext = zeros(size(threeState,1),N);
variCLE = zeros(size(threeState,1),N);
lagCLEext = zeros(size(threeState,1),N);
detr_lagCLEext = zeros(size(threeState,1),N);
lagCLE = zeros(size(threeState,1),N);

dt = 0.05;
totTime = 400;

%with extrinsic noise
for i0 = 1:size(threeState,1);
    j0 = 0;
    param(2) = x0list(i0);
    inv = max(allHL(i0,:));
    while(j0 <=N )
        [t_CLEext,Y_CLEext] = CLEselfActi(inv,param,totTime);
        if(i0>=min(inx_h) && i0<=max(inx_h)) % b we need selecet data mannually
            plot(t_CLEext,Y_CLEext)
            choice = input('the choice:\n','s');      %input the choice
            if(strcmp(choice,'m')) %mannually selection
                [x,~] = ginput(2);
                [~,inx1] = sort(abs(t_CLEext-x(1)));
                [~,inx2] = sort(abs(t_CLEext-x(2)));
                dataSelect = Y_CLEext(inx1(1):inx2(1),1);   %select data
                t_s = t_CLEext(inx1(1):inx2(1));
                detrendData = dataDetrend(dataSelect,t_s,100);  %detrending methods
            elseif(strcmp(choice,'y'))  %use all the data
                dataSelect = Y_CLEext(:,1);
                t_s = t_CLEext;
                detrendData = dataDetrend(dataSelect,t_s,100);  %detrending methods
            else
                continue      %try another running
            end
        else
            dataSelect = Y_CLEext(:,1);
            t_s = t_CLEext;
            detrendData = dataDetrend(dataSelect,t_s,100);  %detrending methods
        end
            j0 = j0 +1;
            variCLEext(i0,j0) = var(dataSelect);
            lagCLEext(i0,j0) = lagOne(dataSelect,dt);
            detr_variCLEext(i0,j0) = var(detrendData);
            detr_lagCLEext(i0,j0) = lagOne(detrendData,dt);
    end

end

%without extrinsic nosie
% param(5) = 0;
% for i0 = 1:size(threeState,1);
%     j0 = 0;
%     inv = max(allHL(i0,:));
%     param(2) = x0list(i0);
%     while(j0 <=N )
%         [t_CLE,Y_CLE] = CLEselfActi(inv,param,totTime);
%         if(i0>=min(inx_h) && i0<=max(inx_h)) % b we need selecet data mannually
%             plot(t_CLE,Y_CLE)
%             choice = input('the choice:\n','s');      %input the choice
%             if(strcmp(choice,'m')) %mannually selection
%                 [x,~] = ginput(2);
%                 [~,inx1] = sort(abs(t_CLE-x(1)));
%                 [~,inx2] = sort(abs(t_CLE-x(2)));
%                 dataSelect = Y_CLE(inx1(1):inx2(1),1);   %select data
%             elseif(strcmp(choice,'y'))  %use all the data
%                 dataSelect = Y_CLE(:,1);
%             else
%                 continue      %try another running
%             end
%         else
%             dataSelect = Y_CLE(:,1);
%         end
%             j0 = j0 +1;
%             variCLE(i0,j0) = var(dataSelect);
%             lagCLE(i0,j0) = lagOne(dataSelect,dt);
%             
%     end
% 
% end

%plot the results
figure(1)
d1 = detr_variCLEext(1:end-6,:);
d2 = variCLE(1:end-6,:);
b = zeros(size(d1,1),2,2);
b(:,:,1) = [std(d1,0,2),std(d1,0,2)];
b(:,:,2) = [std(d2,0,2),std(d2,0,2)];
boundedline(x0list(1:end-6),[mean(d1,2),mean(d2,2)],b,'alpha')
errorbar(x0list',mean(variCLEext,2),std(variCLEext,0,2))
hold on
errorbar(x0list',mean(detr_variCLEext,2),std(detr_variCLEext,0,2))
errorbar(x0list',mean(variCLE,2),std(variCLE,0,2))
hold off

figure(2)
hold on
b = zeros(size(lagCLEext,1),2,2);
b(:,:,1) = [std(lagCLEext,0,2),std(lagCLEext,0,2)];
b(:,:,2) = [std(detr_lagCLEext,0,2),std(detr_lagCLEext,0,2)];
b(:,:,3) = [std(lagCLE,0,2),std(lagCLE,0,2)];
boundedline(x0list,[mean(lagCLEext,2),mean(detr_lagCLEext,2),mean(lagCLE,2)],b,'alpha')
% errorbar(x0list',mean(lagCLEext,2),std(lagCLEext,0,2))
% hold on
% errorbar(x0list',mean(detr_lagCLEext,2),std(detr_lagCLEext,0,2))
% errorbar(x0list',mean(lagCLE,2),std(lagCLE,0,2))
% hold off
end

function autoCorr = lagOne(data,dt)
Ns= round(1/dt);
d1 = data(1:end-Ns);
d2 = data(Ns+1:end);
c = corrcoef(d1,d2);
autoCorr = c(1,2);
end
function FPE_inxtrinsic()
    %this function directly calculate the steady state Fokker-Planck
    %equation without extrinsic noise
    
    %first,symbolically 
    X_point = -0.2:0.01:1;
    a0 = 0.05;
    x0 = 0.5;
    N = 300;
% %     X_point = 0.8;
% %     y = potential(a0,x0,X_point);
% %     B = @(y) N*(a0 + (1-a0)*y.^2./(y.^2+x0^2*N^2))+y;
    y = steady(a0,x0);
    d1 = secoDiffPoten(a0,x0,y(3));   %second order at point a
    d2 = secoDiffPoten(a0,x0,y(2));     %second order at point b
    meanPassTime = pi/(sqrt(abs(d2*d1)))*exp(N*(potential(a0,x0,y(2))-potential(a0,x0,y(3))));
    

    for i0 = 1:length(X_point);
        U(i0) = potential(a0,x0,X_point(i0));
%     pusi = exp(2*quadl(@(x)myfun1(x,a0,x0,N),0,X_point(i0)));
%     ps = pusi/B(X_point)/quadl(@(x)myfun2(x,pusi,a0,x0),0,2);
%     ps(i0) = pusi/B(X_point(i0));
    end
    
    
    syms x y a0 x0
    A =@(y) a0 + (1-a0)*y^2/(y^2+x0^2)-y;
    
    pusi = @(x) exp(2*int(A(y)/B(y),y,0,x));
    normalization = int(pusi(x)/B(y),y,0,inf);
    steady_distri = pusi(x)/B(x)/normalization;
    ps = double(subs(steady_distri,{x,a0,x0},{0.8,0.05,0.5}));
    
    function f1 = myfun1(y,a0,x0,N)
        f1 = (N*(a0 + (1-a0)*y.^2./(y.^2+x0^2*N^2))-y)./(N*(a0 + (1-a0)*y.^2./(y.^2+x0^2*N^2))+y);
    end

    function f2 = myfun2(y,pusi,a0,x0)
            f2 = pusi/(a0 + (1-a0)*y.^2./(y.^2+x0^2)+y);
    end
    
    function y = potential(a0,x0,X)
             myfun = @(x,a0,x0) -a0 - (1-a0)*x.^2./(x.^2+x0^2)+x;
             y = integral(@(x)myfun(x,a0,x0),0,X);
%              syms x a0 x0
             
    end
    
    

            
end
function y = steadyDeter(a0,x0)
%              x = [0 1];  %find the root in this interval
         func = @(x,a0,x0) a0+(1-a0)*x^2/(x^2+x0^2)-x;
         f = @(x) func(x,a0,x0);
%          y0 = [0 1];
         for i0 = 1:3;
             y0 = [0,0.3,0.7];
             y(i0) = fzero(f,y0(i0));
         end
%          y = unique(z);
end

function d1 = secoDiffPoten(p1,p2,p3)
         % p1 = a0,p2 = x0,p3 = X;
         syms x a0 x0
         z= diff(-a0-(1-a0)*x^2/(x^2+x0^2)+x,x);
         d1 = double(subs(z,{x,a0,x0},{p3,p1,p2}));
end

function ui = potenInner(y)
        syms 
        int(pusi(x)/B(y),y,0,inf)
end
function [var_in,var_ex] = var_WKB(a0,x0,N,sigmae)
%this fucntion calculate the variance with  and without 
%extrinsic when the apprach is adiabatic
%the formula is form Assaf's 2015 paper
fixedValue = max(steadyDeter(a0,x0));
% x = sym('x');
% a = sym('a');
% b = sym('b');
[x,a,b]= deal([]);
syms x a b
z= diff(a+(1-a)*x^2/(x^2+b^2),x);
diff1 = double(subs(z,{x,a,b},{fixedValue,a0,x0}));

var_in = N*fixedValue/(1-diff1);
var_ex = var_in*(1+N*sigmae^2*fixedValue/(1-diff1));
end
function plotSwiTime()
foldername = 'swiTimeAlltau1000-0.1';
% foldername = 'swiAll_tau10_0.1';
% foldername = 'swiTimeAll_ne';
fileName = 'swiTime-tau1000-0.10_';
% fileName = 'swiTime_ne_';
t_swiAll = [];
for i = 1:20;
%     if(i ~=1 && i ~=4 &&i ~=5 && i ~=17)
        load(fullfile(foldername,[fileName,num2str(i-1),'.mat']));
        t_swiAll = [t_swiAll;t_swi];
%     end
end

tswi_mean = zeros(1,size(t_swiAll,2));
tsw_var = zeros(1,size(t_swiAll,2));
for j = 1:size(t_swiAll,2);
    tswi_mean(j) = mean(t_swiAll(t_swiAll(:,j)>0,j));
end
plot(x0,tswi_mean,'ro--','LineWidth',2,'MarkerSize',15)
xlabel('x0','FontSize',24,'FontWeight','Bold')
ylabel('\tau high\rightarrow low','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',2)
legend('simu-intri','WKB-intri','simu \sigma_{E}=0.1,\tau_{c} = 1000','simu \sigma_{E}=0.1, \tau_{c} = 10')

end
function WKBswiRate(a0,N)
              
        x0list = 0.43:0.005:0.52;
        for i = 1:length(x0list);
            x0 = x0list(i);
            y = steady(a0,x0);
            f = @(x) log(x./(a0 + (1-a0)*x.^2./(x.^2+x0^2)));
            actionDiff = integral(f,max(y),y(2));
            logTau(i) = N*actionDiff;
        end
        plot(x0list,13*exp(logTau))

end
function  bifurcation(a0,OMIGA)
x0list = 0.3:0.001:0.7;
steadyTheory = nan(length(x0list),3);
    for i = 1:length(x0list);
        y = steadyDeter(a0,x0list(i));
        if(max(y)-min(y)>0.2)
            steadyTheory(i,:) = sort(y,'descend');
        elseif(max(y)>0.4)
            steadyTheory(i,1) = max(y);
        elseif(max(y<0.4))
            steadyTheory(i,3) = max(y);
        end
    end
    plot(x0list',OMIGA*steadyTheory)
end 

function numericialSelfPromCLEExtrinsic(param,runTime,fileName,N)
load(fileName);  %load the file, which contains x0list, steady states value

I = 1:1:size(steady,1);
saddleIdex = I(~isnan(steady(:,2)));
saddle = steady(:,2)/300*param(3);
high = steady(:,1)/300*param(3);

meanVal  = zeros(saddleIdex(end-6),N);
variance = zeros(saddleIdex(end-6),N);
lagAuto = zeros(saddleIdex(end-6),N);

dt= 0.02;
MAXTRY = 2*N;
x0listSelect = x0list(1:saddleIdex(end-6));
for i = 1:saddleIdex(end-6);
    param(2) = x0list(i);
    n = 0;
    flag = 0;
    while(n<N && flag<MAXTRY)
        rng('shuffle')
        [time,Y] = CLESelfExtrinsic(runTime,high(i),dt,param);
%         [time,Y]=YGlpExtri(round(high(i)),300,0.05,param(2),runTime,param(5),param(4));
          if(i<saddleIdex(1) ||i>=saddleIdex(1) && Y(end,1)<saddle(i))
            n = n+1;
            meanVal(i,n) = mean(Y(:,1));        
            variance(i,n) = var(Y(:,1));
            C3 = corrcoef(Y(1:end-round(1/dt),1),Y(round(1/dt)+1:end,1));
            lagAuto(i,n) = C3(1,2);

          elseif(i>=saddleIdex(1))
            plot(time,Y)
            choice = input('the choice:\n','s');      %input the choice
            if(strcmp(choice,'y'))  %use all the data
                n = n+1;
                meanVal(i,n) = mean(Y(:,1));        
                variance(i,n) = var(Y(:,1));
                C3 = corrcoef(Y(1:end-round(1/dt),1),Y(round(1/dt)+1:end,1));
                lagAuto(i,n) = C3(1,2);
            else
                continue      %try another running
            end
            
            flag = flag + 1;
         end
    end     
end
saveFile = ['selfSwiTime_tau',num2str(param(4)),'_se',num2str(param(5)),'_',date,'.mat'];
save(saveFile,'t_swi','x0listSelect')


    figure(1)
    hold on
    errorbar(x0listSelect',mean(meanVal,2),std(meanVal,0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('mean expression level','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    hold off
    
    figure(2)
    hold on
    errorbar(x0listSelect',mean(variance,2),std(variance,0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('variance','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    
    
    figure(3)
    hold on
    errorbar(x0listSelect',mean(lagAuto,2),std(lagAuto,0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('lag one auto correlation','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
end
function [time,Y] = CLESelfExtrinsic(totTime,iniVal,dt,param)
a0 = param(1);
x0 = param(2);
OMIGA = param(3);
tauc = param(4);
ss = param(5);
sigmae = sqrt(2*ss^2/tauc);

N = round(totTime/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = [iniVal;0];         %revised on Dec 31,2015
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);

for i1 = 1:N;
    determin = OMIGA*(a0 + (1-a0)*Xtemp(1)^2/((x0*OMIGA)^2+Xtemp(1)^2))-(1+Xtemp(2))*Xtemp(1);
    stoch = [sqrt(OMIGA*(a0 + (1-a0)*Xtemp(1)^2/((x0*OMIGA)^2+Xtemp(1)^2))),-sqrt(abs((1+Xtemp(2))*Xtemp(1)))];

    Xtemp(1) = Xtemp(1) + determin*dt + stoch*randn(2,1)*sqrt(dt);
    Xtemp(2) = Xtemp(2)*mu+ Dext*randn;
    if (Xtemp(2) <-1)
        Xtemp(2) = -1;
    end
    
    X(i1,:) = Xtemp;
end

time = (0:dt:N*dt)';
Y = [[iniVal;0]';X];

end