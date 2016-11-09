%-----------------------------------------------------------------------
% PROGRAM NAME: toggCrossCorr.m
% DESCIPTION:  
%   this program uses LNA to calculate the cross correlation coefficient
%   of a toggle switch model subject to extrinsic noise
% writen by Shanshan Qin@Tanglab, PKU
% last revised on Jan.21,2015
%----------------------------------------------------------------------

function toggCrossCorr()
% close all
clear
clc

a1list = 4:1:18;
sigmae = 0.1;
tauc = 1;
param = [15;10;2;2;0.03;tauc;sigmae];
OMIGA = 20;        %system size
NUM = 500;
Nswi = 50;       %number of running when calculate the switching rate
sigma = [0.1,0.1];    %for LE
Yini = [0,10];

extriLE_CLE_noiseAmpli(param,OMIGA,sigma,sigmae,Yini);

end

function [t,y]=odeRun(param,totTime,Yini)

tspan = [0 totTime];
[t,y] = ode45(@toggleSwitch,tspan,Yini,[],param);
% plot(t,y,'LineWidth',2)
% legend('gene1','gene2')

end
function dydt = toggleSwitch(t,y,parameter)

% sythestic toggle switch model adopted form James J.Collins' 2000 Nature
% paper

a1 = parameter(1);          %Gene1 sythesis rate
a2 = parameter(2);          %Gene2 sythesis rate
beta = parameter(3);        %hill coefficient
gamma = parameter(4);       %hill coefficient
a0 = parameter(5);          %basal level synthesis rate

dydt = zeros(size(y));
u = y(1);
v = y(2);
dudt = a1*(a0+(1-a0)/(1+v^beta)) - u;
dvdt = a2*(a0 + (1-a0)/(1+u^gamma)) - v;
dydt = [dudt;dvdt];
end
function [Ac,Ace] = jacMatrx(param,vec,OMIGA)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
% modified the reaction with basal synthesis rate
syms x y E a1 a2 a0 n1 n2 b N

J3 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - y],[x y]);
Ac = double(subs(J3,{x,y,a1 a2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));
% J4 = jacobian([a1*N/(1+E+(y/N)^n1) - x;a2*N/(1+E+(x/N)^n2) - y;-b*E],[x y E]);
J4 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - (1+E)*x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - (1+E)*y;-b*E],[x y E]);
% J4 = jacobian([N*a1*(a0 + (1+E)*(1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1+E)*(1-a0)/(1+(x/N)^n2)) - y;-b*E],[x y E]);
% J4 = jacobian([N*a1*(a0 + (1-a0)/(1+E+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+E+(x/N)^n2)) - y;-b*E],[x y E]);
Ace = double(subs(J4,{x,y,E,a1 a2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,0,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

end
function extriLE_CLE_noiseAmpli(param,OMIGA,sigma,sigmae,Yini)
%this function comparing the correlation coefficient, variance, and lag one
%auto correlation for LE and CLE with extrinsic noise
% sigmae define diffent extrinsic noise amplitude

a1_s = 4;
a1_e = 16;
N_a1 = 40;
runTime = 1000;
param(6) = 1/param(6);

ss= [1e-3,1e-2,1e-1,0.25];
sigmae = 2*ss.^2*param(6);
corrcoef_ext_cle_theory = zeros(N_a1,1);
variance_ext_cle_theory = zeros(N_a1,length(sigmae));
lagOneCLEext_theory = zeros(N_a1,length(sigmae));
variance_ext_c_theory = zeros(N_a1,1);
corrcoef_ext_c_theory = zeros(N_a1,1);
lagOneC_theory = zeros(N_a1,1);
lagCrossCorr = cell(N_a1,length(sigmae)); %store cross correlation coefficient
lagCrossCorrExt = cell(N_a1,length(sigmae));

lagTimeList = -10:0.2:10;  %used to calculate lag crosscorrelation
for j1 = 1:N_a1;
    a1_list(j1) = a1_s + (a1_e - a1_s)*j1/N_a1;
    param(1)= a1_list(j1);
    [td,Yd]=odeRun(param,runTime,Yini);
    steadyODE(j1,:)=Yd(end,:);
    for j2 = 1:length(sigmae)

        [Ac,Acle] = jacMatrx(param,steadyODE(j1,:),OMIGA);
       
        Dc = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+OMIGA*steadyODE(j1,1),0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+OMIGA*steadyODE(j1,2)];
        Dcle = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+OMIGA*steadyODE(j1,1),0,0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+OMIGA*steadyODE(j1,2),0;0,0,sigmae(j2)];
        
        C3 = lyap(Ac,Dc);
        for k0 = 1:length(lagTimeList);
            if(lagTimeList(k0)>=0)
                crossTemp = C3*expm(Ac'*lagTimeList(k0));
            else
                crossTemp = C3*expm(Ac*abs(lagTimeList(k0)));
            end
            lagCrossCorr{j1,j2}(k0) = crossTemp(1,2)/(sqrt(C3(1,1)*C3(2,2)));
        end
        corrcoef_ext_c_theory(j1) = C3(1,2)/(sqrt(C3(1,1)*C3(2,2)));
%         variance_ext_c_theory(j1) = C3(2,2)/(OMIGA*steadyODE(j1,2));
        variance_ext_c_theory(j1) = C3(2,2)/(OMIGA*steadyODE(j1,2));
        G3 = C3*expm(Ac');
        lagOneC_theory(j1) = G3(2,2)/C3(2,2);
        
        C4 = lyap(Acle,Dcle);
        for k0 = 1:length(lagTimeList);
            if(lagTimeList(k0)>=0)
                crossTemp = C4*expm(Acle'*lagTimeList(k0));
            else
                crossTemp = expm(Acle*abs(lagTimeList(k0)))*C4;
            end
            lagCrossCorrExt{j1,j2}(k0) = crossTemp(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        end        
%         corrcoef_ext_cle_theory(j1,j2) = C4(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        variance_ext_cle_theory(j1,j2) = C4(2,2)/(OMIGA*steadyODE(j1,2));
        variance_ext_cle_theory(j1,j2) = C4(2,2);
        G4 = C4*expm(Acle');
        lagOneCLEext_theory(j1,j2) = G4(2,2)/C4(2,2);

    end
end

%plot the theoretic results
colorSet = [[0;0.3;0.6;0.9],zeros(4,2)];
figure(1)
% plot(a1_list,corrcoef_theory,'r-','Linewidth',3)
hold on
plot(a1_list,corrcoef_ext_c_theory,'k-','Linewidth',3)
for j3 = 1:length(sigmae);
%     plot(a1_list,corrcoef_ext_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:))
    plot(a1_list,corrcoef_ext_cle_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:),'LineStyle','--')
end
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
hold off

figure(2)
% plot(a1_list,variance_theory,'r-','Linewidth',3)
hold on
% plot(a1_list,OMIGA^2*variance_ext_c_theory,'k-','Linewidth',3)
plot(a1_list,variance_ext_c_theory,'k-','Linewidth',3)
for j4 = 1:length(sigmae);
%     plot(a1_list,variance_ext_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:))
%     plot(a1_list,OMIGA^2*variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--')
    plot(a1_list,variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--')
    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('variacne of gene B','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

figure(3)
% plot(a1_list,lagOneLe_theory,'r-','LineWidth',3)
hold on
plot(a1_list,lagOneC_theory,'k-','Linewidth',3)
for j5 = 1:length(sigmae);
%     plot(a1_list,lagOneLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:))
    plot(a1_list,lagOneCLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:),'LineStyle','--')    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

%plot the cross correlation 
figure(4)
sInx = [1,20,30,40];
for i = 1:4;
subplot(2,2,i)
hold on
for j = 1:4;
    plot(lagTimeList,lagCrossCorrExt{sInx(i),j})
end
xlabel('time','FontSize',16,'FontWeight','Bold')
ylabel('cross correlation','FontSize',16,'FontWeight','Bold')
hold off
end

end
function simuDistri(param,NUM,a1list,OMIGA)

%this fucntion

allIni = lhsdesign(NUM,2).*[ones(NUM,1)*OMIGA*10,ones(NUM,1)*OMIGA*10];
t_sim = 1e4;
N_s = 1e3;    % the transient time which the system has not relax to steady state
dt = 0.1;
Data1  = zeros(NUM*N_s,length(a1list));
Data2  = zeros(NUM*N_s,length(a1list));
for j0 = 1:length(a1list);
    param(1) = a1list(j0);
    for i0 = 2:NUM;
        [time,Y]=CLEtoggleExtrinsic(t_sim,allIni(i0,:)',dt,param,OMIGA);
        Data1((i0-1)*N_s+1:i0*N_s,j0) = Y(end-N_s+1:end,1);
        Data2((i0-1)*N_s+1:i0*N_s,j0) = Y(end-N_s+1:end,2);
    end

end
    save('toggleDistri_tau1e3.mat','Data1','Data2','a1list')
end
function plotDistri(a1list)

%first find all the minimum and maximum of distribution
%then calculate the variance at each metastate

%      load('/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData/selfDistri-1220.mat')
%     load('selfDistri-1e-2-1227.mat'); 
    load('toggleDistri_tau1e3.mat');
    steady = nan(length(a1list),3);
    for i0 = 1:length(a1list);
        [f,xi]=ksdensity(Data2(:,i0),'npoints',200,'support',[-10 400]);
        [maxTab, minTab] = peakdet(f,1e-4,xi);
        semilogy(xi,f)
        hold on
        plot(maxTab(:,1),maxTab(:,2),'ro')

        if (size(maxTab,1)==1&& maxTab(1,1) > 30)
            steady(i0,1) = maxTab(1,1);
        elseif(size(maxTab,1)==1&& maxTab(1,1) < 30)
            steady(i0,3) = maxTab(1,1);
        elseif(size(maxTab,1)==2)
            steady(i0,1) = max(maxTab(:,1));
            steady(i0,3) = min(maxTab(:,1));
        end
        if(~isempty(minTab))
            plot(minTab(:,1),minTab(:,2),'g*')
            steady(i0,2) = minTab(:,1);
        end
        
    end
%     save('stochSteady_tau_10_0.1.mat','x0list','steady')
%     figure(2)
%     plot(a1list',steady,'ro')
    
    inx_h = ~isnan(steady(:,2));
    high = steady( ~isnan(steady(:,1)),1);  %all the high state
    saddle = steady(:,2);
    saddle(isnan(saddle))=0;
    var_high = nan(length(a1list),1);
%     fileB = 'toggSwi-fixedB.csv';
%     dataB = csvread(fileB);
%     saddle = dataB(:,3)*20;
    for j0 = 1:length(high);
        var_high(j0) = var(Data2(Data2(:,j0)>saddle(j0),j0));
    end
%     figure(3)
    plot(a1list',var_high,'r--o')
    xlabel('x0','FontSize',24,'FontWeight','Bold')
    ylabel('variance at high state','FontSize',24,'FontWeight','Bold')
    set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',2)
    
end
function fixedPoint(param,a1list)
        
        guessIni = [0,10;[lhsdesign(20,1)*15,lhsdesign(20,1)*8]];
        low = nan(length(a1list),2);
        saddle = nan(length(a1list),2);
        high = nan(length(a1list),2);
        for j = 1:length(a1list)
            param(1) = a1list(j);
            for i = 1:size(guessIni,1);
                fixed(i,:) = round(1e3*fsolve(@(y)eqns(y,param), guessIni(i,:)))/1e3; 
            end
            uniValue = unique(fixed,'rows');
            if(size(uniValue,1)==3);
              low(j,:) = uniValue(1,:);
              saddle(j,:) = uniValue(2,:);
              high(j,:) = uniValue(3,:);
            elseif(size(uniValue,1)==1)
                if(uniValue(1)>uniValue(2))
                    high(j,:) = uniValue;
                else
                    low(j,:) = uniValue;
                end
            elseif(size(uniValue,1)==2)
                low(j,:) = uniValue(1,:);
                high(j,:) = uniValue(2,:);
            end
        end
end
function myfun = eqns(y,parameter)
        a1 = parameter(1);        
        a2 = parameter(2);        
        beta = parameter(3);       
        gamma = parameter(4);       
        a0 = parameter(5);          

        myfun = zeros(size(y));
        u = y(1);
        v = y(2);
        dudt = a1*(a0+(1-a0)/(1+v^beta)) - u;
        dvdt = a2*(a0 + (1-a0)/(1+u^gamma)) - v;
        myfun = [dudt;dvdt];
end
