%********************************************************************
% PROGRAM NAME: mutualActive_LNA.m
% DESCRIPTION: 
%    this program calculate the variance, correlation and lag one aouto
%    correlation using linear noise approximation methods
% last revised on Jan 4,2016
%**********************************************************************
function mutualActive_LNA()

% close all
clear
clc

k1list = 0.1:0.01:0.55;
sigmae = 0.05;
tauc = 10;
param = [0.3;0.5;3;3;0.05;tauc;sigmae];
OMIGA = 500;        %system size

% extriLE_CLE_timescale(param,OMIGA,sigma,sigmae,Yini);
extriLE_CLE_noiseAmpli(param,OMIGA);
end
function [Ac,Ace] = jacMatrx(param,vec,OMIGA)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
% modified the reaction with basal synthesis rate
syms x y E k1 k2 a0 n1 n2 b N

J3 = jacobian([N*(a0 + (1-a0)*(y/N)^n1/(k1^n1+(y/N)^n1)) - x;N*(a0 + (1-a0)*(x/N)^n2/(k2^n2+(x/N)^n2)) - y],[x y]);
Ac = double(subs(J3,{x,y,k1 k2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

% J4 = jacobian([N*(a0 + (1-a0)*(y/N)^n1/(k1^n1+(y/N)^n1)) - (1+E)*x;N*(a0 + (1-a0)*(x/N)^n2/(k2^n2+(x/N)^n2)) - (1+E)*y;-b*E],[x y E]);
J4 = jacobian([N*(a0 + (1+0)*(1-a0)*(y/N)^n1/(k1^n1+(y/N)^n1))-x;N*(a0 + (1+E)*(1-a0)*(x/N)^n2/(k2^n2+(x/N)^n2))-y;-b*E],[x y E]);
Ace = double(subs(J4,{x,y,E,k1 k2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,0,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

end
function extriLE_CLE_timescale(param,OMIGA,sigma,sigmae,Yini)
%this function comparing the correlation coefficient, variance, and lag one
%auto correlation for LE and CLE with extrinsic noise
% beta define different time scale of extrinsic noise

a1_s = 4;
a1_e = 16;
N_a1 = 40;
runTime = 1000;
% ss = 2*sigmae^2/param(6);  %the extrinsic noise O-U process

beta = [10,1/10,1/100,1/1000];
% beta = 1/param(6);
corrcoef_theory = zeros(N_a1,1);
corrcoef_ext_theory = zeros(N_a1,length(beta));
variance_theory = zeros(N_a1,1);
variance_ext_theory = zeros(N_a1,length(beta));
corrcoef_ext_cle_theory = zeros(N_a1,1);
variance_ext_cle_theory = zeros(N_a1,length(beta));
lagOneLe_theory = zeros(N_a1,1);
lagOneLEext_theory = zeros(N_a1,length(beta));
lagOneCLEext_theory = zeros(N_a1,length(beta));
variance_ext_c_theory = zeros(N_a1,1);
corrcoef_ext_c_theory = zeros(N_a1,1);
lagOneC_theory = zeros(N_a1,1);
for j1 = 1:N_a1;
    a1_list(j1) = a1_s + (a1_e - a1_s)*j1/N_a1;
    param(1)= a1_list(j1);
    [td,Yd]=odeRun(param,runTime,Yini);
    steadyODE(j1,:)=Yd(end,:);
    a0 = param(5);
    for j2 = 1:length(beta)
        param(6) = beta(j2);
        ss = 2*sigmae^2*beta(j2);
        [A,Ae,Ac,Acle] = jacMatrx(param,steadyODE(j1,:),OMIGA);
        D = [sigma(1)^2,0;0,sigma(2)^2];
        De = [sigma(1)^2,0,0;0,sigma(2)^2,0;0,0,ss];
        Dc = [a1_list(j1)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,2)^param(3)))+steadyODE(j1,1)*OMIGA,0;...
                0,param(2)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,1)^param(4)))+steadyODE(j1,2)*OMIGA];
        Dcle = [a1_list(j1)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,2)^param(3)))+steadyODE(j1,1)*OMIGA,0,0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+steadyODE(j1,2)*OMIGA,0;0,0,ss];
        C1 = lyap(A,D);
        corrcoef_theory(j1) = C1(1,2)/(sqrt(C1(1,1)*C1(2,2)));
        variance_theory(j1) = C1(2,2);
        G1 = C1*expm(A');
        lagOneLe_theory(j1) = G1(2,2)/C1(2,2);
        
        C2 = lyap(Ae,De);
        corrcoef_ext_theory(j1,j2) = C2(1,2)/(sqrt(C2(1,1)*C2(2,2)));
        variance_ext_theory(j1,j2) = C2(2,2);
        G2 = C2*expm(Ae');
        lagOneLEext_theory(j1,j2) = G2(2,2)/C2(2,2);
        
        C3 = lyap(Ac,Dc);            
        corrcoef_ext_c_theory(j1) = C3(1,2)/(sqrt(C3(1,1)*C3(2,2)));
        variance_ext_c_theory(j1) = C3(2,2)/(OMIGA^2);
        G3 = C3*expm(Ac');
        lagOneC_theory(j1) = G3(2,2)/C3(2,2);
        
        C4 = lyap(Acle,Dcle);
        corrcoef_ext_cle_theory(j1,j2) = C4(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        variance_ext_cle_theory(j1,j2) = C4(2,2)/(OMIGA^2);
        G4 = C4*expm(Acle');
        lagOneCLEext_theory(j1,j2) = G4(2,2)/C4(2,2);

    end
end

%plot the theoretic resujicolorSet = [ones(4,1)*60/256,[0;0.3;0.6;0.9],ones(4,1)*189/256];
colorSet = [[0;0.3;0.6;0.9],zeros(4,2)];
figure(1)
% plot(a1_list,corrcoef_theory,'r-','Linewidth',3)
hold on
plot(a1_list,corrcoef_ext_c_theory,'k-','Linewidth',3)
for j3 = 1:length(beta);
%     plot(a1_list,corrcoef_ext_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:))
    plot(a1_list,corrcoef_ext_cle_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:),'LineStyle','--')
end
% legend('LE','CLE','LE-extr-1','CLE-extr-1','LE-extr-10','CLE-extr-10','LE-extr-50','CLE-extr-50','LE-extr-500','CLE-extr-500')
legend('CLE-intrinsic','CLE-extr-\tau_{E}=10^{-1}','CLE-extr-\tau_{E}=10','CLE-extr-\tau_{E}=10^{2}','CLE-extr-\tau_{E}=10^{3}')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
% legend('simu-intri','simu-\sigma_{E}=0.2,\tau = 10^{3}','simu-\sigma_{E}=0.01,\tau = 10^{3}','simu-\sigma_{E}=0.1,\tau = 1',...
%     'LNA-\sigma_{E}=0.2,\tau = 10^{3}','LNA-\sigma_{E}=0.1,\tau = 1','LNA-\sigma_{E}=0.01,\tau = 10^{3}','LNA-intrin')
% hold off

figure(2)
% plot(a1_list,variance_theory,'r-','Linewidth',3)
hold on
plot(a1_list,variance_ext_c_theory*OMIGA^2,'k-','Linewidth',3)
for j4 = 1:length(beta);
%     plot(a1_list,variance_ext_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:))
    plot(a1_list,OMIGA^2*variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--')
    
end
% hold off
% legend('LE','CLE','LE-extr-1','CLE-extr-1','LE-extr-10','CLE-extr-10','LE-extr-50','CLE-extr-50','LE-extr-500','CLE-extr-500')
legend('CLE-intrinsic','CLE-extr-\tau_{E}=10^{-1}','CLE-extr-\tau_{E}=10','CLE-extr-\tau_{E}=10^{2}','CLE-extr-\tau_{E}=10^{3}')

xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('variacne of gene B','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

figure(3)
% plot(a1_list,lagOneLe_theory,'r-','LineWidth',3)
hold on
plot(a1_list,lagOneC_theory,'k-','Linewidth',3)
for j5 = 1:length(beta);
%     plot(a1_list,lagOneLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:))
    plot(a1_list,lagOneCLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:),'LineStyle','--')    
end
% hold off
% legend('LE','CLE','LE-extr-1','CLE-extr-1','LE-extr-10','CLE-extr-10','LE-extr-50','CLE-extr-50','LE-extr-500','CLE-extr-500')
legend('CLE-intrinsic','CLE-extr-\tau_{E}=10^{-1}','CLE-extr-\tau_{E}=10','CLE-extr-\tau_{E}=10^{2}','CLE-extr-\tau_{E}=10^{3}')

xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
end
function extriLE_CLE_noiseAmpli(param,OMIGA)

fileA = 'mutualActiveBifurA.csv';
fileB = 'mutualActiveBifurB.csv';
dataA = csvread(fileA);
dataB = csvread(fileB);

k1_s = 0.1;
k1_e = 0.7;

a0 = param(1);
k2 = param(3);
n1 = param(4);
n2 = param(5);
param(6) = 1/param(6);


ss= [1e-3,1e-2,1e-1,0.25];  %amplitute of extrinsic noise
sigmae = 2*ss.^2*param(6);

k1_list = k1_s:0.01:k1_e;
N_k1 = length(k1_list);
for i = 1:N_k1;
    indexofk1(i) = find(dataA(:,1)==round(k1_list(i)*100)/100);
    steadyValue(i,:) = [dataA(indexofk1(i),4),dataB(indexofk1(i),4)];
end


corrcoef_ext_cle_theory = zeros(N_k1,1);
variance_ext_cle_theory = zeros(N_k1,length(sigmae));

lagOneCLEext_theory = zeros(N_k1,length(sigmae));
variance_CLE_theory = zeros(N_k1,1);
corrcoef_CLE_theory = zeros(N_k1,1);
lagOneCLE_theory = zeros(N_k1,1);

for j1 = 1:N_k1;
    param(1)= k1_list(j1);

    for j2 = 1:length(sigmae)

        [Ac,Acle] = jacMatrx(param,steadyValue(j1,:),OMIGA);

        Dc = [OMIGA*(a0 + (1-a0)*steadyValue(j1,2)^n1/(k1_list(j1)^n1 + steadyValue(j1,2)^n1))+OMIGA*steadyValue(j1,1),0;...
                0,OMIGA*(a0 + (1-a0)*steadyValue(j1,1)^n2/(k2^n2+steadyValue(j1,1)^n2))+OMIGA*steadyValue(j1,1)];
            
        Dcle = [OMIGA*(a0 + (1-a0)*steadyValue(j1,2)^n1/(k1_list(j1)^n1+steadyValue(j1,2)^n1))+OMIGA*steadyValue(j1,1),0,0;...
                0,OMIGA*(a0 + (1-a0)*steadyValue(j1,1)^n2/(k2^n2+steadyValue(j1,1)^n2))+OMIGA*steadyValue(j1,1),0;0,0,sigmae(j2)]; 
                   
        C3 = lyap(Ac,Dc);            
        corrcoef_CLE_theory(j1) = C3(1,2)/(sqrt(C3(1,1)*C3(2,2)));
        variance_CLE_theory(j1) = C3(2,2)/(OMIGA*steadyValue(j2,2));
        G3 = C3*expm(Ac');
        lagOneCLE_theory(j1) = G3(2,2)/C3(2,2);
        
        C4 = lyap(Acle,Dcle);
        corrcoef_ext_cle_theory(j1,j2) = C4(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        variance_ext_cle_theory(j1,j2) = C4(2,2)/(OMIGA*steadyValue(j2,2));
        G4 = C4*expm(Acle');
        lagOneCLEext_theory(j1,j2) = G4(2,2)/C4(2,2);

    end
end
%plot the theoretic results
 colorSet = [[0;0.3;0.6;0.9],zeros(4,2)];
figure(1)

hold on
plot(k1_list,corrcoef_CLE_theory,'k-','Linewidth',3)
for j3 = 1:length(sigmae);
    plot(k1_list,corrcoef_ext_cle_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:),'LineStyle','--')
end
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
hold off

figure(2)
hold on
plot(k1_list,variance_CLE_theory,'k-','Linewidth',3)
for j4 = 1:length(sigmae);
    plot(k1_list,variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--')
    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('variacne of gene B','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

figure(3)

hold on
plot(k1_list,lagOneCLE_theory,'k-','Linewidth',3)
for j5 = 1:length(sigmae);
    plot(k1_list,lagOneCLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:),'LineStyle','--')    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
end