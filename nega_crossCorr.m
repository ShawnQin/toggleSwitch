%********************************************************************
% PROGRAM NAME: nega_crossCorr.m
% DESCRIPTION: 
%    this program calculate the variance, correlation and lag one aouto
%    correlation using linear noise approximation methods
% last revised on Jan 4,2016
%**********************************************************************
function nega_crossCorr()

% close all
clear
clc

% k1list = 0.1:0.01:0.55;
sigmae = 0.05;
tauc = 1;
param = [0.3;0.5;3;3;0.05;tauc;sigmae];
OMIGA = 500;        %system size

% extriLE_CLE_timescale(param,OMIGA,sigma,sigmae,Yini);
extriLE_CLE_noiseAmpli(param,OMIGA);
end
function [Ac,Ace] = jacMatrx(param,vec,OMIGA)

syms x y E k1 k2 a0 n1 n2 b N

J3 = jacobian([N*(a0 + (1-a0)*k1^n1/(k1^n1+(y/N)^n1)) - x;N*(a0 + (1-a0)*(x/N)^n2/(k2^n2+(x/N)^n2)) - y],[x y]);
Ac = double(subs(J3,{x,y,k1 k2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

J4 = jacobian([N*(a0 + (1-a0)*k1^n1/(k1^n1+(y/N)^n1)) - (1+E)*x;N*(a0 + (1-a0)*(x/N)^n2/(k2^n2+(x/N)^n2)) - (1+E)*y;-b*E],[x y E]);
% J4 = jacobian([N*(a0 + (1+E)*(1-a0)*k1^n1/(k1^n1+(y/N)^n1))-x;N*(a0 + (1+E)*(1-a0)*(x/N)^n2/(k2^n2+(x/N)^n2))-y;-b*E],[x y E]);
Ace = double(subs(J4,{x,y,E,k1 k2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,0,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

end
function extriLE_CLE_noiseAmpli(param,OMIGA)



iniVal = [0.5;0.5];

k1_s = 0.1;
k1_e = 2;

a0 = param(1);
k2 = param(3);
n1 = param(4);
n2 = param(5);
param(6) = 1/param(6);


ss= [1e-2,5e-2,1e-1,0.25];  %amplitute of extrinsic noise
sigmae = 2*ss.^2*param(6);
k1_list = k1_s:0.02:k1_e;
N_k1 = length(k1_list);


tspan = [0 1000];
corrcoef_ext_cle_theory = zeros(N_k1,1);
variance_ext_cle_theory = zeros(N_k1,length(sigmae));

lagOneCLEext_theory = zeros(N_k1,length(sigmae));
variance_CLE_theory = zeros(N_k1,1);
corrcoef_CLE_theory = zeros(N_k1,1);
lagOneCLE_theory = zeros(N_k1,1);
lagCrossCorr = cell(N_k1,length(sigmae)); %store cross correlation coefficient
lagCrossCorrExt = cell(N_k1,length(sigmae));

lagTimeList = -10:0.2:10;  %used to calculate lag crosscorrelation
for j1 = 1:N_k1;
    param(1)= k1_list(j1);
    [t,y] = ode15s(@negativeODE,tspan,iniVal,[],param);
    steadyValue(j1,:) = y(end,:);
    for j2 = 1:length(sigmae)

        [Ac,Acle] = jacMatrx(param,steadyValue(j1,:),OMIGA);

        Dc = [OMIGA*(a0 + (1-a0)*k1_list(j1)^n1/(k1_list(j1)^n1 + steadyValue(j1,2)^n1))+OMIGA*steadyValue(j1,1),0;...
                0,OMIGA*(a0 + (1-a0)*steadyValue(j1,1)^n2/(k2^n2+steadyValue(j1,1)^n2))+OMIGA*steadyValue(j1,1)];
            
        Dcle = [OMIGA*(a0 + (1-a0)*k1_list(j1)^n1/(k1_list(j1)^n1+steadyValue(j1,2)^n1))+OMIGA*steadyValue(j1,1),0,0;...
                0,OMIGA*(a0 + (1-a0)*steadyValue(j1,1)^n2/(k2^n2+steadyValue(j1,1)^n2))+OMIGA*steadyValue(j1,1),0;0,0,sigmae(j2)]; 
                   
        C3 = lyap(Ac,Dc);            
        corrcoef_CLE_theory(j1) = C3(1,2)/(sqrt(C3(1,1)*C3(2,2)));
        variance_CLE_theory(j1) = C3(2,2)/steadyValue(j1,2)/OMIGA;   % CV shoud be used to characterize the variance ,revised on Jan 07,2016
        
        for k0 = 1:length(lagTimeList);
            if(lagTimeList(k0)>=0)
                crossTemp = C3*expm(Ac'*lagTimeList(k0));
            else
                crossTemp = expm(Ac*abs(lagTimeList(k0)))*C3;
            end
            lagCrossCorr{j1,j2}(k0) = crossTemp(1,2)/(sqrt(C3(1,1)*C3(2,2)));
        end        
        
        G3 = C3*expm(Ac');
        lagOneCLE_theory(j1) = G3(2,2)/C3(2,2);
        
        C4 = lyap(Acle,Dcle);
        corrcoef_ext_cle_theory(j1,j2) = C4(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        variance_ext_cle_theory(j1,j2) = C4(2,2)/steadyValue(j1,2)/OMIGA;
        for k0 = 1:length(lagTimeList);
            if(lagTimeList(k0)>=0)
                crossTemp = C4*expm(Acle'*lagTimeList(k0));
            else
                crossTemp = expm(Acle*abs(lagTimeList(k0)))*C4;
            end
            lagCrossCorrExt{j1,j2}(k0) = crossTemp(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        end             
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

figure(4)
sInx = [1,30,70,90];
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

function dydt = negativeODE(t,y,param)
    k1 = param(1);
    k2 = param(2);
    n1 = param(3);
    n2 = param(4);
    a0 = param(5);
    
    dydt = [a0 + (1-a0)*k1^n1/(k1^n1 + y(2)^n2) - y(1); a0 + (1-a0)*y(1)^n2/(k2^n2 + y(1)^n2) - y(2)];
end