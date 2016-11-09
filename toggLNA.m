%-----------------------------------------------------------------------
% PROGRAM NAME: toggLNA.m
% DESCIPTION:  
%   this program use Linear Noise Approximation to calculate the
%   early-worning signals of a toggle switch model
% writen by Shanshan Qin@Tanglab, PKU
% last revised on Jan.24,2016
%----------------------------------------------------------------------

%----------------------------------------------------------

function toggLNA()
close all
clear
clc

a1list = 4:1:18;
sigmae = 0.1;
tauc1 = 100;
tauc2 = 100;
param = [15;10;2;2;0.03;tauc1;tauc2;sigmae];
OMIGA = 50;        %system size
NUM = 500;
Nswi = 50;       %number of running when calculate the switching rate
sigma = [0.1,0.1];    %for LE
Yini = [0,10];

% noiseAmpli_same(param,OMIGA,sigma,sigmae,Yini);
noiseAmpli_diff(param,OMIGA,sigma,Yini);

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

function [A,Ae,Ac,Ace] = jacMatrx(param,vec,OMIGA)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
% modified the reaction with basal synthesis rate
syms x y E a1 a2 a0 n1 n2 b N
J1 = jacobian([a1*(a0 + (1-a0)/(1+y^n1)) - x;a2*(a0 + (1-a0)/(1+x^n2)) - y],[x,y]);
A = double(subs(J1,{x,y,a1 a2 n1 n2 a0},{vec(1),vec(2),param(1),param(2),param(3),param(4),param(5)}));
% J2 = jacobian([a1*(a0 + (1-a0)/(1+y^n1)) - x + E;a2*(a0 + (1-a0)/(1+x^n2)) - y + E;-b*E],[x y E]);
% J2 = jacobian([a1*(a0 + (1-a0)/(1+y^n1)) - (1+E)*x;a2*(a0 + (1-a0)/(1+x^n2)) - (1+E)*y;-b*E],[x y E]);
% % % J2 = jacobian([a1/(1+E+y^n1) - x;a2/(1+E+x^n2) - y;-b*E],[x y E]);
% J2 = jacobian([a1*(1+E)*(a0 + (1-a0)/(1+y^n1)) - x;a2*(1+E)*(a0 + (1-a0)/(1+x^n2)) - y;-b*E],[x y E]);
J2 = jacobian([a1*(a0 + (1-a0)/(1+y^n1)) - x*(1+E);a2*(a0 + (1-a0)/(1+x^n2)) - (1+E)*y;-b*E],[x y E]);
Ae = double(subs(J2,{x,y,E,a1 a2 n1 n2  a0 b},{vec(1),vec(2),0,param(1),param(2),param(3),param(4),param(5),param(6)}));

J3 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - y],[x y]);
Ac = double(subs(J3,{x,y,a1 a2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));
% J4 = jacobian([a1*N/(1+E+(y/N)^n1) - x;a2*N/(1+E+(x/N)^n2) - y;-b*E],[x y E]);
J4 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - (1+E)*x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - (1+E)*y;-b*E],[x y E]);
% J4 = jacobian([N*a1*(a0 + (1+E)*(1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1+E)*(1-a0)/(1+(x/N)^n2)) - y;-b*E],[x y E]);
% J4 = jacobian([N*a1*(a0 + (1-a0)/(1+E+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+E+(x/N)^n2)) - y;-b*E],[x y E]);
Ace = double(subs(J4,{x,y,E,a1 a2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,0,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

end
function [Ac,Ace] = jacMatrx_diff(param,vec,OMIGA)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
% modified the reaction with basal synthesis rate
syms x y E1 E2 a1 a2 a0 n1 n2 b1 b2 N

J3 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - y],[x y]);
Ac = double(subs(J3,{x,y,a1 a2 n1 n2 a0 N},{vec(1)*OMIGA,vec(2)*OMIGA,param(1),param(2),param(3),param(4),param(5),OMIGA}));

% J4 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - (1+E1)*x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - (1+E2)*y;-b1*E1;-b2*E2],[x y E1 E2]);
J4 = jacobian([N*a1*(a0 + (1-a0)*(1+E1)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - y;-b1*E1;-b2*E2],[x y E1 E2]);
% J4 = jacobian([N*a1*(a0 + (1-a0)*(1+E1)/(1+(y/N)^n1)) - (1+E2)*x;N*a2*(a0 + (1-a0)*(1+E1)/(1+(x/N)^n2)) - (1+E2)*y;-b1*E1;-b2*E2],[x y E1 E2]);
% J4 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)*(1+E2)/(1+(x/N)^n2)) - y;-b1*E1;-b2*E2],[x y E1 E2]);
% J4 = jacobian([N*a1*(a0 + (1-a0)*exp(E1)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)*exp(E2)/(1+(x/N)^n2)) - y;-b1*E1;-b2*E2],[x y E1 E2]);
Ace = double(subs(J4,{x,y,E1,E2,a1 a2 n1 n2 a0 b1 b2 N},{vec(1)*OMIGA,vec(2)*OMIGA,0,0,param(1),param(2),param(3),param(4),param(5),param(6),param(7),OMIGA}));

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
function noiseAmpli_same(param,OMIGA,sigma,sigmae,Yini)
%this function comparing the correlation coefficient, variance, and lag one
%auto correlation for LE and CLE with extrinsic noise
% sigmae define diffent extrinsic noise amplitude
% the extrinsic noise for the two varibles are the same

a1_s = 4;
a1_e = 16;
N_a1 = 40;
runTime = 1000;
param(6) = 1/param(6);

 ss= [1e-3,1e-2,1e-1,0.25];
sigmae = 2*ss.^2*param(6);
corrcoef_theory = zeros(N_a1,1);
corrcoef_ext_theory = zeros(N_a1,length(sigmae));
variance_theory = zeros(N_a1,1);
variance_ext_theory = zeros(N_a1,length(sigmae));
corrcoef_ext_cle_theory = zeros(N_a1,1);
variance_ext_cle_theory = zeros(N_a1,length(sigmae));
lagOneLe_theory = zeros(N_a1,1);
lagOneLEext_theory = zeros(N_a1,length(sigmae));
lagOneCLEext_theory = zeros(N_a1,length(sigmae));
variance_ext_c_theory = zeros(N_a1,1);
corrcoef_ext_c_theory = zeros(N_a1,1);
lagOneC_theory = zeros(N_a1,1);

for j1 = 1:N_a1;
    a1_list(j1) = a1_s + (a1_e - a1_s)*j1/N_a1;
    param(1)= a1_list(j1);
    [td,Yd]=odeRun(param,runTime,Yini);
    steadyODE(j1,:)=Yd(end,:);
    for j2 = 1:length(sigmae)

        [A,Ae,Ac,Acle] = jacMatrx(param,steadyODE(j1,:),OMIGA);
        D = [sigma(1)^2,0;0,sigma(2)^2];
        De = [sigma(1)^2,0,0;0,sigma(2)^2,0;0,0,sigmae(j2)];
        Dc = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+OMIGA*steadyODE(j1,1),0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+OMIGA*steadyODE(j1,2)];
        Dcle = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+OMIGA*steadyODE(j1,1),0,0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+OMIGA*steadyODE(j1,2),0;0,0,sigmae(j2)];
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
%         variance_ext_c_theory(j1) = C3(2,2)/(OMIGA*steadyODE(j1,2));
        variance_ext_c_theory(j1) = C3(2,2)/(OMIGA*steadyODE(j1,2));
        G3 = C3*expm(Ac');
        lagOneC_theory(j1) = G3(2,2)/C3(2,2);
        
        C4 = lyap(Acle,Dcle);
        corrcoef_ext_cle_theory(j1,j2) = C4(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        variance_ext_cle_theory(j1,j2) = C4(2,2)/(OMIGA*steadyODE(j1,2));
%         variance_ext_cle_theory(j1,j2) = C4(2,2);
        G4 = C4*expm(Acle');
        lagOneCLEext_theory(j1,j2) = G4(2,2)/C4(2,2);

    end
end
%plot the theoretic results
% colorSet = [ones(4,1)*60/256,[0;0.3;0.7;1],ones(4,1)*189/256];


% colorSet = [ones(3,1)*60/256,[0;0.5;1],ones(3,1)*189/256];
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
end
function noiseAmpli_diff(param,OMIGA,sigma,Yini)
%this function comparing the correlation coefficient, variance, and lag one
%auto correlation for LE and CLE with extrinsic noise
% sigmae define diffent extrinsic noise amplitude
% the extrinsic noise for the two varibles are different or independent

a1_s = 4;
a1_e = 16;
N_a1 = 40;
runTime = 1000;
param(6) = 1/param(6);
param(7) = 1/param(7);

ss= [1e-3,1e-2,1e-1,0.25];
sigmae = 2*ss.^2*param(6);


corrcoef_ext_cle_theory = zeros(N_a1,1);
variance_ext_cle_theory = zeros(N_a1,length(sigmae));
lagOneCLEext_theory = zeros(N_a1,length(sigmae));
variance_ext_c_theory = zeros(N_a1,1);
corrcoef_ext_c_theory = zeros(N_a1,1);
lagOneC_theory = zeros(N_a1,1);

for j1 = 1:N_a1;
    a1_list(j1) = a1_s + (a1_e - a1_s)*j1/N_a1;
    param(1)= a1_list(j1);
    [td,Yd]=odeRun(param,runTime,Yini);
    steadyODE(j1,:)=Yd(end,:);
    for j2 = 1:length(sigmae)

        [Ac,Acle] = jacMatrx_diff(param,steadyODE(j1,:),OMIGA);
        Dc = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+OMIGA*steadyODE(j1,1),0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+OMIGA*steadyODE(j1,2)];
        Dcle = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+OMIGA*steadyODE(j1,1),0,0,0,;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+OMIGA*steadyODE(j1,2),0,0;...
                0,0,sigmae(j2),0;0,0,0,sigmae(j2)];
        
        C3 = lyap(Ac,Dc); 
        corrcoef_ext_c_theory(j1) = C3(1,2)/(sqrt(C3(1,1)*C3(2,2)));
        variance_ext_c_theory(j1) = C3(2,2)/(OMIGA*steadyODE(j1,2));
        G3 = C3*expm(Ac');
        lagOneC_theory(j1) = G3(2,2)/C3(2,2);
        
        C4 = lyap(Acle,Dcle);
        corrcoef_ext_cle_theory(j1,j2) = C4(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        variance_ext_cle_theory(j1,j2) = C4(2,2)/(OMIGA*steadyODE(j1,2));
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
    plot(a1_list,corrcoef_ext_cle_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:),'LineStyle','--')
end
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
hold off

figure(2)

hold on
plot(a1_list,variance_ext_c_theory,'k-','Linewidth',3)
for j4 = 1:length(sigmae);
    plot(a1_list,variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--')
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('variacne of gene B','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

figure(3)
hold on
plot(a1_list,lagOneC_theory,'k-','Linewidth',3)
for j5 = 1:length(sigmae);
    plot(a1_list,lagOneCLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:),'LineStyle','--')    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
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
