%-----------------------------------------------------------------------
% PROGRAM NAME: intriExtrinMutualActive.m
% DESCIPTION:  
%   this program using different method to simulate the stochastic
% process of a double acitve model. We use Langevin equation, Chmeical
% Langevin equation, Gillespie to compare the effect of extrinsic noise
% when the system is tuned at criticality
% writen by Shanshan Qin@Tanglab, PKU
% last revised on Nov.16,2015
%----------------------------------------------------------------------
% the toggle swithch model
% deterministic
% du = a0 + (1-a0)*v^n1/(K1^n1+v^n1) - u
% dv = a0 + (1-a0)*u^n2/(K2^n+u^n2) - v
%parameters
% K1,K2,n1,n2,a0
% noise strength for the LE part
% extrinsic noise is model as O-U prcess
% dE = -taoe*E*dt + c1*dW(t)

%----------------------------------------------------------

function intriExtrinMutualActive()
close all
clear
clc

param = [0.1;0.5;3;3;0.05];
sigma = [0.01;0.01];      %noise amplitude
sigmae = 0.001;
OMIGA = 200;        %system size
dt = 0.05;          %time step
tspan = 200;        %integration time
Yini = [1;1];      %initial value

%scan parameter
K1_s = 0.01;
K1_e = 0.7;
N_a1 = 40;
beta = 1/100;
num_run = 20;  %20 times repeat
Nsample = 3500;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%this is the theoretic calculation part
%theoretic calculation
extriLE_CLE_timescale(param,OMIGA,sigma,sigmae,tspan,Yini);
% extriLE_CLE_noiseAmpli(param,OMIGA,sigma,beta,tspan,Yini);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%this is the numeric part
% corr_var_LE_CLE_numeri(param,OMIGA,K1_s,K1_e,N_a1,beta,num_run,sigma,sigmae,dt,tspan,Yini,Nsample);
    
       



%save all the data
% save('allData-11-24-1.mat')
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%CLE
% [timeCLE,YCLE] = CLEtoggle(tspan,Yini,dt,param,OMIGA);

%CLE with extrinsic noise
% beta = 1/10;
% sigmae = 10;
% paramCLEExt = [param;beta;sigmae];
% iniCLEExt = OMIGA*Yini;
% [timeCLEExt,YCLEExt] = CLEtoggleExtrinsic(tspan,iniCLEExt,dt,paramCLEExt,OMIGA);

%Gillepie model
% iniGlp = round(OMIGA*Yini');
% ParamGlp = [param;1;1];
% toggleGillespe(tspan,iniGlp,ParamGlp,OMIGA)

%Gillepie model with extrinsic noise
% iniGlp = round(OMIGA*Yini');
% ParamGlp = [param;1/10;0.005];
% [t_GE,Y_GE]=toggleGillespeExt(tspan,iniGlp,ParamGlp,OMIGA);
% plot(t_GE,Y_GE)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%plot the results
%check the stability of stochastic simulation


end

function [t,y]=odeRun(param,totTime,Yini)

tspan = [0 totTime];
[t,y] = ode45(@mutualActi,tspan,Yini,[],param);
% plot(t,y,'LineWidth',2)
% legend('gene1','gene2')

end

function dydt= mutualActi(t,y,parameter)
    K1 = parameter(1);          
    K2 = parameter(2);          
    n1 = parameter(3);        %hill coefficient
    n2 = parameter(4);       %hill coefficient
    a0 = parameter(5);      %basal trate

dydt = zeros(size(y));
u = y(1);
v = y(2);
dudt = a0 + (1-a0)*v^n1/(K1^n1+v^n1) - u;
dvdt = a0 + (1-a0)*u^n2/(K2^n1+u^n2) - v;
dydt = [dudt;dvdt];
end

function [time,Y]=LEtoggleSwitch(param,sigma,tspan,dt,Yini)

%parameters
%param is a 4 by 1 vector,and sigma is a 2 by 1 vector
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
sigma1 = sigma(1);
sigma2 = sigma(2);
% a1 = 5.8;
% a2 = 10;
% n1 = 2;
% n2 = 2;
N = round(tspan/dt);       %total simulation steps
X = zeros(N,2);
% dW = sqrt(dt)*randn(N,2);  %a N by 2 matrix with normal distribution elements
% W = cumsum(dW);            %cummulative add of dW
Xtemp = Yini;
dT = [dt;dt];
SIGMA = [sigma1,0;0,sigma2];
for i1 = 1:N;
    TEMP = [a1/(1+Xtemp(2)^n1)-Xtemp(1);a2/(1+Xtemp(1)^n2)-Xtemp(2)];
    Xtemp = Xtemp + TEMP.*dT + SIGMA*(sqrt(dT).*randn(2,1));
    X(i1,:) = Xtemp;
end
% [variance, lagAutocorrelation ,corrFluc] =covAutoCorr( X(:,1),X(:,2) );
%plot the stochastic trajactories
time = (0:dt:N*dt)';
Y = [Yini';X];
% plot(time,Y,'LineWidth',2)
% legend('gene1','gene2')
% xlabel('time','FontSize',24,'FontName','calibri')
% ylabel('gene expression','FontSize',24,'FontName','calibri')

% using deterministic method to calculate 

%usign Euler method to calculate
% XtempEu = Xzero;
% for i2 = 1:N;
%     stepTemp = [a1/(1+XtempEu(2)^n1)-XtempEu(1);a2/(1+XtempEu(1)^n2)-XtempEu(2)];
%     XtempEu = XtempEu + stepTemp.*dT;
%     Xeuler(i2,:) = XtempEu;
% end
% plot((0:dt:T)',[Xzero';Xeuler])
% hold off
end
function [time,Y]=LEtoggleExtrinsic(param,sigma,tspan,dt,Yini)

%parameters
%param is a 4 by 1 vector,and sigma is a 2 by 1 vector
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
beta = param(5);      %the time scale of extrinsic noise
sigma1 = sigma(1);
sigma2 = sigma(2);
sigmae = sigma(3);    %the white noise term of the extrinsic noise

N = round(tspan/dt);       %total simulation steps
X = zeros(N,3);
Xtemp = [Yini;0];
dT = [dt;dt;dt];
SIGMA = [sigma1,0,0;0,sigma2,0;0,0,sigmae];
for i1 = 1:N;
%     TEMP = [a1/(1+Xtemp(2)^n1)-Xtemp(1)+Xtemp(3);...
%         a2/(1+Xtemp(1)^n2)-Xtemp(2)+Xtemp(3);...
%         -beta*Xtemp(3)];
    TEMP = [a1/(1+Xtemp(2)^n1)-(1+Xtemp(3))*Xtemp(1);...
        a2/(1+Xtemp(1)^n2)-(1+Xtemp(3))*Xtemp(2);...
        -beta*Xtemp(3)];
    
    Xtemp = Xtemp + TEMP.*dT + SIGMA*(sqrt(dT).*randn(3,1));
    X(i1,:) = Xtemp';
end
time = (0:dt:N*dt)';
Y = [[Yini;0]';X];
% plot(time,Y,'LineWidth',2)
% legend('gene1','gene2')
% xlabel('time','FontSize',24,'FontName','calibri')
% ylabel('gene expression','FontSize',24,'FontName','calibri')

end
function [time,Y]=toggleGillespe(totTime,inV,param,OMIGA)

%basic Gillespie alogrithm

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
beta1 = param(5);
beta2 = param(6);

NUM = 1000000;                %total number point
Y(1,:) = inV;
t = 0;
i = 1;

while(t < totTime && i<= NUM+1)
    a = [OMIGA*a1/(1+(Y(i,2)/OMIGA)^n1); beta1*Y(i,1);...
        OMIGA*a2/(1+(Y(i,1)/OMIGA).^n2); beta2*Y(i,2)];
    a0 = sum(a);   
    t = t + 1/a0*log(1/rand);
    time(i)= t;
    reference = rand*a0;
    if(reference<a(1))
        Y(i+1,1) = Y(i,1) + 1;
        Y(i+1,2) = Y(i,2);
        i = i +1;
    elseif(reference > a(1) && (reference <= sum(a(1:2))))
        Y(i+1,1) = Y(i,1) -1;
        Y(i+1,2) = Y(i,2);
        i = i+1;
    elseif(reference > sum(a(1:2)) && (reference <= sum(a(1:3))))
        Y(i+1,2) = Y(i,2) + 1;
        Y(i+1,1) = Y(i,1);
        i =i +1;
    elseif(reference > sum(a(1:3)) && (reference <= sum(a(1:4))))
        Y(i+1,2) = Y(i,2) - 1;
        Y(i+1,1) = Y(i,1);
        i =i +1;       
    end
end
time = [0;time'];
% plot(time,Y)
% hold on
%--------------------------------------------------------------------------
%ode
% tspan = [0 time(length(time))];
% [t,y] = ode45(@GapGeneNetwork,tspan,Y(1,:)',[],Km,tao,n,OMIGA);
% plot(t,y)
% xlabel('Time(a.u)','FontSize',30,'FontName','Times New Roman')
% ylabel('Gene Expression(a.u)','FontSize',30,'FontName','Times New Roman')
% set(gca,'FontSize',20)
end
function [time,Y]=toggleGillespeExt(totTime,inV,param,OMIGA)

%basic Gillespie alogrithm

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
beta = param(5);
sigmae = param(6);

NUM = 1000000;                %total number point
Y(1,:) = inV;
YE(1) = 0;
t = 0;
i = 1;

while(t < totTime && i<= NUM+1)
    a = [OMIGA*a1/(1+(Y(i,2)/OMIGA)^n1); (1+YE(i))*Y(i,1);...
        OMIGA*a2/(1+(Y(i,1)/OMIGA).^n2); (1+YE(i))*Y(i,2)];
    a0 = sum(a);
    dt = 1/a0*log(1/rand);
    t = t + dt;
    time(i)= t;
    reference = rand*a0;
    if(reference<a(1))
        Y(i+1,1) = Y(i,1) + 1;
        Y(i+1,2) = Y(i,2);
        i = i +1;
    elseif(reference > a(1) && (reference <= sum(a(1:2))))
        Y(i+1,1) = Y(i,1) -1;
        Y(i+1,2) = Y(i,2);
        i = i+1;
    elseif(reference > sum(a(1:2)) && (reference <= sum(a(1:3))))
        Y(i+1,2) = Y(i,2) + 1;
        Y(i+1,1) = Y(i,1);
        i =i +1;
    elseif(reference > sum(a(1:3)) && (reference <= sum(a(1:4))))
        Y(i+1,2) = Y(i,2) - 1;
        Y(i+1,1) = Y(i,1);
        i =i +1;       
    end
    
    %update the extrinsic noise
    mu = exp(-beta*dt);
    Dext = sqrt(abs(sigmae/(2*beta)*(1-mu^2)));
    YE(i) = YE(i-1)*mu + Dext*randn;
end
time = [0;time'];
% plot(time,Y)
% hold on
% plot(time,YE)
% hold off
% %--------------------------------------------------------------------------
%ode
% tspan = [0 time(length(time))];
% [t,y] = ode45(@GapGeneNetwork,tspan,Y(1,:)',[],Km,tao,n,OMIGA);
% plot(t,y)
% xlabel('Time(a.u)','FontSize',30,'FontName','Times New Roman')
% ylabel('Gene Expression(a.u)','FontSize',30,'FontName','Times New Roman')
% set(gca,'FontSize',20)
end
function [time,Y] = CLEtoggle(totTime,inV,dt,param,OMIGA)
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);

N = round(totTime/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = inV;
for i1 = 1:N;
    determin = [a1*OMIGA/(1+(Xtemp(2)/OMIGA)^n1)-Xtemp(1);...
        a2*OMIGA/(1+(Xtemp(1)/OMIGA)^n2)-Xtemp(2)];
    stoch = [sqrt(abs(a1*OMIGA/(1+(Xtemp(2)/OMIGA)^n1))),-sqrt(abs(Xtemp(1))),0,0;...
            0,0,sqrt(a2*OMIGA/(1+(Xtemp(1)/OMIGA)^n2)),-sqrt(abs(Xtemp(2)))];        
    Xtemp = Xtemp + determin*dt + stoch*randn(4,1)*sqrt(dt);
    X(i1,:) = Xtemp;
end
% [variance, lagAutocorrelation ,corrFluc] =covAutoCorr( X(:,1),X(:,2) );
%plot the stochastic trajactories
time = (0:dt:N*dt)';
Y = [inV';X];
% plot(time,Y,'LineWidth',2)
% legend('gene1','gene2')
% xlabel('time','FontSize',24,'FontName','calibri')
% ylabel('gene expression','FontSize',24,'FontName','calibri')

end
function [time,Y] = CLEtoggleExtrinsic(totTime,inV,dt,param,OMIGA)
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
beta = param(5);
sigmae = param(6);  %the amplitude of extrinsic noise

N = round(totTime/dt);       %total simulation steps
X = zeros(N,3);
Xtemp = [inV;0];
for i1 = 1:N;
    determin = [a1*OMIGA/(1+(Xtemp(2)/OMIGA)^n1)-(1+Xtemp(3))*Xtemp(1);...
        a2*OMIGA/(1+(Xtemp(1)/OMIGA)^n2)-(1+Xtemp(3))*Xtemp(2);...
        -beta*Xtemp(3)];
    stoch = [sqrt(a1*OMIGA/(1+(Xtemp(2)/OMIGA)^n1)),-sqrt(abs((1+Xtemp(3))*Xtemp(1))),0,0,0;...
            0,0,sqrt(a2*OMIGA/(1+(Xtemp(1)/OMIGA)^n2)),-sqrt(abs((1+Xtemp(3))*Xtemp(2))),0;...
            0,0,0,0,sigmae];

    Xtemp = Xtemp + determin*dt + stoch*randn(5,1)*sqrt(dt);
    X(i1,:) = Xtemp;
end

time = (0:dt:N*dt)';
Y = [[inV;0]';X];
% plot(time,Y,'LineWidth',2)
% legend('gene1','gene2')
% xlabel('time','FontSize',24,'FontName','calibri')
% ylabel('gene expression','FontSize',24,'FontName','calibri')

end
function crossCorr= crossCorrCoef(X,Y,tau)
%X and Y are vector that have the same size
%dt is the time step size
%tau is the lag time
X_res = X - mean(X);
Y_res = Y - mean(Y);
X_std = std(X_res);
Y_std = std(Y_res);
LEN = length(X);
inxList = 1:1:LEN-abs(tau);
if(tau>=0)
    crossCorr = sum(X_res(tau+inxList).*Y_res(inxList))/(LEN-abs(tau))/sqrt(X_std*Y_std);
else
    crossCorr = sum(Y_res(abs(tau)+inxList).*X_res(inxList))/(LEN-abs(tau))/sqrt(X_std*Y_std);
end

end
function [time,Y]=LEindependent(param,sigma,tspan,dt,Yini)

% this function using Langevin equation to simulate two gene without
% interaction to show the effect of extrinsic noise 
% the simple model is 
% dA/dt = a1 - b1*A;
% dB/dt = a2 - b2*B;
%parameters
%param is a 4 by 1 vector,and sigma is a 2 by 1 vector
a1 = param(1);
a2 = param(2);
b1 = param(3);
b2 = param(4);
beta = param(5);      %the time scale of extrinsic noise
sigma1 = sigma(1);
sigma2 = sigma(2);
sigmae = sigma(3);    %the white noise term of the extrinsic noise

N = round(tspan/dt);       %total simulation steps
X = zeros(N,3);
Xtemp = [Yini;0];
dT = [dt;dt;dt];
SIGMA = [sigma1,0,0;0,sigma2,0;0,0,sigmae];
for i1 = 1:N;
%     TEMP = [a1/(1+Xtemp(2)^n1)-Xtemp(1)+Xtemp(3);...
%         a2/(1+Xtemp(1)^n2)-Xtemp(2)+Xtemp(3);...
%         -beta*Xtemp(3)];
    TEMP = [a1/(1+Xtemp(2)^n1)-(1+Xtemp(3))*Xtemp(1);...
        a2/(1+Xtemp(1)^n2)-(1+Xtemp(3))*Xtemp(2);...
        -beta*Xtemp(3)];
    
    Xtemp = Xtemp + TEMP.*dT + SIGMA*(sqrt(dT).*randn(3,1));
    X(i1,:) = Xtemp';
end
time = (0:dt:N*dt)';
Y = [[Yini;0]';X];
% plot(time,Y,'LineWidth',2)
% legend('gene1','gene2')
% xlabel('time','FontSize',24,'FontName','calibri')
% ylabel('gene expression','FontSize',24,'FontName','calibri')

end
function [A,Ae,Ac,Ace] = jacMatrx(param,vec,OMIGA)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
syms x y E K1 K2 n1 n2 a0 b N
J1 = jacobian([a0 + (1-a0)*y^n1/(K1^n1+y^n1) - x;a0 + (1-a0)*x^n2/(K2^n2+x^n2) - y],[x,y]);
A = double(subs(J1,{x,y,K1 K2 n1 n2 a0},{vec(1),vec(2),param(1),param(2),param(3),param(4),param(5)}));
% J2 = jacobian([a1/(1+y^n1) - x + E;a2/(1+x^n2) - y + E;-b*E],[x y E]);
J2 = jacobian([a0 + (1-a0)*y^n1/(K1^n1+y^n1) - (1+E)*x;a0 + (1-a0)*x^n2/(K2^n2+x^n2) - (1+E)*y;-b*E],[x y E]);
% J2 = jacobian([a1/(1+E+y^n1) - x;a2/(1+E+x^n2) - y;-b*E],[x y E]);
% J2 = jacobian([a1*(1+E)/(1+y^n1) - x;a2*(1+E)/(1+x^n2) - y;-b*E],[x y E]);
Ae = double(subs(J2,{x,y,E,K1 K2 n1 n2 a0 b},{vec(1),vec(2),0,param(1),param(2),param(3),param(4),param(5),param(6)}));

J3 = jacobian([a0*N + (1-a0)*N*y^n1/((K1*N)^n1+y^n1) - x;a0*N + (1-a0)*N*x^n2/((K2*N)^n2+x^n2) - y],[x y]);
Ac = double(subs(J3,{x,y,K1 K2 n1 n2 a0 N},{vec(1)*OMIGA,vec(2)*OMIGA,param(1),param(2),param(3),param(4),param(5),OMIGA}));
% J4 = jacobian([a1*N/(1+E+(y/N)^n1) - x;a2*N/(1+E+(x/N)^n2) - y;-b*E],[x y E]);
J4 = jacobian([a0*N + (1-a0)*N*y^n1/((K1*N)^n1+y^n1) - (1+E)*x;a0*N + (1-a0)*N*x^n2/((K2*N)^n2+x^n2) - (1+E)*y;-b*E],[x y E]);
% J4 = jacobian([a1*(1+E)*N/(1+(y/N)^n1) - x;a2*(1+E)*N/(1+(x/N)^n2) - y;-b*E],[x y E]);
Ace = double(subs(J4,{x,y,E,K1 K2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,0,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

end
function extriLE_CLE_timescale(param,OMIGA,sigma,sigmae,tspan,Yini)
%this function comparing the correlation coefficient, variance, and lag one
%auto correlation for LE and CLE with extrinsic noise
% beta define different time scale of extrinsic noise

K1_s = 0.01;
K1_e = 0.7;
N_K1 = 50;
K2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);

beta = [1,1/10,1/50,1/500];
corrcoef_theory = zeros(N_K1,1);
corrcoef_ext_theory = zeros(N_K1,length(beta));
variance_theory = zeros(N_K1,1);
variance_ext_theory = zeros(N_K1,length(beta));
corrcoef_ext_cle_theory = zeros(N_K1,1);
variance_ext_cle_theory = zeros(N_K1,length(beta));
lagOneLe_theory = zeros(N_K1,1);
lagOneLEext_theory = zeros(N_K1,length(beta));
lagOneCLEext_theory = zeros(N_K1,length(beta));
variance_ext_c_theory = zeros(N_K1,1);
corrcoef_ext_c_theory = zeros(N_K1,1);
lagOneC_theory = zeros(N_K1,1);
for j1 = 1:N_K1;
    K1_list(j1) = K1_s + (K1_e - K1_s)*j1/N_K1;
    param(1)= K1_list(j1);
    [td,Yd]=odeRun(param,tspan,Yini);
    steadyODE(j1,:)=Yd(end,:);
    for j2 = 1:length(beta)
        [A,Ae,Ac,Acle] = jacMatrx([param;beta(j2)],steadyODE(j1,:),OMIGA);
        D = [sigma(1)^2,0;0,sigma(2)^2];
        De = [sigma(1)^2,0,0;0,sigma(2)^2,0;0,0,sigmae^2];
        Dc = [a0*OMIGA + (1-a0)*OMIGA*steadyODE(j1,2)^n1/((K1_list(j1)*OMIGA)^n1+steadyODE(j1,2)^n1)+steadyODE(j1,1),0;...
                0,a0*OMIGA + (1-a0)*OMIGA*steadyODE(j1,1)^n2/((param(2)*OMIGA)^n2+steadyODE(j1,1)^n2)+steadyODE(j1,2)];
        Dcle = [a0*OMIGA + (1-a0)*OMIGA*steadyODE(j1,2)^n1/((K1_list(j1)*OMIGA)^n1+steadyODE(j1,2)^n1)+steadyODE(j1,1),0,0;...
                0,a0*OMIGA + (1-a0)*OMIGA*steadyODE(j1,1)^n2/((param(2)*OMIGA)^n2+steadyODE(j1,1)^n2)+steadyODE(j1,2),0;0,0,sigmae^2];
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
%plot the theoretic results
colorSet = [ones(4,1)*60/256,[0;0.3;0.6;0.9],ones(4,1)*189/256];
% colorSet = [[0;0.3;0.6;0.9],zeros(4,2)];
figure(1)
plot(K1_list,corrcoef_theory,'r-','Linewidth',3)
hold on
plot(K1_list,corrcoef_ext_c_theory,'k-','Linewidth',3)
for j3 = 1:length(beta);
    plot(K1_list,corrcoef_ext_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:))
    plot(K1_list,corrcoef_ext_cle_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:),'LineStyle','--')
end
legend('LE','CLE','LE-extr-1','CLE-extr-1','LE-extr-10','CLE-extr-10','LE-extr-50','CLE-extr-50','LE-extr-500','CLE-extr-500')
xlabel('K1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
hold off

figure(2)
plot(K1_list,variance_theory,'r-','Linewidth',3)
hold on
plot(K1_list,variance_ext_c_theory,'k-','Linewidth',3)
for j4 = 1:length(beta);
    plot(K1_list,variance_ext_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:))
    plot(K1_list,variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--')
    
end
hold off
legend('LE','CLE','LE-extr-1','CLE-extr-1','LE-extr-10','CLE-extr-10','LE-extr-50','CLE-extr-50','LE-extr-500','CLE-extr-500')
xlabel('K1','FontSize',24,'FontWeight','Bold')
ylabel('variacne of gene B','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

figure(3)
plot(K1_list,lagOneLe_theory,'r-','LineWidth',3)
hold on
plot(K1_list,lagOneC_theory,'k-','Linewidth',3)
for j5 = 1:length(beta);
    plot(K1_list,lagOneLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:))
    plot(K1_list,lagOneCLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:),'LineStyle','--')    
end
hold off
legend('LE','CLE','LE-extr-1','CLE-extr-1','LE-extr-10','CLE-extr-10','LE-extr-50','CLE-extr-50','LE-extr-500','CLE-extr-500')
xlabel('K1','FontSize',24,'FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
end
function extriLE_CLE_noiseAmpli(param,OMIGA,sigma,b,tspan,Yini)
%this function comparing the correlation coefficient, variance, and lag one
%auto correlation for LE and CLE with extrinsic noise
% sigmae define diffent extrinsic noise amplitude

a1_s = 2;
a1_e = 50;
N_a1 = 50;


sigmae = [1e-3,1e-2,0.05,0.1];
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
    [td,Yd]=odeRun(param,tspan,Yini);
    steadyODE(j1,:)=Yd(end,:);
    for j2 = 1:length(sigmae)
        [A,Ae,Ac,Acle] = jacMatrx([param;b],steadyODE(j1,:),OMIGA);
        D = [sigma(1)^2,0;0,sigma(2)^2];
        De = [sigma(1)^2,0,0;0,sigma(2)^2,0;0,0,sigmae(j2)^2];
        Dc = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+steadyODE(j1,1),0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+steadyODE(j1,2)];
        Dcle = [a1_list(j1)*OMIGA/(1+steadyODE(j1,2)^param(3))+steadyODE(j1,1),0,0;...
                0,param(2)*OMIGA/(1+steadyODE(j1,1)^param(4))+steadyODE(j1,2),0;0,0,sigmae(j2)^2];
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
%plot the theoretic results
colorSet = [ones(4,1)*60/256,[0;0.3;0.7;1],ones(4,1)*189/256];
% colorSet = [[0;0.3;0.6;0.9],zeros(4,2)];
figure(1)
plot(a1_list,corrcoef_theory,'r-','Linewidth',3)
hold on
plot(a1_list,corrcoef_ext_c_theory,'k-','Linewidth',3)
for j3 = 1:length(sigmae);
    plot(a1_list,corrcoef_ext_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:))
    plot(a1_list,corrcoef_ext_cle_theory(:,j3),'Linewidth',3,'Color',colorSet(j3,:),'LineStyle','--')
end
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
hold off

figure(2)
plot(a1_list,variance_theory,'r-','Linewidth',3)
hold on
plot(a1_list,variance_ext_c_theory,'k-','Linewidth',3)
for j4 = 1:length(sigmae);
    plot(a1_list,variance_ext_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:))
    plot(a1_list,variance_ext_cle_theory(:,j4),'Linewidth',3,'Color',colorSet(j4,:),'LineStyle','--')
    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('variacne of gene B','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)

figure(3)
plot(a1_list,lagOneLe_theory,'r-','LineWidth',3)
hold on
plot(a1_list,lagOneC_theory,'k-','Linewidth',3)
for j5 = 1:length(sigmae);
    plot(a1_list,lagOneLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:))
    plot(a1_list,lagOneCLEext_theory(:,j5),'Linewidth',3,'Color',colorSet(j5,:),'LineStyle','--')    
end
hold off
legend('LE','CLE','LE-extr-0.001','CLE-extr-0.001','LE-extr-0.01','CLE-extr-0.01','LE-extr-0.05','CLE-extr-0.05','LE-extr-0.1','CLE-extr-0.1')
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',2)
end
function corr_var_LE_CLE_numeri(param,OMIGA,a1_s,a1_e,N_a1,beta,num_run,sigma,sigmae,dt,tspan,Yini,Nsample)
a1_list = zeros(N_a1,1);
steadyODE = zeros(N_a1,2);
steadyLE = zeros(N_a1,num_run,2);
% steadyLEExt = zeros(N_a1,num_run,2);
steadyCLE = zeros(N_a1,num_run,2);
steadyCLEext = zeros(N_a1,num_run,2);

variLE = zeros(N_a1,num_run,2);
% variLE_Ext = zeros(N_a1,num_run,2);
variCLE = zeros(N_a1,num_run,2);
variCLEext = zeros(N_a1,num_run,2);

coef_corr_LE = zeros(N_a1,num_run);
% coef_corr_LE_ext = zeros(N_a1,num_run);
coef_corr_CLE =  zeros(N_a1,num_run);
coef_corr_CLEext = zeros(N_a1,num_run);

lagAuto_corr_LE = zeros(N_a1,num_run,2);
% lagAuto_corr_LE_Ext = zeros(N_a1,num_run,2);
lagAuto_corr_CLE = zeros(N_a1,num_run,2);
lagAuto_corr_CLEext = zeros(N_a1,num_run,2);
% Nsample = 3500;  %number of data used
for i0 = 1:N_a1;
     
% deterministic model
    a1_list(i0) = a1_s + (a1_e - a1_s)*i0/N_a1;
    param(1)= a1_list(i0);
    [td,Yd]=odeRun(param,tspan,Yini);
    steadyODE(i0,:)=Yd(end,:);
    
    [td,Y_ref]=odeRun(param,tspan,[steadyODE(i0,2);steadyODE(i0,1)]);
    REF = Y_ref(end,:)*OMIGA;
    j0 = 1;

% simulte wiht LE equation
    Yini_LE = transpose(steadyODE(i0,:));
    while(j0<=num_run);
        [t_LE,Y_LE]=LEtoggleSwitch(param,sigma,tspan,dt,Yini_LE);
        if(max(Y_LE(:,1)) - min(Y_LE(:,1))<15)
            steadyLE(i0,j0,:) = mean(Y_LE,1);
            variLE(i0,j0,:) = var(Y_LE,0,1);
            [~, lagAuto_corr_LE(i0,j0,:) ,coef_corr_LE(i0,j0)] =covAutoCorr(Y_LE(:,1),Y_LE(:,2),Nsample,dt);
            j0 = j0 +1;
        end
    end
% LE with extrinsic noise
%     j0 =1;
%     while(j0 <=num_run )
%         SigmaExt = [sigma;sigmae];
%         beta = 1/50;  %the fluctuation time scale of extrinsic noise
%         Yini_LEExt = transpose(steadyODE(i0,:));
%         ParamExt = [param;beta];
%         [t_LEExt,Y_LEExt]=LEtoggleExtrinsic(ParamExt,SigmaExt,tspan,dt,Yini_LEExt);
%         if(max(Y_LEExt(:,1))- min(Y_LEExt(:,1))<15)
%             steadyLEExt(i0,j0,:) = mean(Y_LEExt(:,1:2),1);
%             variLE_Ext(i0,j0,:) = var(Y_LEExt(:,1:2),0,1);
%             [~,lagAuto_corr_LE_Ext(i0,j0,:),coef_corr_LE_ext(i0,j0)] = covAutoCorr(Y_LEExt(:,1),Y_LEExt(:,2),Nsample,dt);
%             j0 = j0 +1;
%         end
%     end

% CLE
    j0 =1;
    Yini_CLE = transpose(round(OMIGA*steadyODE(i0,:)));
    while(j0 <=num_run )
        [t_CLE,Y_CLE]=CLEtoggle(tspan,Yini_CLE,dt,param,OMIGA);
        if(abs(Y_CLE(end,1) - REF(1))>100 && abs(Y_CLE(end,2) - REF(2))>100)
            steadyCLE(i0,j0,:) = mean(Y_CLE(:,1:2)/OMIGA,1);
            variCLE(i0,j0,:) = var(Y_CLE(:,1:2)/OMIGA,0,1);
            [~,lagAuto_corr_CLE(i0,j0,:),coef_corr_CLE(i0,j0)] = covAutoCorr(Y_CLE(:,1),Y_CLE(:,2),Nsample,dt);
            j0 = j0 +1;
        end
    end
    
%CLE with Extrinsic noise
    j0 =1;
    while(j0 <=num_run )
        paramCLEext = [param;beta;sigmae];
        [t_CLEext,Y_CLEext]= CLEtoggleExtrinsic(tspan,Yini_CLE,dt,paramCLEext,OMIGA);
%         plot(t_CLEext,Y_CLEext)
        if(abs(Y_CLEext(end,1) - REF(1))>100 && abs(Y_CLEext(end,2) - REF(2))>100)
            steadyCLEext(i0,j0,:) = mean(Y_CLEext(:,1:2)/OMIGA,1);
            variCLEext(i0,j0,:) = var(Y_CLEext(:,1:2)/OMIGA,0,1);
            [~,lagAuto_corr_CLEext(i0,j0,:),coef_corr_CLEext(i0,j0)] = covAutoCorr(Y_CLEext(:,1),Y_CLEext(:,2),Nsample,dt);
            j0 = j0 +1;
        end
    end 
end
le = mean(steadyLE,2);
mean_LE = [le(:,:,1),le(:,:,2)];
sle = std(steadyLE,0,2);
std_LE = [sle(:,:,1),sle(:,:,2)];

cle = mean(steadyCLE,2);
mean_CLE = [cle(:,:,1),cle(:,:,2)];
scle = std(steadyCLE,0,2);
std_CLE = [scle(:,:,1),scle(:,:,2)];

clex = mean(steadyCLEext,2);
mean_CLEext = [clex(:,:,1),clex(:,:,2)];
sclex = std(steadyCLEext,0,2);
std_CLEext = [sclex(:,:,1),sclex(:,:,2)];
% 
% ge = mean(steadyGE,2);
% mean_GE = [ge(:,:,1),ge(:,:,2)];
% sge = std(steadyGE,0,2);
% std_GE = [sge(:,:,1),sge(:,:,2)];

mean_corrcoef_LE = mean(coef_corr_LE,2);
std_corrcoef_LE = std(coef_corr_LE,0,2);
mean_corrcoef_CLE = mean(coef_corr_CLE,2);
std_corrcoef_CLE = std(coef_corr_CLE,0,2);
mean_corrcoef_CLEext = mean(coef_corr_CLEext,2);
std_corrcoef_CLEext = std(coef_corr_CLEext,0,2);
% mean_corrcoef_GP = mean(coef_corr_GP,2);
% std_corrcoef_GP = std(coef_corr_GP,0,2);
% mean_corrcoef_GE = mean(coef_corr_GE,2);
% std_corrcoef_GE = std(coef_corr_GE,0,2);


mle = mean(variLE,2);
mean_vari_LE = [mle(:,:,1),mle(:,:,2)];
stdle = std(variLE,0,2);
std_vari_LE = [stdle(:,:,1),stdle(:,:,2)];

mcle = mean(variCLE,2);
mean_vari_CLE = [mcle(:,:,1),mcle(:,:,2)];
stdcle = std(variCLE,0,2);
std_vari_CLE = [stdcle(:,:,1),stdcle(:,:,2)];

mclex = mean(variCLEext,2);
mean_vari_CLEext = [mclex(:,:,1),mclex(:,:,2)];
stdvari_clex = std(variCLEext,0,2);
std_vari_CLEext = [stdvari_clex(:,:,1),stdvari_clex(:,:,2)];

% mGE = mean(variGE,2);
% mean_vari_GE = [mGE(:,:,1),mGE(:,:,2)];
% stdGE = std(variGE,0,2);
% std_vari_GE = [stdGE(:,:,1),stdGE(:,:,2)];


mlagle = mean(lagAuto_corr_LE,2);
mean_lagLE = [mlagle(:,:,1),mlagle(:,:,2)];
stdlagle = std(lagAuto_corr_LE,0,2);
std_LogLE = [stdlagle(:,:,1),stdlagle(:,:,2)];

mlagcle = mean(lagAuto_corr_CLE,2);
mean_lagCLE= [mlagcle(:,:,1),mlagcle(:,:,1)];
stdlagcle = std(lagAuto_corr_CLE,0,2);
std_logCLE= [stdlagcle(:,:,1),stdlagcle(:,:,2)];

mlagclex = mean(lagAuto_corr_CLEext,2);
mean_lagCLEext= [mlagclex(:,:,1),mlagclex(:,:,1)];
stdlagclex = std(lagAuto_corr_CLEext,0,2);
std_logCLEext= [stdlagclex(:,:,1),stdlagclex(:,:,2)];

% save('numeric_LE-CLE-11-25.mat')

% figure(1)
% hold on
% errorbar([a1_list,a1_list],mean_LE,std_LE,'ro-','LineWidth',2)
% errorbar([a1_list,a1_list],mean_CLE,std_CLE,'k^-','LineWidth',2)
% errorbar([a1_list,a1_list],mean_CLEext,std_CLEext,'k^-','LineWidth',2)
% xlabel('a1','FontSize',24,'FontWeight','Bold')
% ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',20,'LineWidth',2)
% hold off
% 
% 
% %correlation of fluction
% figure(2)
% hold on
% errorbar(a1_list,mean_corrcoef_LE,std_corrcoef_LE,'ro-','LineWidth',2)
% errorbar(a1_list,mean_corrcoef_CLE,std_corrcoef_CLE,'k^-','LineWidth',2)
% errorbar(a1_list,mean_corrcoef_CLEext,std_corrcoef_CLEext,'k^-','LineWidth',2)
% xlabel('a1','FontSize',24,'FontWeight','Bold')
% ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',20,'LineWidth',2)
% hold off
% 
% %variance
% figure(3)
% errorbar([a1_list,a1_list],mean_vari_LE,std_vari_LE,'r-o','LineWidth',2,'MarkerSize',6)
% hold on
% errorbar([a1_list,a1_list],mean_vari_CLE,std_vari_CLE,'k-^','LineWidth',2,'MarkerSize',6)
% errorbar([a1_list,a1_list],mean_vari_CLEext,std_vari_CLEext,'k-^','LineWidth',2,'MarkerSize',6)
% xlabel('a1','FontSize',24,'FontWeight','Bold')
% ylabel('variance','FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',20,'LineWidth',2)
% hold off
% 
% %lag auto correlation coefficient
% figure(4)
% errorbar([a1_list,a1_list],mean_lagLE,std_LogLE,'r-o','LineWidth',2,'MarkerSize',6)
% hold on
% errorbar([a1_list,a1_list],mean_lagCLE,std_logCLE,'k-^','LineWidth',2,'MarkerSize',6)
% errorbar([a1_list,a1_list],mean_lagCLEext,std_logCLEext,'k-^','LineWidth',2,'MarkerSize',6)
% ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',20,'LineWidth',2)
% hold off

end