%-----------------------------------------------------------------------
% PROGRAM NAME: intrinsicExtrinsicToggleSwitch.m
% DESCIPTION:  
%   this program using different method to simulate the stochastic
% process of a toggle swithch model. We use Langevin equation, Chmeical
% Langevin equation, Gillespie to compare the effect of extrinsic noise
% when the system is tuned at criticality
% writen by Shanshan Qin@Tanglab, PKU
% last revised on Dec.2,2015
%----------------------------------------------------------------------
% the toggle swithch model
% deterministic
% du = a1/1+v^n1 - u
% dv = a2/(1+ u^n2) - v
%
% noise strength for the LE part
% extrinsic noise is model as O-U prcess
% dE = -taoe*E*dt + c1*dW(t)
% version one: the extrinsic noise is added to the degradation term
% du = a1/1+v^n1 - (1+E)u
% dv = a2/(1+ u^n2) - (1+E)v

%version two: the extrinsic noise is added to the synthesis rate
% du = a1(1+E)/1+v^n1 - u
% dv = a2(1+E)/(1+ u^n2) - v

%version three: the extrinsic noise is added to the basal level synthesis rate
% du = a1/1+v^n1 - u + E
% dv = a2/(1+ u^n2) - v + E
% and time scale for the extrinsic noise is much slower than the intrinsic time scale

%----------------------------------------------------------

function intrinsicExtrinsicToggleSwitch()
% close all
clear
clc


% a1list = 4:1:18;
% sigmae = 0.1;
% tauc = 100;
% param = [15;10;2;2;0.03;tauc;sigmae];
% OMIGA = 20;        %system size
% NUM = 500;
% Nswi = 50;       %number of running when calculate the switching rate
% sigma = [0.1,0.1];    %for LE
% Yini = [0,10];
% % runTime = max([20*tauc,1e3]);
% runTime = 10*tauc;
%__________________________________________________________________________
%plot the bifurcation diagram
% plotBifur()

%_________________________________________________________________________
%this is the theoretic calculation part
%theoretic calculation
% extriLE_CLE_timescale(param,OMIGA,sigma,sigmae,Yini);
% extriLE_CLE_noiseAmpli(param,OMIGA,sigma,sigmae,Yini);

%_________________________________________________________________________

%this is the numeric part
% corr_var_LE_CLE_numeri(param,OMIGA,a1_s,a1_e,N_a1,beta,num_run,sigma,sigmae,dt,tspan,Yini,Nsample);
    
%_________________________________________________________________________

%simulation the distribution using CLE method
% simuDistri(param,NUM,a1list,OMIGA)
%_________________________________________________________________________

% plot the distribution of simulation
% plotDistri(a1list) 

%_________________________________________________________________________
%bifurcation diagram
% fixedPoint(param,a1list)
% [timeCLE,YCLE] = CLEtoggle(tspan,Yini,dt,param,OMIGA);
%_________________________________________________________________________
%plot the switching time
plotSwiTime()

%__________________________________________________________________________

%numeric variance, coefficient of correlation, and lag one auto correlation
% a1list = 4:0.5:15;
% varCorrCLEext(param,a1list,NUM,runTime,OMIGA)


%__________________________________________________________________________
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
function [time,Y]=LEtoggleSwitch(param,sigma,tspan,dt,Yini)

%parameters
%param is a 4 by 1 vector,and sigma is a 2 by 1 vector
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
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
    TEMP = [a1*(a0 +(1-a0 )/(1+Xtemp(2)^n1))-Xtemp(1);a2*(a0 + (1-a0)/(1+Xtemp(1)^n2))-Xtemp(2)];
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
    TEMP = [a1*(a0 + (1-a0)/(1+Xtemp(2)^n1))-(1+Xtemp(3))*Xtemp(1);...
        a2*(a0 + (1-a0)/(1+Xtemp(1)^n2))-(1+Xtemp(3))*Xtemp(2);...
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
function [time,Y]=toggleGillespeBasal(totTime,inV,param,OMIGA)

%this is a first-reaction Gillespie method to simulate the modified 
%toggle switch model(with basal level expression rate)

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);

NUM = 10000000;                %total number point
Y(1,:) = inV;
t = 0;
i = 1;
S = [1,-1,0,0;0,0,1,-1];
while(t < totTime && i<= NUM+1)
    a = [a0*a1*OMIGA+OMIGA*a1*(1-a0)/(1+(Y(i,2)/OMIGA)^n1);Y(i,1);a0*a2*OMIGA+OMIGA*a2*(1-a0)/(1+(Y(i,1)/OMIGA)^n2); Y(i,2)];
    [tau,inx] = min(-1./a.*log(rand(4,1)));
    Y(i+1,:) = Y(i,:) + S(:,inx)';  
    t = t + tau;
    time(i)= t;
    i = i+1;
end
time = [0;time'];
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
a0 = param(5);
N = round(totTime/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = inV;
for i1 = 1:N;
    determin = [a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-Xtemp(1);...
        a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-Xtemp(2)];
    stoch = [sqrt(abs(a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs(Xtemp(1))),0,0;...
            0,0,sqrt(a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))),-sqrt(abs(Xtemp(2)))];        
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
%         variance_ext_cle_theory(j1,j2) = C4(2,2)/(OMIGA*steadyODE(j1,2));
        variance_ext_cle_theory(j1,j2) = C4(2,2);
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

% maxTry = 2000;  %the maximum number of try
for i0 = 6:N_a1;
     
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
        if(a1_list(i0) >=12 && a1_list(i0)<=16) % near the transition point, we need selecet data mannually
            plot(t_LE,Y_LE)
            choice = input('the choice:\n','s');
            if(strcmp(choice,'m')) %mannually selection
                [x,~] = ginput(2);
                [~,inx1] = sort(abs(t_LE-x(1)));
                [~,inx2] = sort(abs(t_LE-x(2)));
                dataSelect = Y_LE(inx1(1):inx2(1),1:2);   %select data
            elseif(strcmp(choice,'y'))  %use all the data
                dataSelect = Y_LE;
            else
                continue      %try another running
            end
        else
            dataSelect = Y_LE(:,1:2);
        end
            steadyLE(i0,j0,:) = mean(dataSelect,1);
            variLE(i0,j0,:) = var(dataSelect,0,1);
            [~, lagAuto_corr_LE(i0,j0,:) ,coef_corr_LE(i0,j0)] =covAutoCorr(dataSelect(:,1),dataSelect(:,2),dt);
            j0 = j0 +1;
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
        if(a1_list(i0) >=10 && a1_list(i0)<=16) % near the transition point, we need selecet data mannually
            plot(t_CLE,Y_CLE)
            choice = input('the choice:\n','s');      %input the choice
            if(strcmp(choice,'m')) %mannually selection
                [x,~] = ginput(2);
                [~,inx1] = sort(abs(t_CLE-x(1)));
                [~,inx2] = sort(abs(t_CLE-x(2)));
                dataSelect = Y_CLE(inx1(1):inx2(1),1:2);   %select data
            elseif(strcmp(choice,'y'))  %use all the data
                dataSelect = Y_CLE;
            else
                continue      %try another running
            end
        else
            dataSelect = Y_CLE;
        end
            steadyCLE(i0,j0,:) = mean(dataSelect,1);
            variCLE(i0,j0,:) = var(dataSelect,0,1);
            [~, lagAuto_corr_CLE(i0,j0,:) ,coef_corr_CLE(i0,j0)] =covAutoCorr(dataSelect(:,1),dataSelect(:,2),dt);
            j0 = j0 +1;
    end
    
%CLE with Extrinsic noise
    j0 =1;
    while(j0 <=num_run )
        paramCLEext = [param;beta;sigmae];
        [t_CLEext,Y_CLEext]= CLEtoggleExtrinsic(tspan,Yini_CLE,dt,paramCLEext,OMIGA);
        if(a1_list(i0) >=10 && a1_list(i0)<=16) % near the transition point, we need selecet data mannually
            plot(t_CLEext,Y_CLEext)
            choice = input('the choice:\n','s');
            if(strcmp(choice,'m')) %mannually selection
                [x,~] = ginput(2);
                [~,inx1] = sort(abs(t_CLEext-x(1)));
                [~,inx2] = sort(abs(t_CLEext-x(2)));
                dataSelect = Y_CLEext(inx1(1):inx2(1),1:2);   %select data
            elseif(strcmp(choice,'y'))  %use all the data
                dataSelect = Y_CLEext(:,1:2);
            else
                continue      %try another running
            end
        else
            dataSelect = Y_CLEext(:,1:2);
        end
            steadyCLEext(i0,j0,:) = mean(dataSelect,1);
            variCLEext(i0,j0,:) = var(dataSelect,0,1);
            [~, lagAuto_corr_CLEext(i0,j0,:) ,coef_corr_CLEext(i0,j0)] =covAutoCorr(dataSelect(:,1),dataSelect(:,2),dt);
            j0 = j0 +1;
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

save('numeric_LE-CLE-1211.mat')

figure(1)
hold on
errorbar([a1_list,a1_list],mean_LE,std_LE,'ro-','LineWidth',2)
errorbar([a1_list,a1_list],mean_CLE,std_CLE,'k^-','LineWidth',2)
errorbar([a1_list,a1_list],mean_CLEext,std_CLEext,'k^-','LineWidth',2)
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',2)
hold off


%correlation of fluction
figure(2)
hold on
errorbar(a1_list,mean_corrcoef_LE,std_corrcoef_LE,'ro-','LineWidth',2)
errorbar(a1_list,mean_corrcoef_CLE,std_corrcoef_CLE,'k^-','LineWidth',2)
errorbar(a1_list,mean_corrcoef_CLEext,std_corrcoef_CLEext,'k^-','LineWidth',2)
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',2)
hold off

%variance
figure(3)
errorbar([a1_list,a1_list],mean_vari_LE,std_vari_LE,'r-o','LineWidth',2,'MarkerSize',6)
hold on
errorbar([a1_list,a1_list],mean_vari_CLE,std_vari_CLE,'k-^','LineWidth',2,'MarkerSize',6)
errorbar([a1_list,a1_list],mean_vari_CLEext,std_vari_CLEext,'k-^','LineWidth',2,'MarkerSize',6)
xlabel('a1','FontSize',24,'FontWeight','Bold')
ylabel('variance','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',2)
hold off

%lag auto correlation coefficient
figure(4)
errorbar([a1_list,a1_list],mean_lagLE,std_LogLE,'r-o','LineWidth',2,'MarkerSize',6)
hold on
errorbar([a1_list,a1_list],mean_lagCLE,std_logCLE,'k-^','LineWidth',2,'MarkerSize',6)
errorbar([a1_list,a1_list],mean_lagCLEext,std_logCLEext,'k-^','LineWidth',2,'MarkerSize',6)
ylabel('lag one autocorrelation','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',2)
hold off

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
function plotBifur()
    fileA = 'toggSwiBifurAll_A.csv';
    fileB = 'toggSwiBifurAll_B.csv';
    dataA = csvread(fileA);
    dataA(dataA==0)=nan;
    dataB = csvread(fileB);
    dataB(dataB==0)=nan;
    plot(dataB(:,1),20*dataB(:,2:4))
    
end
function plotSwiTime()
foldername = 'allST_N20_t100_s0.1_synth';
fileName = 'N20_t100_s0.1_';
x0 = 12:0.5:15;
t_swiAll = [];
for i = 1:20;
        load(fullfile(foldername,[fileName,num2str(i-1),'.mat']));
        t_swiAll = [t_swiAll;t_swi_h2l];
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
legend('intrinsic','\sigma_{E}=0.1,\tau_{c} = 10^{3}','\sigma_{E}=0.25, \tau_{c} = 10^{3}')
legend('CLE-\tau_{E} = 10,\sigma_{E}=0.1','CLE-\tau_{E} = 10^{3},\sigma_{E}=0.1','Glp-\tau_{E} = 10^{3},\sigma_{E}=0.25')

end
