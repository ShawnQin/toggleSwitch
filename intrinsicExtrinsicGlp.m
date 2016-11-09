%-----------------------------------------------------------------------
% PROGRAM NAME: intrinsicExtrinsicGlp.m
% DESCIPTION:  
%   this program using Gillespie to compare the effect of extrinsic noise
% writen by Shanshan Qin@Tanglab, PKU
% last revised on Nov.24,2015
%----------------------------------------------------------------------
% the toggle swithch model
% deterministic
% du = a1*(a0 + (1-a0)/1+v^n1) - u
% dv = a2*(a0 + (1-a0)/(1+ u^n2)) - v
%----------------------------------------------------------

function intrinsicExtrinsicGlp()
close all
clear
clc

a1list = 4:1:18;
sigmae = 0;
tauc = 100;
param = [7;10;2;2;0.03;tauc;sigmae];
OMIGA = 20;        %system size
NUM = 3;
Nswi = 50;       %number of running when calculate the switching rate
sigma = [0.1,0.1];    %for LE
Yini = [0,10];

%_________________________________________________________________________
a1list = 11:1:15;
% switchTime(param,a1list,NUM,OMIGA);
% glpDataStatis()
[time,Y]= GillespeBasalExt(1e3,[20;200],param,OMIGA);
% [corr,variance] = GlpCorrVar(param,OMIGA,sigmae,beta);
N_dist = 1e3;
a1_select = [2,5,10,15,18,20,22,25,30];
[endGlp,endGlpExt] = GlpDistri(a1_select,param,beta,sigmae,OMIGA,N_dist);
save('GillespieEndVaule-12-09.mat','endGlp','endGlpExt','a1_select')


% for i0 = 1:N_a1;    
% %deterministic model
%     a1_list(i0) = a1_s + (a1_e - a1_s)*i0/N_a1;
%     param(1)= a1_list(i0);
%     [td,Yd]=odeRun(param,tspan,Yini);
%     steadyODE(i0,:)=Yd(end,:);
%     for j0 = 1:num_run;
        
%Gillespie algorithum
%         iniGlp = round(OMIGA*[Yd(end,1);Yd(end,2)]);
%         ParamGlp = [param;1;1];
%         [t_GP,Y_GP]=toggleGillespe(tspan,iniGlp,ParamGlp,OMIGA);
%         [f,xi] = ksdensity(Y_GP(:,2));
%         plot(xi,f)
%         steadyGlp(i0,j0,:) = mean(Y_GP,1);
%         variGlp(i0,j0,:) = var(Y_GP,0,1);
%         num_select = round(0.1*length(Y_GP(:,1)));
%         C1 = corrcoef(Y_GP(num_select:end,1),Y_GP(num_select:end,2));
%         coef_corr_GP(i0,j0) = C1(1,2);
        
%Gillespie with extrinsic noise
%         ParamGlp = [param;beta;0.1*SigmaExt];
%         [t_GE,Y_GE]=toggleGillespeExt(tspan,iniGlp,ParamGlp,OMIGA);
%         steadyGE(i0,j0,:) = mean(Y_GE,1);
%         variGE(i0,j0,:) = var(Y_GE,0,1);
%         num_select = round(0.2*length(Y_GE(:,1)));
%         C2 = corrcoef(Y_GE(num_select:end,1),Y_GE(num_select:end,2));
%         coef_corr_GE(i0,j0) = C2(1,2);
       
%     end
% end


% gp = mean(steadyGlp,2);
% mean_Glp = [gp(:,:,1),gp(:,:,2)];
% sgp = std(steadyGlp,0,2);
% std_Glp = [sgp(:,:,1),sgp(:,:,2)];

% ge = mean(steadyGE,2);
% mean_GE = [ge(:,:,1),ge(:,:,2)];
% sge = std(steadyGE,0,2);
% std_GE = [sge(:,:,1),sge(:,:,2)];

% mean_corrcoef_GP = mean(coef_corr_GP,2);
% std_corrcoef_GP = std(coef_corr_GP,0,2);
% mean_corrcoef_GE = mean(coef_corr_GE,2);
% std_corrcoef_GE = std(coef_corr_GE,0,2);

% mGlp = mean(variGlp,2);
% mean_vari_Glp = [mGlp(:,:,1),mGlp(:,:,2)];
% stdGlp = std(variGlp,0,2);
% std_vari_Glp = [stdGlp(:,:,1),stdGlp(:,:,2)];

% mGE = mean(variGE,2);
% mean_vari_GE = [mGE(:,:,1),mGE(:,:,2)];
% stdGE = std(variGE,0,2);
% std_vari_GE = [stdGE(:,:,1),stdGE(:,:,2)];

%save all the data
% save('GlpCorrVar-11-30-1.mat','mean_Glp','std_Glp','mean_corrcoef_GP','std_corrcoef_GP','mean_vari_Glp','std_vari_Glp')

end

function [t,y]=odeRun(param,totTime,Yini)

tspan = [0 totTime];
[t,y] = ode45(@toggleSwitch,tspan,Yini,[],param);
% plot(t,y,'LineWidth',2)
% legend('gene1','gene2')

end

function [time,Y]=toggleGillespe(totTime,inV,param,OMIGA)

%basic Gillespie alogrithm

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
beta1 = param(5);
beta2 = param(6);

NUM = 10000000;                %total number point
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
function [time,Y]=tauLeapBasal(totTime,inV,param,OMIGA)

%this is a direct tau leap method to simulate the modified 
%toggle switch model(with basal level expression rate

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
beta1 = param(5);
beta2 = param(6);
a0 = param(7);

% NUM = 10000000;                %total number point
Y(1,:) = inV;
tau = 0.05;
NUM = round(totTime/tau);
S = [1,-1,0,0;0,0,1,-1];

for i=1:NUM
    a = [a0*a1*OMIGA+OMIGA*a1*(1-a0)/(1+(Y(i,2)/OMIGA)^n1); beta1*Y(i,1);a0*a2*OMIGA+OMIGA*a2*(1-a0)/(1+(Y(i,1)/OMIGA)^n2); beta2*Y(i,2)];
    new_value = Y(i,:)';
    for j0=1:4
        Num=poissrnd(a(j0)*tau);
        % Make sure things don't go negative
        Use=min([Num new_value(find(S(:,j0)<0))]);
        new_value=new_value+S(:,j0)*Use;
    end
    Y(i+1,:) = new_value';  
end
   time = (0:tau:NUM*tau)';

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
function [time,Y]=GillespeBasalExt(totTime,inV,param,OMIGA)

%this program uses first-reaction Gillespie method to simulate the modified 
%toggle switch model(with basal level expression rate), with extrinsic
%noise

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc = param(6);  %the time scale of extrinsic noise
s = param(7);
sigmae = sqrt(2*s*s/tauc);


NUM = 10000000;                %total number point
Y(1,:) = inV;
% S = [3,-1,0,0;0,0,3,-1];       %introduce birth size, Dec 31,2015
S = [1,-1,0,0;0,0,1,-1];

dt = 0.1;
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);


Next = round(1.1*totTime/dt);
Xtemp = 0;
Xe = zeros(Next,1);
for i1 = 1:Next;
    Xtemp = Xtemp*mu+ Dext*randn;
    Xe(i1) = Xtemp;
end

a1t = a1;   % a1 is the parameter with extrinsic noise
a2t = a2;
b = 1;  %the basical death rate
t = 0;
te = dt;
i = 1;
j = 0;
while(t < totTime && i<= NUM+1)
    a = [OMIGA*a1t*(a0 + (1-a0)/(1+(Y(i,2)/OMIGA)^n1)); b*Y(i,1);OMIGA*a2t*(a0 + (1-a0)/(1+(Y(i,1)/OMIGA)^n2)); b*Y(i,2)];
    [tau,inx] = min(-1./a.*log(rand(4,1)));
    if(t+tau<te)
        Y(i+1,:) = Y(i,:) + S(:,inx)';  
        t = t + tau;  
        time(i+1)= t;
        i = i+1;
    else
%         a1t = a1*(1+Xe(j+1));
%         a2t = a2*(1+Xe(j+1));
        b = 1+Xe(j+1);
        t = te;
        te = te+dt;   %updata the extrinsic time point
        j = j+1;
    end
    
end
time = [0;time(1:end-1)'];
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
function [A,Ae] = jacMatrx(param,vec)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
syms x y E a1 a2 n1 n2 b
J1 = jacobian([a1/(1+y^n1) - x;a2/(1+x^n2) - y],[x,y]);
A = double(subs(J1,{x,y,a1 a2 n1 n2},{vec(1),vec(2),param(1),param(2),param(3),param(4)}));
J2 = jacobian([a1/(1+y^n1) - x + E;a2/(1+x^n2) - y + E;-b*E],[x y E]);
Ae = double(subs(J2,{x,y,E,a1 a2 n1 n2 b},{vec(1),vec(2),0,param(1),param(2),param(3),param(4),param(5)}));

end
function [swiTime,swiRate] = switchingGlp(param,totTime,OMIGA)

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);

tspan = [0 1000];              %integration time for ode
NUM_RUN = 10000;
t_max = 10^6;
swiTime = zeros(7,3);
swiRate = zeros(7,1);

%the B state
a1_list = [10,20,25,30,31,32,33];
for j2 = 1:length(a1_list);
    a1 = a1_list(j2);
    param(1) = a1;
    [t_o,y_o]= ode45(@toggleSwitch,tspan,[40;1],[],param);
    ref = y_o(end,:)*OMIGA;
    region = [0.1*ref(1),5];

    tspan = [0 totTime];
    [t_ode,y_ode] = ode45(@toggleSwitch,tspan,[1;20],[],param);
    Y(1,:) = y_ode(end,:)*OMIGA;

    count = 0;    
    while(count <=200)        
         t =0;
         Y_temp = round(Y(1,:));
         while(t<t_max)
         a = [OMIGA*a1/(1+(Y_temp(2)/OMIGA)^n1); Y_temp(1);...
             OMIGA*a2/(1+(Y_temp(1)/OMIGA).^n2); Y_temp(2)];
         a0 = sum(a);   
         t = t + 1/a0*log(1/rand);
         reference = rand*a0;
         if(reference<a(1))
            Y_temp(1) = Y_temp(1) + 1;               
         elseif(reference > a(1) && (reference <= sum(a(1:2))))
            Y_temp(1) = Y_temp(1) -1;
         elseif(reference > sum(a(1:2)) && (reference <= sum(a(1:3))))
            Y_temp(2) = Y_temp(2) + 1;
         elseif(reference > sum(a(1:3)) && (reference <= sum(a(1:4))))
            Y_temp(2) = Y_temp(2) - 1;     
         end
         if(abs(Y_temp(2)-ref(2))<region(2) && abs(Y_temp(1)-ref(1))<region(1))
            count = count + 1;
            swiTime(j2,count) = t;
            break
         end
         end
    end
    
    
    %compute switching rate in a fixed time
    count_sr = 0;
    for j1=1:NUM_RUN;        
        Y_temp = round(Y(1,:));
        t=0;
        while(t < totTime)
            a = [OMIGA*a1/(1+(Y_temp(2)/OMIGA)^n1); Y_temp(1);...
                OMIGA*a2/(1+(Y_temp(1)/OMIGA).^n2); Y_temp(2)];
            a0 = sum(a);   
            t = t + 1/a0*log(1/rand);
            reference = rand*a0;
            if(reference<a(1))
                Y_temp(1) = Y_temp(1) + 1;
            elseif(reference > a(1) && (reference <= sum(a(1:2))))
                Y_temp(1) = Y_temp(1) -1;
            elseif(reference > sum(a(1:2)) && (reference <= sum(a(1:3))))
                Y_temp(2) = Y_temp(2) + 1;
            elseif(reference > sum(a(1:3)) && (reference <= sum(a(1:4))))
                Y_temp(2) = Y_temp(2) - 1;     
            end
            if(abs(Y_temp(2)-ref(2))<region(2) && abs(Y_temp(1)-ref(1))<region(1))
                count_sr = count_sr + 1;
                break
            end
        end
    end
    swiRate(j2) = count_sr/NUM_RUN;
end
    
    
end
function glpDataStatis()
    %this function load the data of GLP and Glp with extrinsic noise
    a1_list = 2 + (50-2)*(1:1:50)/50;
    
    allMeanSteady1 = [];
    allMeanSteady2 = [];
    allCorr = [];
    allVar1 = [];
    allVar2 = [];
    allEndGlpExt = cell(9,1);
    allEndGlp = cell(9,1);
    dataFolder = '/home/shan/Documents/MATLAB/criticlity/toggleSwitch/allGlp1202';
%     dataFolderEnd = '/home/shan/Documents/MATLAB/criticlity/toggleSwitch/allGlpEnd-1201';
    for i0 = 0:19;
        fileName = fullfile(dataFolder,['GlpCorrVar-1202_',num2str(i0),'.mat']);
        load(fileName);
        allMeanSteady1 = [allMeanSteady1,mean_Glp(:,1)];
        allMeanSteady2 = [allMeanSteady2,mean_Glp(:,2)];
        allCorr = [allCorr,mean_corrcoef_GP];
        allVar1 = [allVar1,mean_vari_Glp(:,1)];
        allVar2 = [allVar2,mean_vari_Glp(:,2)];
        fileEnd = fullfile(dataFolder,['GillespieEndVaule-1202_',num2str(i0),'.mat']);
        load(fileEnd);
        for i2 = 1:9;
            allEndGlp{i2} = [allEndGlp{i2};endGlp{i2}];
            allEndGlpExt{i2} = [allEndGlpExt{i2};endGlpExt{i2}];
        end
    end
    
    %plot the distribution of protein
    figure(1)
    for j1 = 1:9;
        subplot(3,3,j1)
%         [f,xi] = ksdensity(allEndGlp{j1}(:,2));
%         plot(xi,f)
        hist(allEndGlp{j1}(:,2),50)
    end
    figure(2)
    for j1 = 1:9;
        subplot(3,3,j1)
%         [f,xi] = ksdensity(allEndGlpExt{j1}(:,2));
%         plot(xi,f)
        hist(allEndGlpExt{j1}(:,2),50)
    end
    
    %plot the mean correlation
    figure(3)
    errorbar(a1_list',mean(allCorr,2),std(allCorr,0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('correlation coefficient','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    %plot the mean variation
    figure(4)
    subplot(1,2,1)
    plot(a1_list',mean(allVar1,2))
    subplot(1,2,2)
    plot(a1_list',mean(allVar2,2))
end
function [corr,variance] = GlpCorrVar(param,OMIGA,sigmae,beta)
% param is the parameter used in deterministic model
% OMIGA is the system size
% sigmae    the amplitude of extrinsic noise
% beta     the time scale of extrinsic noise

a1 = [5,8,15,20,22,25,30];

%determinisic initial value
tspan = [0 1e4];
Yini = [1,10];
tGlp = 200;   % simulation time of Gillespie methods
for i =3:length(a1);
    param(1) = a1(i);
    [t,y] = ode15s(@toggleSwitch,tspan,Yini,[],param);
    YinGlp = round(OMIGA*y(end,:));
    j1 = 0;
    while(j1 <11);   %each parameter runs 10 times
%         [time,Y]=toggleGillespeBasal(tGlp,YinGlp,param,OMIGA);
        [time,Y]=GillespeBasalExt(tGlp,YinGlp,[param;beta;sigmae],OMIGA);
        plot(time,Y);  %plot and select a interval to calculate corrcoef and variance
        if(i==4 || i == 5) % near the transition point, we need selecet data mannually
            plot(time,Y)
            if(input('the choice:\n'))
                [x,~] = ginput(2);
                if(isempty(x))
                    dataSelect = Y;
                else
                    [~,inx1] = sort(abs(time-x(1)));
                    [~,inx2] = sort(abs(time-x(2)));
%             [left,right] = [inx1(1),inx2(1)];
                    dataSelect = Y(inx1(1):inx2(1),:);   %select data
                end
                j1=j1+1;
                C1 = corrcoef(dataSelect(:,1),dataSelect(:,2));
                corrGlp(i,j1) = C1(1,2);
                varGlp1(i,j1) = var(dataSelect(:,1));
                varGlp2(i,j1) = var(dataSelect(:,2)); 
            end
        else
            dataSelect = Y;
            C1 = corrcoef(dataSelect(:,1),dataSelect(:,2));
            j1 = j1+1;
            corrGlp(i,j1) = C1(1,2);
            varGlp1(i,j1) = var(dataSelect(:,1));
            varGlp2(i,j1) = var(dataSelect(:,2));             
        end
   
    end
            
end
save('corrGlpExt-1209-3.mat','a1','corrGlp','varGlp1','varGlp2')
plot(a1',mean(corrGlp,2),'ro-')

end
function [endGlp,endGlpExt] = GlpDistri(a1_select,param,beta,sigmae,OMIGA,N_dist)
%this function get the distribution of toggle swithch model with and
%without extrinsic noise,using Gillespie method

endGlp = cell(length(a1_select),1);
endGlpExt = cell(length(a1_select),1);
runT = 15;    %the relax time
%the initial value is dertermined by DOE
tspan = [0 1000]; %long enough to make sure it reaches steady state
Ylow = [1;10];
Yhigh = [30;1];
for i1 = 1:length(a1_select)
    param(1) = a1_select(i1);
    [t,y] = ode15s(@toggleSwitch,tspan,Ylow,[],param);
    YinGlpLow = round([OMIGA*y(end,1).*(1+0.2*randn(N_dist/2,1)),OMIGA*y(end,2)*(1+0.2*randn(N_dist/2,1))]);
    [t,y] = ode15s(@toggleSwitch,tspan,Yhigh,[],param);
    YinGlpHigh = round([OMIGA*y(end,1).*(1+0.2*randn(N_dist/2,1)),OMIGA*y(end,2)*(1+0.2*randn(N_dist/2,1))]);;
    allSample = [YinGlpLow;YinGlpHigh];
    
    for i2 = 1:N_dist;
%         iniGlp = round([allSample(i2,1)*20*OMIGA;allSample(i2,2)*10*OMIGA]);
        iniGlp = allSample(i2,:);
        ParamGlpExt = [param;beta;sigmae];
        [t_GP,Y_GP]=toggleGillespeBasal(runT,iniGlp,param,OMIGA);
        endGlp{i1} = [endGlp{i1};Y_GP(round(0.9*size(Y_GP,1)):end,:)];
        [t_GPExt,Y_GPExt]=GillespeBasalExt(runT,iniGlp,ParamGlpExt,OMIGA);
        endGlpExt{i1} = [endGlpExt{i1};Y_GPExt(round(0.9*size(Y_GPExt,1)):end,:)];
    end
end
% save('GillespieEndVaule-12-09.mat','endGlp','endGlpExt','a1_select')

end
function varCorrCLEext(param,a1list,N,runTime,OMIGA)
    %this function simulate the variance, coefficient of correlation and
    %auto correlation
    

fileA = 'toggSwiBifurAll_A.csv';
fileB = 'toggSwiBifurAll_B.csv';
dataA = csvread(fileA);
dataB = csvread(fileB);

dt = 0.1;
for j1 = 1:length(a1list);
    inx_h(j1) = find(dataB(:,1) ==a1list(j1));
end
high = OMIGA*[dataA(inx_h,2),dataB(inx_h,2)];
low = OMIGA*[dataA(inx_h,4),dataB(inx_h,4)];
saddle = OMIGA*[dataA(inx_h,3),dataA(inx_h,3)];
MAXTRY = 5*N;   % maximum trying times


%declare some varibles
meanVal1 = nan(length(high),N);
meanVal2 = nan(length(high),N);
variance1 = nan(length(high),N); 
variance2 = nan(length(high),N);
corrCLE = nan(length(high),N);
lagAuto1 = nan(length(high),N);
lagAuto2 = nan(length(high),N);
dataSelect = cell(length(high),N);

for i0 = 1:length(a1list); 
    param(1) = a1list(i0);
    n = 0;
    flag = 0;
    while(n<N && flag<MAXTRY)
        rng('shuffle')
        [time,Y] = CLEtoggleExtrinsic(runTime,high(i0,:)',dt,param,OMIGA);
%         if(Y(end,1)<saddle(i0,1) && Y(end,2)>saddle && (max(Y(:,2))-min(Y(:,2)))<0.95*(high(i0,2)-low(i0,2)))
         if(param(1)<7 ||(param(1)>=7 && Y(end,1)<saddle(i0,1) && Y(end,2)>saddle(i0,2)))
            n = n+1;
            dataSelect{i0,n} = Y;

%             meanVal1(i0,n) = mean(Y(:,1));
%             meanVal2(i0,n) = mean(Y(:,2));            
%             variance1(i0,n) = var(Y(:,1));
%             variance2(i0,n) = var(Y(:,2));
%             C1 = corrcoef(Y(:,1),Y(:,2));
%             corrCLE(i0,n) = C1(1,2);
%             
%             C3 = corrcoef(Y(1:end-round(1/dt),1),Y(round(1/dt)+1:end,1));
%             lagAuto1(i0,n) = C3(1,2);
%             C4 = corrcoef(Y(1:end-round(1/dt),2),Y(round(1/dt)+1:end,2));
%             lagAuto2(i0,n) = C4(1,2);

         else
            flag = flag + 1;
         end
    end
end
end