%-----------------------------------------------------------------------
% PROGRAM NAME: toggSwiCLEStaticEnsembleExt.m
% DESCIPTION:  
%   this program using CLE to simulate the toggle switch model with
%   extrinsic noise, different time scale and amplitute is considered
%   here, we use an static ensemble approach to mimic the external fluctuation
%   this is suggesed by Chao
% last revised on 11/03/2016
%----------------------------------------------------------------------
function toggSwiCLEStaticEnsembleExt()
close all
clear
clc

a1list = 10:0.5:13;
sigmae = 0.1;
tauc = 100;
% param = [15;10;2;2;0.03;tauc;sigmae];  %when the extrinsic noise terms are same
%param = [15;10;2;2;0.03;1;sigmae;tauc;sigmae];  %when the extrinsic noise terms are different
param = [15;10;2;2;0.05;1;1];  
OMIGA = 50;        %system size
NUM = 1e2;    %repeats times
%runTime = max([20*tauc,1e3]);  %runing time
whichone = 6;
% test initial value 
%OdePara = [15;10;2;2;0.03;1;1];
%y0 = [15;0];
%intiVal = steadyInitial(OdePara,y0);

varCorrCLEextEnsemble(param,a1list,NUM,OMIGA,sigmae,whichone);
end

%%
function varCorrCLEextEnsemble(param,a1list,N,OMIGA,sigmae,whichone)
%this function simulate the variance, coefficient of correlation and
%auto correlation

% param      a vector contain all the basic parameters and amlitude and correlation
% time scale of extrinisc noise
% a1list     the range of control parameter a1
% N          simulation times at each parameter
% runTime    simulation time
% OMIGA      system size, here we use Van Kamppen system size expantion method
% typeExt    type of extrinsic noise, can be added at the degradation rate or other
% parameters
%==================================================================  

%declare some varibles use to store statistical quantities?
steadyVal1 = nan(10*N,length(a1list));
steadyVal2 = nan(10*N,length(a1list));
meanVal = nan(length(a1list),2);
variance = nan(length(a1list),2); 
corrCLE = nan(length(a1list),1);

dt = 0.02;  %step size
runTime = 20/param(6);  %the runining time is about 10 fold of degra
OdePara = param; %this is the reference parameters
y0 = [0,param(2)];

for i0 = 1:length(a1list); %loop over control parameter
    OdePara(1) = a1list(i0);
    
    %ensemble of extrinsic fluctuation
    %allExt = ParaExtEnsem(N,sigmae,param,whichone);
    MaxTry = 10;  %maximum number of try 
    
    for j0 = 1:N
        count = 0;
        flag = 1;
        while(flag && count < MaxTry)
            OdePara(whichone) = ParaExtEnsem(1,sigmae,param,whichone);
            OdePara(7) = OdePara(whichone);  % parameter 6 and 7 are bothe degradation rates
            SteadyVal = FixedPoint(OdePara);
        if (size(SteadyVal,1)>1 || (size(SteadyVal,1)==1 && SteadyVal(2) > SteadyVal(1)))
            flag = 0;
            if(size(SteadyVal,1)==1 && SteadyVal(2) > SteadyVal(1))  % could be two at the saddle point
                HihgIniVal = SteadyVal*OMIGA;
            else
                HihgIniVal = OMIGA*SteadyVal(1,:);
                LowIniVal = OMIGA*SteadyVal(3,:);    
                SaddleIniVal = OMIGA*SteadyVal(2,:);  %saddle point         
            end
        end
            count = count + 1;
        end
            
        [time,Y] = CLEtoggleExtrinsic_diff(runTime,HihgIniVal*(1+0.05*randn),dt,OdePara,OMIGA);
        plot(time,Y)
        % select 10 points in the trajectory before jump
        if(size(SteadyVal,1)==1)
            deltT = round(0.8*runTime/dt/10);
            selectedINX = round(runTime/dt):-deltT:round(0.2*runTime/dt); %index of selected 10 points
        elseif(Y(end,2) < SaddleIniVal(2) &&  Y(end,1) > SaddleIniVal(1)) %jumped to althernative state
            TurningPoint = find(Y(:,2) > SaddleIniVal(2), 1, 'last' );
            selectedINX = randperm(TurningPoint,10);  %randomly select 10 data points before jumping
        elseif(Y(end,2) > SaddleIniVal(2) &&  Y(end,1) < SaddleIniVal(1)) % did not jumped to althernative state 
            deltT = round(0.8*runTime/dt/10);
            selectedINX = round(runTime/dt):-deltT:round(0.2*runTime/dt); %index of selected 10 points
        end
        
        Ys = Y(selectedINX,:);
        steadyVal1((j0-1)*10+1:j0*10,i0) = Ys(1:10,1);
        steadyVal2((j0-1)*10+1:j0*10,i0) = Ys(1:10,2);        
    end
    meanVal(i0,:) = [mean(steadyVal1(:,i0)),mean(steadyVal2(:,i0))];
    variance(i0,:) = [std(steadyVal1(:,i0)),std(steadyVal2(:,i0))];
    tempC = corrcoef(steadyVal1(:,i0),steadyVal2(:,i0));
    corrCLE(i0) = tempC(2,1);
end


 saveFile = ['togg_varCorr_CLE_','tau-',num2str(whichone),'_se',num2str(sigmae),'_N',num2str(OMIGA),'_',date,'.mat'];
 save(saveFile,'a1list','steadyVal1','steadyVal2','meanVal','variance','corrCLE')

%plot the results
%     folder = '/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData';
%     file = fullfile(folder,'simu_varCoefLag_tau1e2_se0.1-0101.mat');
%     load(file);
%     notNanInx = ~isnan(mean(meanVal2,2));
%     figure(1)
% %     plot(dataB(:,1),OMIGA*dataB(:,2))
%     hold on
%     errorbar(a1list(notNanInx)',mean(meanVal2(notNanInx,:),2),std(meanVal2(notNanInx,:),0,2))
%     xlabel('a1','FontSize',24,'FontWeight','Bold')
%     ylabel('mean expression level','FontSize',24,'FontWeight','Bold')
%     set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
%     hold off
    
% %     for i = 1:size(variance2,2)
% %         variance2New(:,i) = variance2(notNanInx,i)./high(notNanInx,2);
% %     end
%     figure(2)
%     hold on
%     errorbar(a1list(notNanInx)',mean(variance2(notNanInx,:),2),std(variance2(notNanInx,:),0,2))
% %     errorbar(a1list(notNanInx)',mean(variance2New(notNanInx,:),2),std(variance2New(notNanInx,:),0,2))
% %     boundedline(a1list(notNanInx)',mean(variance2(notNanInx,:),2),[std(variance2(notNanInx,:),0,2),std(variance2(notNanInx,:),0,2)],'alpha')
%     xlabel('a1','FontSize',24,'FontWeight','Bold')
%     ylabel('variance','FontSize',24,'FontWeight','Bold')
%     set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
%     
%     figure(3)
%     hold on
%     errorbar(a1list(notNanInx)',mean(corrCLE(notNanInx,:),2),std(corrCLE(notNanInx,:),0,2))
% %     boundedline(a1list(notNanInx)',mean(corrCLE(notNanInx,:),2),[std(corrCLE(notNanInx,:),0,2),std(corrCLE(notNanInx,:),0,2)],'alpha')
%     xlabel('a1','FontSize',24,'FontWeight','Bold')
%     ylabel('correlation coefficient of noise','FontSize',24,'FontWeight','Bold')
%     set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
%     
%     figure(4)
%     hold on
%     errorbar(a1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),std(lagAuto2(notNanInx,:),0,2))
% %     boundedline(a1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),[std(lagAuto2(notNanInx,:),0,2),std(lagAuto2(notNanInx,:),0,2)],'alpha')
%     xlabel('a1','FontSize',24,'FontWeight','Bold')
%     ylabel('lag one auto correlation','FontSize',24,'FontWeight','Bold')
%     set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
end

%%
function [time,Y] = CLEtoggleExtrinsic(totTime,inV,dt,param,OMIGA)
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc = param(6);
s = param(7);  %the amplitude of extrinsic noise

sigmae = sqrt(2*s*s/tauc);

N = round(totTime/dt);       %total simulation steps
X = zeros(N,3);
% Xtemp = [inV+0.1*inV.*randn(2,1);normrnd(0,s)];  %revised on Dec 31,2015
Xtemp = [inV;0];         %revised on Dec 31,2015
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);
for i1 = 1:N;
    determin = [a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-(1+Xtemp(3))*Xtemp(1);...
        a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-(1+Xtemp(3))*Xtemp(2)];
    stoch = [sqrt(abs(a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs((1+Xtemp(3))*Xtemp(1))),0,0;...
            0,0,sqrt(abs(a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs((1+Xtemp(3))*Xtemp(2)))];

% determin = [a1*OMIGA*(a0 + (1+Xtemp(3))*(1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-Xtemp(1);...
%         a2*OMIGA*(a0 + (1+Xtemp(3))*(1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-Xtemp(2)];
% stoch = [sqrt(abs(a1*OMIGA*(1+Xtemp(3))*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs(Xtemp(1))),0,0;...
%             0,0,sqrt(abs(a2*OMIGA*(1+Xtemp(3))*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs(Xtemp(2)))];

    Xtemp([1 2]) = Xtemp([1 2]) + determin*dt + stoch*randn(4,1)*sqrt(dt);
    Xtemp(3) = Xtemp(3)*mu+ Dext*randn;
    if(Xtemp(3)<-1)
      Xtemp(3) = -1;
    end
    X(i1,:) = Xtemp;
end

time = (0:dt:N*dt)';
Y = [[inV;0]';X];

end
%%
function [time,Y] = CLEtoggleExtrinsic_diff(totTime,inV,dt,param,OMIGA)
% this function simulate the situation, also considering the situation when 
% independent extrinsic noise are added to each of the gene
% by specifying the input arguments number
% this is actually much like an intrinsic situation

%total number of arguments
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc1 = param(6);
tauc2 = param(7);

N = round(totTime/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = inV';         %revised on Dec 31,2015

for i1 = 1:N;
    determin = [a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-tauc1*Xtemp(1);...
    a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-tauc2*Xtemp(2)];
    stoch = [sqrt(abs(a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs(tauc1*Xtemp(1))),0,0;...
            0,0,sqrt(a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))),-sqrt(abs(tauc2*Xtemp(2)))];        
    Xtemp = Xtemp + determin*dt + stoch*randn(4,1)*sqrt(dt);
    X(i1,:) = Xtemp';
end
Y = [inV;X];
time = (0:dt:N*dt)';

end
%%
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
    if(Xtemp<-1)
        Xtemp=-1;  %avoid unreasonable value of extrinsic noise
    end
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
%     a = [OMIGA*a1t*(a0 + (1-a0)/(1+(Y(i,2)/OMIGA)^n1)); b*Y(i,1);OMIGA*a2t*(a0 + (1-a0)/(1+(Y(i,1)/OMIGA)^n2)); b*Y(i,2)];
    a = [OMIGA*a1*(a0 + a1t*(1-a0)/(1+(Y(i,2)/OMIGA)^n1)); b*Y(i,1);OMIGA*a2*(a0 + a2t*(1-a0)/(1+(Y(i,1)/OMIGA)^n2)); b*Y(i,2)];
    [tau,inx] = min(-1./a.*log(rand(4,1)));
    if(t+tau<te)
        Y(i+1,:) = Y(i,:) + S(:,inx)';  
        t = t + tau;  
        time(i+1)= t;
        i = i+1;
    else
        a1t = a1*(1+Xe(j+1));
        a2t = a2*(1+Xe(j+1));
%         b = 1+Xe(j+1);
        t = te;
        te = te+dt;   %updata the extrinsic time point
        j = j+1;
    end
    
end
time = [0;time(1:end-1)'];
end
%%
function detrendData = dataDetrend(oldData,time,h,N_smooth)
 [ROW,COL] = size(oldData);
 detrendData = zeros(length(time),COL);
 %using kernel smoothing regreesion method to smooth and detrend the data
 for i0 = 1:COL;
    data = ksr(time,oldData(:,i0),h,N_smooth);
    allData = interp1(data.x,data.f,time,'spline');
    plot(time,oldData)
    hold on
    plot(data.x,data.f)
    detrendData(:,i0) = allData- oldData(:,i0);
 end
end

%%
function iniVal = steadyInitial(OdePara,y0)
% this function solve the ode and return the steady steady state value
tspan = [0 1e4];  %long enough
options = odeset('Events',@OdeStop);
[~,y] = ode15s(@toggleSwitch,tspan,y0,options,OdePara);
iniVal = y(end,:);
end
%%
function dydt = toggleSwitch(t,y,param)

% sythestic toggle switch model adopted form James J.Collins' 2000 Nature
% paper, modified version

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc1 = param(6);
tauc2 = param(7);

dydt = zeros(size(y));
u = y(1);
v = y(2);
dudt = a1*(a0 + (1-a0)/(1+v^n1)) - tauc1*u;
dvdt = a2*(a0 + (1-a0)/(1+u^n2)) - tauc2*v;
dydt = [dudt;dvdt];
end

%%
function [value, isterminal,direction] = OdeStop(t,y,param)
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc1 = param(6);
tauc2 = param(7);

dydt = zeros(size(y));
u = y(1);
v = y(2);
dudt = a1*(a0 + (1-a0)/(1+v^n1)) - tauc1*u;
dvdt = a2*(a0 + (1-a0)/(1+u^n2)) - tauc2*v;
dydt = [dudt;dvdt];

value = norm(dydt) - 1.0e-6;
isterminal = 1;
direction = 0;

end

function para = ParaExtEnsem(N,sigmae,OdePara,whichone)
% this function return an ensemble of extrinsic fluctuation
% the fluctuation can be either of the ode parameter
   p0 = OdePara(whichone);
   %more than one parameters need to be set
   if (length(p0) >1 && length(sigmae) > 1)
       para = zeros(N,length(p0));
       for i0 = 1:length(p0)
           para(:,i0) = abs(normrnd(p0(i0),p0(i0)*sigmae(i0),N,1)); 
       end
   else
       para = abs(normrnd(p0,p0*sigmae,N,1));
   end
   
end

function SteadyVal = FixedPoint(param)

% this function find all the fixed points of the equations by finding the
% cross point of nullclines

    a1 = param(1);
    a2 = param(2);
    n1 = param(3);
    n2 = param(4);
    a0 = param(5);
    tauc1 = param(6);
    tauc2 = param(7);
    f1 = @(u,v) a1*(a0 + (1-a0)./(1+v.^n1)) - tauc1*u;
    f2 = @(u,v) a2*(a0 + (1-a0)./(1+u.^n2)) - tauc2*v;
    PlotRange = [0,a1,0,a2];
    g1=ezplot(f1,PlotRange);
    DataPoints1 = get(g1,'ContourMatrix');
    X1 = DataPoints1(1,2:end);
    Y1 = DataPoints1(2,2:end);
    
    g2=ezplot(f2,PlotRange);
    DataPoints2 = get(g2,'ContourMatrix');
    X2 = DataPoints2(1,2:end);
    Y2 = DataPoints2(2,2:end);
    
    %find the cross points of the two nullclines
    [X0,Y0] = intersections(X1,Y1,X2,Y2);
    SteadyVal = [X0,Y0];
end
