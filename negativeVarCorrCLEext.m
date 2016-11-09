%-----------------------------------------------------------------------
% PROGRAM NAME: negativeVarCorrCLEext.m
% DESCIPTION:  
%   this program using CLE to simulate mutual activation model with
%   extrinsic noise, different time scale and amplitute is considered
% last revised on Dec.2,2015
%----------------------------------------------------------------------
function negativeVarCorrCLEext()
% close all
clear
clc

k1list = 0.1:0.2:2;
sigmae = 0.1;
tauc = 100;
param = [0.3;0.5;3;3;0.05;tauc;sigmae];
OMIGA =500;        %system size
NUM = 100;
runTime = max([20*tauc,5e3]);
% runTime = 500;

varCorrCLEext(param,k1list,NUM,runTime,OMIGA);
end

function varCorrCLEext(param,k1list,N,runTime,OMIGA)

iniODE = [0.5;0.5];  %since this system has only one stable state, the initial value doesn't matter
tspan = [0 1e3];
dt = 0.05;

%declare some varibles
meanVal1 = nan(length(k1list),N);
meanVal2 = nan(length(k1list),N);
variance1 = nan(length(k1list),N); 
variance2 = nan(length(k1list),N);
corrCLE = nan(length(k1list),N);
lagAuto1 = nan(length(k1list),N);
lagAuto2 = nan(length(k1list),N);

for i0 = 1:length(k1list); 
    param(1) = k1list(i0);
    [t,y] = ode15s(@negativeODE,tspan,iniODE,[],param);
    iniVal(i0,:) = OMIGA*y(end,:);

   for j0 = 1:N;
        rng('shuffle')
        [time,Yraw] = CLEtoggleExtrinsic(runTime,iniVal(i0,:)',dt,param,OMIGA);
%         [time,Y] = GillespeBasalExt(runTime,round(iniVal(i0,:)'),param,OMIGA);
%             Y = Yraw; %no detrending
            h = 20;  %detrending window
            num_smooth = round(size(Yraw,1)/50);  %number of point used to smooth
            Y = [dataDetrend(Yraw(:,1),time,h,num_smooth),dataDetrend(Yraw(:,2),time,h,num_smooth)];  %use the detrended data,the average window h=30 can be tuned
            meanVal1(i0,j0) = mean(Y(:,1));
            meanVal2(i0,j0) = mean(Y(:,2));             
            variance1(i0,j0) = var(Y(:,1))/mean(Yraw(:,1));  % actuall CV, modified 01/20/2016
            variance2(i0,j0) = var(Y(:,2))/mean(Yraw(:,2));  % actuall CV, modified 01/20/2016
            C1 = corrcoef(Y(:,1),Y(:,2));
            corrCLE(i0,j0) = C1(1,2);
            
            C3 = corrcoef(Y(1:end-round(1/dt),1),Y(round(1/dt)+1:end,1));
            lagAuto1(i0,j0) = C3(1,2);
            C4 = corrcoef(Y(1:end-round(1/dt),2),Y(round(1/dt)+1:end,2));
            lagAuto2(i0,j0) = C4(1,2);
    end
end


saveFile = ['negative_varCorrLag_CLE_','tau-',num2str(param(6)),'_se',num2str(param(7)),'_N',num2str(OMIGA),'_',date,'.mat'];
save(saveFile,'k1list','meanVal1','meanVal2','variance1','variance2','corrCLE','lagAuto1','lagAuto2')

    notNanInx = ~isnan(mean(meanVal2,2));
    

    figure(1)
%     plot(dataB(:,1),OMIGA*dataB(:,2))
    hold on
    errorbar(k1list(notNanInx)',mean(meanVal2(notNanInx,:),2),std(meanVal2(notNanInx,:),0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('mean expression level','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    hold off
    
    
%     for i = 1:size(variance2,2)
%         variance2New(:,i) = variance2(notNanInx,i)./mean(meanVal2(notNanInx,:),2);
%     end
    figure(2)
    hold on
    errorbar(k1list(notNanInx)',mean(variance2(notNanInx,:),2),std(variance2(notNanInx,:),0,2))
%     errorbar(k1list(notNanInx)',mean(variance2New(notNanInx,:),2),std(variance2New(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(variance2(notNanInx,:),2),[std(variance2(notNanInx,:),0,2),std(variance2(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('variance','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    
    figure(3)
    hold on
    errorbar(k1list(notNanInx)',mean(corrCLE(notNanInx,:),2),std(corrCLE(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(corrCLE(notNanInx,:),2),[std(corrCLE(notNanInx,:),0,2),std(corrCLE(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('correlation coefficient of noise','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    
    figure(4)
    hold on
    errorbar(k1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),std(lagAuto2(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),[std(lagAuto2(notNanInx,:),0,2),std(lagAuto2(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('lag one auto correlation','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
end

function [time,Y] = CLEtoggleExtrinsic(totTime,inV,dt,param,OMIGA)
k1 = param(1);
k2 = param(2);
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
    determin = [OMIGA*(a0 + (1-a0)*k1^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1))-(1+Xtemp(3))*Xtemp(1);...
        OMIGA*(a0 + (1-a0)*(Xtemp(1)/OMIGA)^n2/(k2^n2+(Xtemp(1)/OMIGA)^n2))-(1+Xtemp(3))*Xtemp(2)];
    stoch = [sqrt(abs(OMIGA*(a0 + (1-a0)*k1^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs((1+Xtemp(3))*Xtemp(1))),0,0;...
            0,0,sqrt(abs(OMIGA*(a0 + (1-a0)*(Xtemp(1)/OMIGA)^n2/(k2^n2+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs((1+Xtemp(3))*Xtemp(2)))];

%     determin = [OMIGA*(a0 + (1+Xtemp(3))*(1-a0)*k1^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1))-Xtemp(1);...
%         OMIGA*(a0 + (1+Xtemp(3))*(1-a0)*(Xtemp(1)/OMIGA)^n2/(k2^n2+(Xtemp(1)/OMIGA)^n2))-Xtemp(2)];
%     stoch = [sqrt(abs(OMIGA*(a0 + (1+Xtemp(3))*(1-a0)*k1^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs(Xtemp(1))),0,0;...
%             0,0,sqrt(abs(OMIGA*(a0 + (1+Xtemp(3))*(1-a0)*(Xtemp(1)/OMIGA)^n2/(k2^n2+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs(Xtemp(2)))];

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
function [time,Y]=GillespeBasalExt(totTime,inV,param,OMIGA)

%this program uses first-reaction Gillespie method to simulate the modified 
%toggle switch model(with basal level expression rate), with extrinsic
%noise

k1 = param(1);
k2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc = param(6);
s = param(7);  %the amplitude of extrinsic noise

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

a1t = 1;   % a1 is the parameter with extrinsic noise
a2t = 1;
b = 1;  %the basical death rate
t = 0;
te = dt;
i = 1;
j = 0;
while(t < totTime && i<= NUM+1)
    a = [OMIGA*(a0 + (1-a0)*k1^n1/(k1^n1+(Y(i,2)/OMIGA)^n1)); b*Y(i,1);OMIGA*(a0 + (1-a0)*(Y(i,1)/OMIGA)^n2/(k2^n2+(Y(i,1)/OMIGA)^n2)); b*Y(i,2)];
%     a = [OMIGA*(a0 + a1t*(1-a0)*k1^n1/(1+(Y(i,2)/OMIGA)^n1)); b*Y(i,1);OMIGA*(a0 + a2t*(1-a0)/(1+(Y(i,1)/OMIGA)^n2)); b*Y(i,2)];
    [tau,inx] = min(-1./a.*log(rand(4,1)));
    if(t+tau<te)
        Y(i+1,:) = Y(i,:) + S(:,inx)';  
        t = t + tau;  
        time(i+1)= t;
        i = i+1;
    else
%         a1t = a1*(1+Xe(j+1));
%         a2t = a2*(1+Xe(j+1));
    if(1+Xe(j+1)>-1)
        b = 1+Xe(j+1);
    else
        b = -1;
    end
        t = te;
        te = te+dt;   %updata the extrinsic time point
        j = j+1;
    end
    
end
time = [0;time(1:end-1)'];
end

function dydt = negativeODE(t,y,param)
    k1 = param(1);
    k2 = param(2);
    n1 = param(3);
    n2 = param(4);
    a0 = param(5);
    
    dydt = [a0 + (1-a0)*k1^n1/(k1^n1 + y(2)^n2) - y(1); a0 + (1-a0)*y(1)^n2/(k2^n2 + y(1)^n2) - y(2)];
end
function detrendData = dataDetrend(oldData,time,h,num_smooth)
 [ROW,COL] = size(oldData);
 detrendData = zeros(length(time),COL);
 %using kernel smoothing regreesion method to smooth and detrend the data
 for i0 = 1:COL;
    data = ksr(time,oldData(:,i0),h,num_smooth);
    allData = interp1(data.x,data.f,time,'spline');
%     plot(time,oldData)
%     hold on
%     plot(data.x,data.f)
    detrendData(:,i0) = allData- oldData(:,i0);
 end
end
