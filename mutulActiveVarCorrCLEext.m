%-----------------------------------------------------------------------
% PROGRAM NAME: toggSwiVarCorrCLEext.m
% DESCIPTION:  
%   this program using CLE to simulate mutual activation model with
%   extrinsic noise, different time scale and amplitute is considered
% last revised on Dec.2,2015
%----------------------------------------------------------------------
function mutulActiveVarCorrCLEext()
close all
clear
clc

k1list = 0.1:0.05:0.55;
sigmae = 0.05;
tauc = 100;
param = [0.3;0.5;3;3;0.05;tauc;sigmae];
OMIGA =500;        %system size
NUM = 100;
runTime = max([2*tauc,1e3]);

varCorrCLEext(param,k1list,NUM,runTime,OMIGA);
end

function varCorrCLEext(param,k1list,N,runTime,OMIGA)
    %this function simulate the variance, coefficient of correlation and
    %auto correlation
    

fileA = 'mutualActiveBifurA.csv';
fileB = 'mutualActiveBifurB.csv';
dataA = csvread(fileA);
dataB = csvread(fileB);

dt = 0.05;
for j1 = 1:length(k1list);
    inx_h(j1) = find(dataB(:,1) ==round(k1list(j1)*100)/100);
end
high = OMIGA*[dataA(inx_h,4),dataB(inx_h,4)];
% low = OMIGA*[dataA(inx_h,4),dataB(inx_h,4)];
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
% dataSelect = cell(length(high),N);

for i0 = 7:length(k1list); 
    param(1) = k1list(i0);
    n = 0;
    flag = 0;
    while(n<N && flag<MAXTRY)
        rng('shuffle')
        [time,Y] = CLEtoggleExtrinsic(runTime,high(i0,:)',dt,param,OMIGA);
%         [time,Y] = GillespeBasalExt(runTime,round(high(i0,:)'),param,OMIGA);
         if(param(1)<0.15 ||(param(1)>=0.15 && Y(end,1)>saddle(i0,1) && Y(end,2)>saddle(i0,2)))
            n = n+1;
%             dataSelect{i0,n} = Y;
            meanVal1(i0,n) = mean(Y(:,1));
            meanVal2(i0,n) = mean(Y(:,2));            
            variance1(i0,n) = var(Y(:,1));
            variance2(i0,n) = var(Y(:,2));
            C1 = corrcoef(Y(:,1),Y(:,2));
            corrCLE(i0,n) = C1(1,2);
            
            C3 = corrcoef(Y(1:end-round(1/dt),1),Y(round(1/dt)+1:end,1));
            lagAuto1(i0,n) = C3(1,2);
            C4 = corrcoef(Y(1:end-round(1/dt),2),Y(round(1/dt)+1:end,2));
            lagAuto2(i0,n) = C4(1,2);

         else
            flag = flag + 1;
         end
    end
end


saveFile = ['mutuActive_varCorrLag_CLE_','tau-',num2str(param(6)),'_se',num2str(param(7)),'_N',num2str(OMIGA),'_',date,'.mat'];
save(saveFile,'a1list','meanVal1','meanVal2','variance1','variance2','corrCLE','lagAuto1','lagAuto2')

%plot the results
%     folder = '/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData';
%     file = fullfile(folder,'simu_varCoefLag_tau1e2_se0.1-0101.mat');
%     load(file);
    notNanInx = ~isnan(mean(meanVal2,2));
    figure(1)
%     plot(dataB(:,1),OMIGA*dataB(:,2))
    hold on
    errorbar(k1list(notNanInx)',mean(meanVal2(notNanInx,:),2),std(meanVal2(notNanInx,:),0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('mean expression level','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    hold off
    
%     variance2CV = zeros(sum(notNanInx),size(variance2,2));
%     meanVAl2 = mean(meanVal2(notNanInx,:),2);
%     for i = 1:size(variance2,2);
%         variance2CV(:,i) = variance2(notNanInx,i)./meanVal2(notNanInx,i);
%     end
    figure(2)
    hold on
    errorbar(k1list(notNanInx)',mean(variance2(notNanInx,:),2),std(variance2(notNanInx,:),0,2))
%         errorbar(k1list(notNanInx)',mean(variance2CV(notNanInx,:),2),std(variance2CV(notNanInx,:),0,2))

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
    determin = [OMIGA*(a0 + (1-a0)*(Xtemp(2)/OMIGA)^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1))-(1+Xtemp(3))*Xtemp(1);...
        OMIGA*(a0 + (1-a0)*(Xtemp(1)/OMIGA)^n2/(k2^n2+(Xtemp(1)/OMIGA)^n2))-(1+Xtemp(3))*Xtemp(2)];
    stoch = [sqrt(abs(OMIGA*(a0 + (1-a0)*(Xtemp(2)/OMIGA)^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs((1+Xtemp(3))*Xtemp(1))),0,0;...
            0,0,sqrt(abs(OMIGA*(a0 + (1-a0)*(Xtemp(1)/OMIGA)^n2/(k2^n2+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs((1+Xtemp(3))*Xtemp(2)))];

%     determin = [OMIGA*(a0 + (1+Xtemp(3))*(1-a0)*(Xtemp(2)/OMIGA)^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1))-Xtemp(1);...
%         OMIGA*(a0 + (1+Xtemp(3))*(1-a0)*(Xtemp(1)/OMIGA)^n2/(k2^n2+(Xtemp(1)/OMIGA)^n2))-Xtemp(2)];
%     stoch = [sqrt(abs(OMIGA*(a0 + (1+Xtemp(3))*(1-a0)*(Xtemp(2)/OMIGA)^n1/(k1^n1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs(Xtemp(1))),0,0;...
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
    a = [OMIGA*(a0 + (1-a0)*(Y(i,2)/OMIGA)^n1/(k1^n1+(Y(i,2)/OMIGA)^n1)); b*Y(i,1);OMIGA*(a0 + (1-a0)*(Y(i,1)/OMIGA)^n2/(k2^n2+(Y(i,1)/OMIGA)^n2)); b*Y(i,2)];
%     a = [OMIGA*(a0 + a1t*(1-a0)/(1+(Y(i,2)/OMIGA)^n1)); b*Y(i,1);OMIGA*a2*(a0 + a2t*(1-a0)/(1+(Y(i,1)/OMIGA)^n2)); b*Y(i,2)];
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