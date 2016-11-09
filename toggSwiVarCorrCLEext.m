%-----------------------------------------------------------------------
% PROGRAM NAME: toggSwiVarCorrCLEext.m
% DESCIPTION:  
%   this program using CLE to simulate the toggle switch model with
%   extrinsic noise, different time scale and amplitute is considered
%   can also specify if use detrending method befor calculating EWS

% last revised on 11/05/2016
%----------------------------------------------------------------------
function toggSwiVarCorrCLEext()
close all
clear
clc

a1list = 4:0.5:15;
sigmae = 0.1;
tauc = 1000;
% param = [15;10;2;2;0.03;tauc;sigmae];  %when the extrinsic noise terms are same
param = [15;10;2;2;0.03;1;sigmae;tauc;sigmae];  %when the extrinsic noise terms are different
OMIGA = 50;        %system size
NUM = 10;    %repeats times
runTime = max([20*tauc,1e3]);  %runing time
detrendOrNot = 0;  %0 for original, 1 for detrending

varCorrCLEext(param,a1list,NUM,runTime,OMIGA,detrendOrNot);
end

function varCorrCLEext(param,a1list,N,runTime,OMIGA,detrendOrNot)
    %this function simulate the variance, coefficient of correlation and
    %auto correlation
    
%direct load the steady state of corresponding ode system
fileA = 'toggSwiBifurAll_A.csv';
fileB = 'toggSwiBifurAll_B.csv';
dataA = csvread(fileA);
dataB = csvread(fileB);

dt = 0.05;
for j1 = 1:length(a1list);
    inx_h(j1) = find(dataB(:,1) ==a1list(j1));
end
high = OMIGA*[dataA(inx_h,2),dataB(inx_h,2)];  %state of high
% low = OMIGA*[dataA(inx_h,4),dataB(inx_h,4)];
saddle = OMIGA*[dataA(inx_h,3),dataB(inx_h,3)]; %saddle point
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

for i0 = 15:length(a1list); 
    param(1) = a1list(i0);
    n = 0;
    flag = 0;
    while(n<=N && flag<MAXTRY)
        rng('shuffle')
        [time,Y] = CLEtoggleExtrinsic_diff(runTime,high(i0,:)',dt,param,OMIGA);
%         [time,Y] = GillespeBasalExt(runTime,round(high(i0,:)'),param,OMIGA);
         if(param(1)<7 ||(param(1)>=7 && Y(end,1)<saddle(i0,1) && Y(end,2)>saddle(i0,2)))
            n = n+1;
%             dataSelect{i0,n} = Y;
            h = 20;
            N_smooth = round(size(Y,1)/50);
            if(detrendOrNot)  %using detrending method
                Y_detr = [dataDetrend(Y(:,1),time,h,N_smooth),dataDetrend(Y(:,2),time,h,N_smooth)];
            else
                Y_detr = Y;
            end
%detrending Jan 19th,2016
            meanVal1(i0,n) = mean(Y_detr(:,1));
            meanVal2(i0,n) = mean(Y_detr(:,2));            
            variance1(i0,n) = var(Y_detr(:,1));
            variance2(i0,n) = var(Y_detr(:,2));
            C1 = corrcoef(Y_detr(:,1),Y_detr(:,2));
            corrCLE(i0,n) = C1(1,2);
            
            C3 = corrcoef(Y_detr(1:end-round(1/dt),1),Y_detr(round(1/dt)+1:end,1));
            lagAuto1(i0,n) = C3(1,2);
            C4 = corrcoef(Y_detr(1:end-round(1/dt),2),Y_detr(round(1/dt)+1:end,2));
            lagAuto2(i0,n) = C4(1,2);

         else
            flag = flag + 1;
         end
    end
end


% saveFile = ['togg_varCorrLag_CLE_','tau-',num2str(param(6)),'_se',num2str(param(7)),'_N',num2str(OMIGA),'_',date,'.mat'];
% save(saveFile,'a1list','meanVal1','meanVal2','variance1','variance2','corrCLE','lagAuto1','lagAuto2')

%plot the results
%     folder = '/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData';
%     file = fullfile(folder,'simu_varCoefLag_tau1e2_se0.1-0101.mat');
%     load(file);
    notNanInx = ~isnan(mean(meanVal2,2));
    figure(1)
%     plot(dataB(:,1),OMIGA*dataB(:,2))
    hold on
    errorbar(a1list(notNanInx)',mean(meanVal2(notNanInx,:),2),std(meanVal2(notNanInx,:),0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('mean expression level','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    hold off
    
%     for i = 1:size(variance2,2)
%         variance2New(:,i) = variance2(notNanInx,i)./high(notNanInx,2);
%     end
    figure(2)
    hold on
    errorbar(a1list(notNanInx)',mean(variance2(notNanInx,:),2),std(variance2(notNanInx,:),0,2))
%     errorbar(a1list(notNanInx)',mean(variance2New(notNanInx,:),2),std(variance2New(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(variance2(notNanInx,:),2),[std(variance2(notNanInx,:),0,2),std(variance2(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('variance','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    
    figure(3)
    hold on
    errorbar(a1list(notNanInx)',mean(corrCLE(notNanInx,:),2),std(corrCLE(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(corrCLE(notNanInx,:),2),[std(corrCLE(notNanInx,:),0,2),std(corrCLE(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('correlation coefficient of noise','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    
    figure(4)
    hold on
    errorbar(a1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),std(lagAuto2(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),[std(lagAuto2(notNanInx,:),0,2),std(lagAuto2(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('lag one auto correlation','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
end

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
function [time,Y] = CLEtoggleExtrinsic_diff(totTime,inV,dt,param,OMIGA)
%this function simulate the situation when independent extrinsic noise are
%added to each of the gene

a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc1 = param(6);
s1 = param(7);  %the amplitude of extrinsic noise
tauc2 = param(8);
s2 = param(9);  %the amplitude of extrinsic noise

sigmae1 = sqrt(2*s1*s1/tauc1);
sigmae2 = sqrt(2*s2*s2/tauc2);

N = round(totTime/dt);       %total simulation steps
X = zeros(N,4);
% Xtemp = [inV+0.1*inV.*randn(2,1);normrnd(0,s)];  %revised on Dec 31,2015
Xtemp = [inV;0;0];         %revised on Dec 31,2015
mu1 = exp(-dt/tauc1);
mu2 = exp(-dt/tauc2);
Dext1 = sigmae1*sqrt((1-mu1^2)*tauc1/2);
Dext2 = sigmae2*sqrt((1-mu2^2)*tauc2/2);
for i1 = 1:N;
%     determin = [a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-(1+Xtemp(3))*Xtemp(1);...
%         a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-(1+Xtemp(4))*Xtemp(2)];
%     stoch = [sqrt(abs(a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs((1+Xtemp(3))*Xtemp(1))),0,0;...
%             0,0,sqrt(abs(a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs((1+Xtemp(4))*Xtemp(2)))];

    determin = [a1*OMIGA*(a0 + (1+Xtemp(3))*(1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-Xtemp(1);...
        a2*OMIGA*(a0 + (1+Xtemp(4))*(1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-Xtemp(2)];
    stoch = [sqrt(abs(a1*OMIGA*(a0 + (1+Xtemp(3))*(1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs(Xtemp(1))),0,0;...
            0,0,sqrt(abs(a2*OMIGA*(a0 + (1+Xtemp(4))*(1-a0)/(1+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs(Xtemp(2)))];
        
    Xtemp([1 2]) = Xtemp([1 2]) + determin*dt + stoch*randn(4,1)*sqrt(dt);
    Xtemp(3) = Xtemp(3)*mu1+ Dext1*randn;
    Xtemp(4) = Xtemp(4)*mu2+ Dext2*randn;
    if(Xtemp(3)<-1)
      Xtemp(3) = -1;
    end
    if(Xtemp(4)<-1)
      Xtemp(4) = -1;
    end
    X(i1,:) = Xtemp;
end

time = (0:dt:N*dt)';
Y = [[inV;0;0]';X];

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
