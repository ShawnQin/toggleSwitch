%-----------------------------------------------------------------------
% PROGRAM NAME: toggSwiVarCorrGLPext.m
% DESCIPTION:  
%   this program using GLP algorithm to simulate the toggle switch model with
%   extrinsic noise, and compared with CLE approximation
%   V2: a new modified Gillespie algorithm named "Extrand" is implementated
%   which gives much faster speed
%   last revised on 11/09/2016
%----------------------------------------------------------------------
function toggSwiVarCorrGLPext()
close all
clear
clc

a1list = 4:0.5:15;
sigmae = 0.1;
tauc = 1;
param = [15;10;2;2;0.03;tauc;sigmae];
OMIGA = 20;        %system size
NUM = 10;
runTime = max([10*tauc,1e2]);


varCorrGLPext(param,a1list,NUM,runTime,OMIGA);
end

function varCorrGLPext(param,a1list,N,runTime,OMIGA)
    %this function simulate the variance, coefficient of correlation and
    %auto correlation
    

% fileA = 'toggSwiBifurAll_A.csv';
% fileB = 'toggSwiBifurAll_B.csv';
fileA = '/lustre1/tangc_pkuhpc/ssqin/toggleSwi/toggSwiBifurAll_A.csv';
fileB = '/lustre1/tangc_pkuhpc/ssqin/toggleSwi/toggSwiBifurAll_B.csv';
dataA = csvread(fileA);
dataB = csvread(fileB);

dt = 0.01;
for j1 = 1:length(a1list)
    inx_h(j1) = find(dataB(:,1) ==a1list(j1));
end
high = round(OMIGA*[dataA(inx_h,2),dataB(inx_h,2)]);
saddle = round(OMIGA*[dataA(inx_h,3),dataA(inx_h,3)]);
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

for i0 = 1:length(a1list) 
    param(1) = a1list(i0);
    n = 0;
    flag = 0;
    while(n<N && flag<MAXTRY)
        rng('shuffle')
        [time,Y] = GillespeBasalExt(runTime,round(high(i0,:)'),param,OMIGA,dt);        
%         [time, Y] = ExtrandeExt(runTime,round(high(i0,:)'),param,OMIGA,dt);

         if(param(1)<7 ||(param(1)>=7 && Y(end,1)<saddle(i0,1) && Y(end,2)>saddle(i0,2)))
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

saveFile = ['togg_varCorrLag_GLP_','tau-',num2str(param(6)),'_se',num2str(param(7)),'_N',num2str(OMIGA),'.mat'];
save(saveFile,'a1list','meanVal1','meanVal2','variance1','variance2','corrCLE','lagAuto1','lagAuto2')

% %plot the results
% %     folder = '/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData';
% %     file = fullfile(folder,'simu_varCoefLag_tau1e2_se0.1-0101.mat');
% %     load(file);
%     notNanInx = ~isnan(mean(meanVal2,2));
%     figure(1)
% %     plot(dataB(:,1),OMIGA*dataB(:,2))
%     hold on
%     errorbar(a1list(notNanInx)',mean(meanVal2(notNanInx,:),2),std(meanVal2(notNanInx,:),0,2))
%     xlabel('a1','FontSize',24,'FontWeight','Bold')
%     ylabel('mean expression level','FontSize',24,'FontWeight','Bold')
%     set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
%     hold off
%     
% %     for i = 1:size(variance2,2)
% %         variance2New(:,i) = variance2(notNanInx,i)./mean(meanVal2(notNanInx,:),2);
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
function [timeY,Y]=GillespeBasalExt(totTime,inV,param,OMIGA,dt)

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


Y(1,:) = inV;
% S = [3,-1,0,0;0,0,3,-1];       %introduce birth size, Dec 31,2015
S = [1,-1,0,0;0,0,1,-1];

mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);
Next = round(1.1*totTime/dt);
Xtemp = 0;
Xe = zeros(Next,1);
for i1 = 1:Next
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
Ytemp = Y(1,:);
timeY = 0;
while(t < totTime)
    a = [OMIGA*a1*(a0 + (1-a0)/(1+(Ytemp(2)/OMIGA)^n1)); b*Ytemp(1);OMIGA*a2*(a0 + (1-a0)/(1+(Ytemp(1)/OMIGA)^n2)); b*Ytemp(2)];
    [tau,inx] = min(-1./a.*log(rand(4,1)));
    if(t+tau<te)
        Ytemp = Ytemp + S(:,inx)';  
        t = t + tau;
%         Y = [Y;Ytemp];
%         timeY = [timeY;t];
%         time(i+1)= t;
    else
%         a1t = a1*(1+Xe(j+1));
%         a2t = a2*(1+Xe(j+1));
        Y = [Y;Ytemp];
        timeY = [timeY;t];
        
        if (Xe(j+1) > -1)
            b = 1+Xe(j+1);
        else
            b = 0;
        end
        j = j+1;
        time(j) = t;
%         Y(j,:) = Ytemp;  %only record part of the data
        t = te;
        te = te+dt;   %updata the extrinsic time point
    end
    
end
time = [0;time(1:end-1)'];
end

function [time, Y] = ExtrandeExt(totTime,inV,param,OMIGA,dt)
% this is an implecation of Extrande algoithm, a modified Gillespie algorithm
% that can simulate time-dependent propensity in a much effective way than exiting
% methods

% the algorithm rely on a virtual reaction, in our case we use OU process to simulate
% time-dependent degradation rate

%parameters
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc = param(6);  %the autocorrelation time scale of extrinsic noise
s = param(7);
sigmae = sqrt(2*s*s/tauc);


Y(1,:) = inV;
% S = [3,-1,0,0;0,0,3,-1];       %introduce birth size, Dec 31,2015
S = [1,-1,0,0;0,0,1,-1];  %chemical stoichemetrics

mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);

Next = round(1.1*totTime/dt);
Xtemp = 0;
b = 1;  %the basical degradation rate


% simulation of extrinsic noise
Xe = zeros(Next,1);
for i1 = 1:Next;
    Xtemp = Xtemp*mu+ Dext*randn;
    if(Xtemp<-1)
        Xtemp=-1;  %avoid unreasonable value of extrinsic noise
    end
    Xe(i1) = Xtemp;
end

t = 0;
Ytemp = Y(1,:);
time = 0;
%implementation of Extrande algorithm
LookTime = min(100*dt,totTime - t);  % look achead time
while (t < totTime)
    inx1 = max(round(t/dt),1); % get rid of 0 index
    inx2 = round((t+LookTime)/dt);
    XeVec = (1+Xe(inx1:inx2))*b;
    PropBound = PropensityBound(param,OMIGA,Ytemp,XeVec);
    tau = - log(rand)/PropBound;
    if(tau > LookTime)
        t = t + LookTime;
    else
        t = t + tau;
        if(floor(t/dt) == 0)
            currXe = 1+ Xe(1)/dt*t;
        else
            currXe = 1+(Xe(ceil(t/dt))-Xe(floor(t/dt)))/dt*(t-ceil(t/dt)*dt);
        end
        
        a = [OMIGA*a1*(a0 + (1-a0)/(1+(Ytemp(2)/OMIGA)^n1));currXe*Ytemp(1);...
            OMIGA*a2*(a0 + (1-a0)/(1+(Ytemp(1)/OMIGA)^n2)); currXe*Ytemp(2)];
        threshold = PropBound*rand;
        if(sum(a)>= threshold)
            whihchReaction = find(cumsum(a)> threshold,1); %which reaction will fire
            Ytemp = Ytemp + S(:,whihchReaction)';  %update number
            Y = [Y;Ytemp];
            time = [time;t];
        end
    end
      
end


end

function propbound = PropensityBound(param,OMIGA,CurrState,XeVec)
% this fucniton return a propensity bound in the Extrande algorithm
% since propensity is a monotonic function of extrinsic noise, in this case it is 
% the fluctuated degradation rate
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);

MaxExt = max(XeVec);
a = [OMIGA*a1*(a0 + (1-a0)/(1+(CurrState(2)/OMIGA)^n1));MaxExt*CurrState(1);...
    OMIGA*a2*(a0 + (1-a0)/(1+(CurrState(1)/OMIGA)^n2)); MaxExt*CurrState(2)];
propbound = sum(a);

end
