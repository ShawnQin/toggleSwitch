%************************************************************************
%PROGRAM NAME: toggTimeDetrend()
%DESCRIPTION:
%    this program show how extrinsic noise make the signature of
%    ciriticality even more doom. Considering the changing paramers
%CHANGE: using Luolan Chen's moment expasion method
%LAST REVISED: June 9,2016
%*************************************************************************
function toggTimeDetrendExtri()
close all
clear
clc

as= 4;   %start control parameter
ae = 20; %end parameter
time_change = 1e3;   %the simulation time, during which the control paramerter changes from as to ae linearly
sigmae = 0;  %std of extrinsic noise
tauc = 10;   %correlation time scale of extrinsic noise
param = [as;10;2;2;0.03;tauc;sigmae];  %parameters
OMIGA = 50;             %system size
NUM = 1;
paraChgCLEext(param,as,ae,NUM,time_change,OMIGA);

% testChenModel(1,time_change);
% JumpTimeChenModel(1e2);

% ChenModelHist()

end

function paraChgCLEext(param,as,ae,N,time_change,OMIGA)


    %if unspecified, use my model
    fileA = 'toggSwiBifurAll_A.csv';
    fileB = 'toggSwiBifurAll_B.csv';
    dataA = csvread(fileA);
    dataB = csvread(fileB);
    inx_iniVal = find(dataB(:,1) == as);
    iniVal = OMIGA*[dataA(inx_iniVal,2),dataB(inx_iniVal,2)];
    dt = 0.01;
    
    % stable states and saddle nodes
%     for j1 = 1:length(a1list);
%         inx_h(j1) = find(dataB(:,1) ==a1list(j1));
%     end
    inx_h = 114;
    high = OMIGA*[dataA(inx_h,2),dataB(inx_h,2)];
    low = OMIGA*[dataA(inx_h,4),dataB(inx_h,4)];
    saddle = OMIGA*[dataA(inx_h,3),dataB(inx_h,3)];
    
%declare some varibles
variance1 = cell(N,1);
variance2 = cell(N,1);
coeffVar1 = cell(N,1);
coeffVar2 = cell(N,1);
FanoFactor1 = cell(N,1);
FanoFactor2 = cell(N,1);
corrCLE = cell(N,1);
lagAuto1 = cell(N,1);
lagAuto2 = cell(N,1);
allTimeSelect = cell(N,1);
alla1 = cell(N,1);
timeJump = nan(N,1);
a1tJump = nan(N,1);
ewsDetr = struct('a1Jump',[],'tJump',[],'var1',[],'var2',[],...
           'cv1',[],'cv2',[],'corr',[],'lag1',[],'lag2',[],'kv1',[],...
           'kv2',[],'kcv1',[],'kcv2',[],'kcurr',[],...
           'klag1',[],'klag2',[],'Fano1',[],'Fano2',[],'kf1',[],'kf2',[],...
           'time',[],'a1chg',[]);
ewsDetr. kv1 = nan(N,1);
ewsDetr. kv2 = nan(N,1);
ewsDetr. kcv1 = nan(N,1);
ewsDetr. kcv2 = nan(N,1);
ewsDetr. kcurr = nan(N,1);
ewsDetr. klag1 = nan(N,1);
ewsDetr. klag2 = nan(N,1);
ewsDetr. kf1 = nan(N,1);
ewsDetr. kf2 = nan(N,1);
% variance_slide_moment = cell(N,1);
% corrCLE_moment = cell(N,1);
% lagAuto_slide_moment = cell(N,1);



for i = 1:N;
    [time,Y,a1t] = CLEtoggleExtrChgPara(param,as,ae,time_change,iniVal',dt,OMIGA);
    inxJump = find(Y(:,1) > 1.3*saddle(1) & Y(:,2) < 0.7*saddle(2), 1 );
    a1tJump(i) = a1t(inxJump);
    timeJump(i) = time(inxJump);
    plot(time,Y)
    hold on
%      hold on
%        line([time(inxJump) time(inxJump)],get(gca,'ylim'))
%      hold off
%      [x,~] = ginput(1);   %select time interval to 
%      hold off
%      [~,inxJump] = sort(abs(time-x));
%      dataSelect = Y(1:inxJump,:);
%      inxJump2 = inxJump - round(10/dt);
     dataSelect = Y(1:inxJump,:);   %select data
     t_s = time(1:inxJump);
     
     %first we don't detrend it
      
     detrendData1 = dataDetrend(dataSelect(:,1),t_s,30);  %detrending methods
     detrendData2 = dataDetrend(dataSelect(:,2),t_s,30);  %detrending methods
     smoothKernel = dataSelect(:,1:2) + [detrendData1,detrendData2];
     plot(t_s,smoothKernel(:,2),'k-','LineWidth',2)
     hold off
     slidWin = round(length(dataSelect)/2);
     GAP = 100;
%      [variance_slide,corrCLE{i},lagAuto_slide,timeSelect] = slidingWidownAnalysis([dataSelect(:,1),dataSelect(:,2)],slidWin,GAP,dt);
     [variance_slide,corrCLE{i},lagAuto_slide,timeSelect,cv_slide,Fano_slide] = slidingWidownAnalysis([detrendData1,detrendData2],slidWin,GAP,dt,smoothKernel);
     variance1{i} = variance_slide(:,1);
     ewsDetr.kv1(i) = corr((1:length(variance1{i}))',variance1{i},'type','Kendall');
     
     variance2{i} = variance_slide(:,2);
     ewsDetr.kv2(i) = corr((1:length(variance1{i}))',variance2{i},'type','Kendall');
     
     coeffVar1{i} = cv_slide(:,1);
     ewsDetr.kcv1(i) = corr((1:length(variance1{i}))',coeffVar1{i},'type','Kendall');
     
     coeffVar2{i} = cv_slide(:,2);
     ewsDetr.kcv2(i) = corr((1:length(variance1{i}))',coeffVar2{i},'type','Kendall');
     
     
     FanoFactor1{i} = Fano_slide(:,1);
     ewsDetr.kf1(i) = corr((1:length(variance1{i}))',FanoFactor1{i},'type','Kendall');
     
     FanoFactor2{i} = Fano_slide(:,2);
     ewsDetr.kf2(i) = corr((1:length(variance1{i}))',FanoFactor2{i},'type','Kendall');
     
     lagAuto1{i} = lagAuto_slide(:,1);
     ewsDetr.klag1(i) = corr((1:length(variance1{i}))',lagAuto1{i},'type','Kendall');
     
     lagAuto2{i} = lagAuto_slide(:,2);
     ewsDetr.klag2(i) = corr((1:length(variance1{i}))',lagAuto2{i},'type','Kendall');
     
     ewsDetr.kcurr(i) = corr((1:length(variance1{i}))',corrCLE{i}(:,1),'type','Kendall');
     
     allTimeSelect{i} = timeSelect;
     alla1{i} = a1t(ismember(time,timeSelect));  % this seems does not work well
    
%      %using moment expansion method proposed by Luolan Chen
%      lagTime = 50;
%      allMoment = monmentFun(dataSelect(:,1:2),lagTime);
%      detrendData3 = dataDetrend(allMoment,t_s,30);
%      slidWin = round(length(detrendData3)/2);
%      GAP = 10;
%      [variance_slide_moment{i},corrCLE_moment{i},lagAuto_slide_moment{i},timeSelect]= slidingWidownAnalysis(detrendData3,slidWin,GAP,dt);
end

%strore all the relevant data into struct
ewsDetr.a1Jump = a1tJump;
ewsDetr.tJump = timeJump;
ewsDetr.var1 = variance1;
ewsDetr.var2 = variance2;
ewsDetr.cv1 = coeffVar1;
ewsDetr.cv2 = coeffVar2;
ewsDetr.corr = corrCLE;
ewsDetr.lag1 = lagAuto1;
ewsDetr.lag2 = lagAuto2;
ewsDetr.Fano1 = FanoFactor1;
ewsDetr.Fano2 = FanoFactor2;
ewsDetr.time = allTimeSelect;
ewsDetr.a1chg = alla1;


% currentFolder = pwd;  
% saveFile = [currentFolder,filesep,'figure and data',filesep,'toggSwiTimeSerialDetrend_N',...
%         num2str(OMIGA),'_se',num2str(param(7)),'_tauc',num2str(param(6)),'.mat'];
% saveFile = ['toggSwiTimeSerialDetrend_N',num2str(OMIGA),'_se',num2str(param(7)),'_tauc',num2str(param(6)),'.mat'];
% save(saveFile,'-struct','ewsDetr')
%plot the results
close all
figure(1)
xlabel('time','FontSize',30,'FontWeight','Bold')
ylabel('variance','FontSize',30,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')
hold all
figure(2)
xlabel('time','FontSize',30,'FontWeight','Bold')
ylabel('correlation coefficient','FontSize',30,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')
hold all
figure(3)
xlabel('time','FontSize',30,'FontWeight','Bold')
ylabel('lag-one autocorrelation','FontSize',30,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')
hold all
for i = 1:N;
    figure(1)
    plot(allTimeSelect{i},variance2{i},'LineWidth',2)
    figure(2)
    plot(allTimeSelect{i},corrCLE{i},'LineWidth',2)
    figure(3)
    plot(allTimeSelect{i},lagAuto2{i},'LineWidth',2)
end

save('exampleTimeDetrIntri.mat','-struct','ewsDetr')
end

function testChenModel(N,time_change)
        %this function test Chen's model, a single variable mdoel
        
        %parameters
        as = -6;
        ae = 6;
        dt = 0.001;
        tauc = 5;
        s = 0.0; %extrinsic noise
        
        
        iniVal = 2.2;
        sigma = 1.5; %intrinsic noise
        
        
        %declare some varibles
        variance1 = cell(N,1);
        corrCLE = cell(N,1);
        lagAuto1 = cell(N,1);

        allTimeSelect = cell(N,1);
        variance_slide_moment = cell(N,1);
        corrCLE_moment = cell(N,1);
        lagAuto_slide_moment = cell(N,1);
        
for i = 1:N;
    [time,Y] = LEChenModel(as,ae,tauc,time_change,iniVal,dt,sigma,s);
     plot(time,Y)
     hold on
     [x,~] = ginput(1);   %select time interval to 
     hold off
     [~,inx1] = sort(abs(time-x));
     dataSelect = Y(1:inx1,:);   %select data
     t_s = time(1:inx1);

     slidWin = round(length(dataSelect)/2);
     GAP = 50;
     [variance_slide,corrCLE{i},lagAuto_slide,timeSelect] = slidingWidownAnalysis([dataSelect(:,1),dataSelect(:,2)],slidWin,GAP,dt);
     variance1{i} = variance_slide(:,1);
     lagAuto1{i} = lagAuto_slide(:,1);
     allTimeSelect{i} = timeSelect;
     
%       plot((1:length(variance1{1}))',variance1{1})
%       plot((1:length(lagAuto1{1}))',lagAuto1{1})
     
     %using moment expansion method proposed by Luolan Chen
     lagTime = 200;
     allMoment = monmentFun(dataSelect(:,1),lagTime);
     plot((1:length(allMoment))',allMoment)
%      detrendData3 = dataDetrend(allMoment,t_s,30);
     detrendData3 = dataDetrend(allMoment,(1:length(allMoment))',30);
%      detrendData3 = allMoment;
     slidWin = round(length(detrendData3)/2);
     GAP = 50;
     [variance_slide_moment{i},corrCLE_moment{i},lagAuto_slide_moment{i},timeSelect2]= slidingWidownAnalysis(detrendData3,slidWin,GAP,dt);
     
     
     %DNB proposed by Chen
     DNBindice = mean(sqrt(variance_slide_moment{1}),2).*abs(corrCLE_moment{1}(:,1));
     %figure one
     figure(1)
     subplot(3,2,1)
     plot(timeSelect,sqrt(variance1{1}))
     
     subplot(3,2,2)
     plot(timeSelect,lagAuto1{1})
     subplot(3,2,3)
     plot(timeSelect2,sqrt(variance_slide_moment{1}))
     subplot(3,2,4)
     plot(timeSelect2,lagAuto_slide_moment{1})
     subplot(3,2,5)
     plot(timeSelect2,corrCLE_moment{1})
     
     
     figure(2)
     plot(timeSelect2,DNBindice)
 end
end

function ChenModelHist()

        timeRun = 100;
        dt = 0.005;
        tauc = 5;
        s = 0.0; %extrinsic noise      
        iniVal = 2.2;
        
        sigma = 1; %intrinsic noise
    p = [0.0,0.6,0.8,0.9];
    selectData = cell(length(p),2); %all and conditional distribution
    figure(1)
    
    for i0 = 1:length(p);
        [time,Y] = LEChenModelFixp(p(i0),tauc,timeRun,iniVal,dt,sigma,s);
        plot(time,Y)
        hold on
        [x,~] = ginput(1);   %select time interval to 
        hold off
        [~,inx1] = sort(abs(time-x));
        selectData{i0,1}(:,1) = time;
        selectData{i0,1}(:,2) = Y(:,1);
        selectData{i0,2}(:,1) = time(1:inx1);
        selectData{i0,2}(:,2) = Y(1:inx1,1);   %select data

%         subplot(1,length(p),i0)
%         hist(Y(:,1),100)   
    end
    
    %plot
    figure(2)
    for i0 = 1:length(p)
        subplot(2,length(p),i0)
        hist(selectData{i0,1}(:,2),100)
        subplot(2,length(p),length(p) + i0)
        hist(selectData{i0,2}(:,2),100,'Color','r')
    end
end

function [time,Y,a1t] = CLEtoggleExtrChgPara(param,as,ae,time_change,inV,dt,OMIGA)

a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc = param(6);
s = param(7);  %the amplitude of extrinsic noise

sigmae = sqrt(2*s*s/tauc);

N = round(time_change/dt);       %total simulation steps
X = zeros(N,3);
Xtemp = [inV;0];         %revised on Dec 31,2015
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);

%a1 changes with time
a1t = (as + (ae-as)*(0:1:round(time_change/dt))*dt/time_change)';

for i1 = 1:round(time_change/dt);
    determin = [a1t(i1)*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-(1+Xtemp(3))*Xtemp(1);...
        a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-(1+Xtemp(3))*Xtemp(2)];
    stoch = [sqrt(abs(a1t(i1)*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs((1+Xtemp(3))*Xtemp(1))),0,0;...
            0,0,sqrt(abs(a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs((1+Xtemp(3))*Xtemp(2)))];
%     determin = [a1t(i1)*OMIGA*(a0 + (1+Xtemp(3))*(1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-Xtemp(1);...
%         a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-(1+Xtemp(3))*Xtemp(2)];
%     stoch = [sqrt(abs(a1t(i1)*OMIGA*(a0 + (1+Xtemp(3))*(1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs(Xtemp(1))),0,0;...
%             0,0,sqrt(abs(a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs((1+Xtemp(3))*Xtemp(2)))];
    Xtemp([1 2]) = Xtemp([1 2]) + determin*dt + stoch*randn(4,1)*sqrt(dt);
    if(Xtemp(3)*mu+ Dext*randn>-1)
        Xtemp(3) = Xtemp(3)*mu+ Dext*randn;
    else
        Xtemp(3) = -1;
    end
    X(i1,:) = Xtemp;
end

time = (0:dt:N*dt)';
Y = [[inV;0]';X];

end
function [variance_slide,corr_slide,lagAuto_slide,timeSelect,cv_slide,Fano_slide] = slidingWidownAnalysis(data,slidWin,GAP,dt,varargin)

% if calculate the Fanofactor

[ROW,COL] = size(data);
numPoint = floor((ROW - slidWin)/GAP+1); %

lagAuto_slide = zeros(numPoint,COL);
variance_slide = zeros(numPoint,COL);

corr_slide = zeros(numPoint,COL*(COL+1)/2);



INITIL = mod((ROW - slidWin),GAP);
timeSelect = zeros(numPoint,1);
for i0 = 1:numPoint;
    Vector = data(INITIL+1+(i0-1)*GAP:(i0-1)*GAP+INITIL+slidWin,:);
    timeSelect(i0) = ((i0-1)*GAP+INITIL+slidWin)*dt;
    variance_slide(i0,:) = var(Vector,0,1);
    
    
    count = 1;
    for j1 = 1:(COL-1);
        for j2 = (j1+1):COL;
            CORR1 = corrcoef(Vector(:,j1),Vector(:,j2));
            corr_slide(i0,count) = CORR1(1,2);
            count = count +1;
        end
    end
        
    for j3 =1:COL;
        CORR2 = corrcoef(Vector(round(1/dt)+1:end,j3),Vector(1:end-round(1/dt),j3));
        lagAuto_slide(i0,j3) = CORR2(1,2);
    end
   
end

if nargin > 4
   Fano_slide = zeros(numPoint,COL);
   cv_slide = zeros(numPoint,COL);
   soomthKernel = varargin{1};
   for i0 = 1:numPoint;
       Vector = soomthKernel(INITIL+1+(i0-1)*GAP:(i0-1)*GAP+INITIL+slidWin,:);
       Fano_slide(i0,:) = variance_slide(i0,:)./mean(Vector,1);
       cv_slide(i0,:) = sqrt(variance_slide(i0,:))./mean(Vector,1);
   end    
end

end
function allMoment = monmentFun(data,lagTime)
    [ROW,COL] = size(data);
    %the length of moment should be the length of original substract lag
    %time, revised on 10/16/2016
    
    temp = ones(lagTime,1)/lagTime;
    %convolution method
    cmethod = 'valid';
    if strcmp(cmethod,'valid')
        allMoment = zeros(ROW-lagTime + 1,COL*(COL+3)/2);
    elseif(strcmp(cmethod,'same'))
        allMoment = zeros(ROW,COL*(COL+3)/2);
    end
    for i1 = 1:COL;
        allMoment(:,i1) = conv(data(:,i1),temp,cmethod);    %only use valid method will it be equivalent to sliding sum    
    end
    

    count = COL + 1;
    for j1 = 1:COL;
        for j2 = j1:COL;
        	allMoment(:,count) = conv(data(:,j1).*data(:,j2),temp,cmethod)-conv(data(:,j1),temp,cmethod).*conv(data(:,j2),temp,cmethod);
            count = count +1;
        end
    end
        
end

function [time,Y] = LEChenModel(as,ae,tauc,time_change,inV,dt,sigma,s)

%model
%dx/dt = -p + 3*(1+E)*x - x^3 + sigma

%the amplitude of extrinsic noise
sigmae = sqrt(2*s*s/tauc);

N = round(time_change/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = [inV;0];         %revised on Dec 31,2015
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);

%a1 changes with time
a1t = as + (ae-as)*(0:1:round(time_change/dt))*dt/time_change;

for i1 = 1:round(time_change/dt);
    determin = -a1t(i1) + 3*(1+Xtemp(2))*Xtemp(1) - Xtemp(1)^3;
    Xtemp(1) = Xtemp(1) + determin*dt + sqrt(sigma)*randn*sqrt(dt);
    if(Xtemp(2)*mu+ Dext*randn>-1)
        Xtemp(2) = Xtemp(2)*mu+ Dext*randn;
    else
        Xtemp(2) = -1;
    end
    X(i1,:) = Xtemp;
    
end

time = (0:dt:N*dt)';
Y = [[inV;0]';X];

end


function [time,Y] = LEChenModelFixp(p,tauc,time,inV,dt,sigma,s)
% Chen's sigle variale model with fixed p
%model
%dx/dt = -p + 3*(1+E)*x - x^3 + sigma

%the amplitude of extrinsic noise
sigmae = sqrt(2*s*s/tauc);

N = round(time/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = [inV;0];         %revised on Dec 31,2015
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);

%a1 changes with time
% a1t = as + (ae-as)*(0:1:round(time/dt))*dt/time;

for i1 = 1:round(time/dt);
    determin = -p + 3*(1+Xtemp(2))*Xtemp(1) - Xtemp(1)^3;
    Xtemp(1) = Xtemp(1) + determin*dt + sqrt(sigma)*randn*sqrt(dt);
    if(Xtemp(2)*mu+ Dext*randn>-1)
        Xtemp(2) = Xtemp(2)*mu+ Dext*randn;
    else
        Xtemp(2) = -1;
    end
    X(i1,:) = Xtemp;
end

time = (0:dt:N*dt)';
Y = [[inV;0]';X];

end

function [mean_jump_time,std_jump_time] = JumpTimeChenModel(N)
   %this function simulate the jump events at large noise
   %and plot the histgram of the corresponding time(instaneous value of parameter)
   p1 = -3;
   p2 = 3;
   tauc = 10;
   time_change = 5e2;
   dt = 2e-3;
   inV = 2.2;
   sigma = [0.01,0.1,1,5];
   s = 0;
   altState = -1;
   
   t_jump = nan(N,length(sigma));
   for j0 = 1:length(sigma);
        for i0 = 1:N;
            t_jump(i0,j0) = LEChenModelJumpTime(p1,p2,tauc,time_change,inV,dt,sigma(j0),s,altState);
        end
   end
   
   mean_jump_time = mean(t_jump,1);
   std_jump_time = std(t_jump,0,1);
   
   figure(1)
   for i0 = 1:length(sigma)
        subplot(1,4,i0)
        hist(t_jump(:,i0),100)
        xlabel('p','FontSize',20,'FontWeight','Bold')
        ylabel('Frequency','FontSize',20,'FontWeight','Bold')
        set(gca,'Xlim',[p1,p2],'LineWidth',1,'FontSize',16)
   end
   
   figure(2)
   errorbar(sigma,mean_jump_time,std_jump_time,'ro','MarkerSize',15,'LineWidth',3)
   xlabel('sigma','FontSize',24,'FontWeight','Bold')
   ylabel('Average jump point','FontSize',24,'FontWeight','Bold')
end

function t_jump = LEChenModelJumpTime(p1,p2,tauc,time_change,inV,dt,sigma,s,altState)

%this funciton return when the state jump to alternative state
    %the amplitude of extrinsic noise
sigmae = sqrt(2*s*s/tauc);

N = round(time_change/dt);       %total simulation steps
X = zeros(N,2);
Xtemp = [inV;0];         %revised on Dec 31,2015
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);
threshold = 0.1; % a threshold used to define if it enters the althernative state

%a1 changes with time
pt = p1 + (p2-p1)*(0:1:round(time_change/dt))*dt/time_change;

for i1 = 1:round(time_change/dt);
    determin = -pt(i1) + 3*(1+Xtemp(2))*Xtemp(1) - Xtemp(1)^3;
    Xtemp(1) = Xtemp(1) + determin*dt + sqrt(sigma)*randn*sqrt(dt);
    if(Xtemp(2)*mu+ Dext*randn>-1)
        Xtemp(2) = Xtemp(2)*mu+ Dext*randn;
    else
        Xtemp(2) = -1;
    end
    X(i1,:) = Xtemp;
    
    if(abs(Xtemp(1)-altState)< threshold)
        t_jump = pt(i1);
        return
    end
    
end

end
function inxJump = toggTimeJump(time,Y,low)
% return when it jumped to alternative state
    inxJump = find(Y(:,1) > 1.05*low(1) && Y(:,2) < 0.95*low(2));

end
