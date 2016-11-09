%-----------------------------------------------------------------------
% PROGRAM NAME: toggleCLESwitime.m
% DESCIPTION:  
%   this program using CLE method to simulate the mean switching time of
%   toggle switch model with extrinsic noise
% writen by Shanshan Qin@Tanglab, PKU
% last revised on Jan,3rd,2016
%----------------------------------------------------------------------
function toggleCLESwitime()
% close all
clear
clc


a1list = 11:1:15;
sigmae = 0.25;
tauc = 1000;
param = [15;10;2;2;0.03;tauc;sigmae];
OMIGA = 20;        %system size
Nswi = 500;       %number of running when calculate the switching rate


switchTime(param,a1list,Nswi,OMIGA);
end

function switchTime(param,a1list,N,OMIGA)

%first load the file contain the steady state information
fileA = 'toggSwiBifurAll_A.csv';
fileB = 'toggSwiBifurAll_B.csv';
dataA = csvread(fileA);
dataB = csvread(fileB);

% inx_h = 8:1:12;
for j1 = 1:length(a1list);
    inx_h(j1) = find(dataB(:,1) ==a1list(j1));
end
high = OMIGA*[dataA(inx_h,2),dataB(inx_h,2)];
low = OMIGA*[dataA(inx_h,4),dataB(inx_h,4)];
min_saddle = OMIGA*[min(dataA(inx_h,3)),min(dataA(inx_h,3))];
% a1 = a1list(inx_h);

t_swi = zeros(N,length(high));
MAXTRY = 2*N;   % maximum trying times
for i0 = 1:length(high); 
    param(1) = a1list(i0);
    % firt from high to low state
    lowRef = [low(i0,1) - 20,low(i0,1)+min([low(i0,1) + 20,min_saddle(1)]);low(i0,2) - 20,low(i0,2)+min([low(i0,2) + 20,min_saddle(2)])]; %the border of lower state
    n = 0;
    flag = 0;
    while(n<N && flag<MAXTRY)
        ts = CLEswitch_ext(high(i0,:)',param,lowRef,OMIGA);
        if(~isnan(ts))
            n = n+1;
            t_swi(n,i0) = ts;
        end
            flag = flag + 1;
    end
end
saveFile = ['swiTime_CLE_','tau',num2str(param(6)),'_se',num2str(param(7)),'_',date,'.mat'];
save(saveFile,'t_swi','a1list')
end

function t_swi = CLEswitch_ext(inV,param,range,OMIGA)
%range define the vicinity of a state
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc = param(6);
s = param(7);  %the amplitude of extrinsic noise

sigmae = sqrt(2*s*s/tauc);

dt = 0.1;
Xtemp = [inV;0];
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);

TMAX = 1e6; %maximun running time
    time = 0;
    while(time < TMAX)
        determin = [a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1))-(1+Xtemp(3))*Xtemp(1);...
        a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2))-(1+Xtemp(3))*Xtemp(2)];
        stoch = [sqrt(abs(a1*OMIGA*(a0 + (1-a0)/(1+(Xtemp(2)/OMIGA)^n1)))),-sqrt(abs((1+Xtemp(3))*Xtemp(1))),0,0;...
            0,0,sqrt(abs(a2*OMIGA*(a0 + (1-a0)/(1+(Xtemp(1)/OMIGA)^n2)))),-sqrt(abs((1+Xtemp(3))*Xtemp(2)))];
        Xtemp([1 2]) = Xtemp([1 2]) + determin*dt + stoch*randn(4,1)*sqrt(dt);
        Xtemp(3) = Xtemp(3)*mu+ Dext*randn;
        if(Xtemp(3)<-1)
            Xtemp(3) = -1; %avoid unreasonable value
        end
    
        time = time + dt;
        if(Xtemp(1)<range(1,2) && Xtemp(1)>range(1,1) && Xtemp(2)<range(2,2) && Xtemp(2)>range(2,1))
            t_swi = time;
            return
        end
    end
    t_swi = nan;   %if no swithing happen, return nan

end