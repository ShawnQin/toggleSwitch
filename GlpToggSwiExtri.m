%-----------------------------------------------------------------------
% PROGRAM NAME: GlpToggSwiExtri.m
% DESCIPTION:  
%   this program using Gillespie algorithum to simulate the mean swithching
%   rate of a toggle switch model
%   writen by Shanshan Qin@Tanglab, PKU
%   last revised on Jan.3,2016
%-------------------------------------------------------------------------
function GlpToggSwiExtri()
close all
clear
clc

sigmae = 0.1;
tauc = 100;
param = [11;10;2;2;0.03;tauc;sigmae];
OMIGA = 20;          %system size
NUM = 50;

a1list = 12:0.5:15;
switchTime(param,a1list,NUM,OMIGA);

end
function t_swi= GlpSwitchExtrinsic(iniVal,param,range,OMIGA)
a1 = param(1);
a2 = param(2);
n1 = param(3);
n2 = param(4);
a0 = param(5);
tauc = param(6);
ss = param(7);  %the amplitude of extrinsic noise

sigmae = sqrt(2*ss^2/tauc);

dt = 0.02;
Ytemp = iniVal;
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);

MAXTIME = 1e6;
Next = round(1.1*MAXTIME/dt);
Xtemp = 0;

Xe = zeros(Next,1);
% for i1 = 1:Next;
%     Xtemp = Xtemp*mu+ Dext*randn;
%     Xe(i1) = Xtemp;
% end

S = [1,-1,0,0;0,0,1,-1];
a1t = a1;   % a1 is the parameter with extrinsic noise
a2t = a2;
b = 1;  %the basical death rate
t = 0;
te = dt;
j = 0;
Z = [];
while(t < MAXTIME)
    a = [OMIGA*a1t*(a0 + (1-a0)/(1+(Ytemp(2)/OMIGA)^n1)); b*Ytemp(1);OMIGA*a2t*(a0 + (1-a0)/(1+(Ytemp(1)/OMIGA)^n2)); b*Ytemp(2)];
    [tau,inx] = min(-1./a.*log(rand(4,1)));
    if(t+tau<te)
        Ytemp = Ytemp + S(:,inx);  
        t = t + tau; 
        Z = [Z;Ytemp'];
    else
%         a1t = a1*(1+Xe(j+1));
%         a2t = a2*(1+Xe(j+1));
        Xtemp = Xtemp*mu+ Dext*randn;
        if(Xtemp>-1)
            b = Xtemp + 1;
        else
            b = -1;
        end
%         b = 1+Xe(j+1);   %extrinsic noise added on the degradation term
        t = te;
        te = te+dt;   %updata the extrinsic time point
        j = j+1;
    end
    
    if(Ytemp(1)<range(1,2) && Ytemp(1)>range(1,1) && Ytemp(2)<range(2,2) && Ytemp(2)>range(2,1))
         t_swi = t;
         return
    end
    
end

t_swi = nan;   %if no swithing happen, return nan

end
function switchTime(param,a1list,N,OMIGA)
%first load the file contain the steady state information

fileA = 'toggSwiBifurAll_A.csv';
fileB = 'toggSwiBifurAll_B.csv';
dataA = csvread(fileA);
dataB = csvread(fileB);

for j1 = 1:length(a1list);
    inx_h(j1) = find(dataB(:,1) ==a1list(j1));
end

high = round(OMIGA*[dataA(inx_h,2),dataB(inx_h,2)]);
low = round(OMIGA*[dataA(inx_h,4),dataB(inx_h,4)]);
min_saddle = round(OMIGA*[min(dataA(inx_h,3)),min(dataA(inx_h,3))]);

t_swi_h2l = zeros(N,length(high));
MAXTRY = 2*N;   % maximum trying times
for i0 = 1:length(high); 
    param(1) = a1list(i0);
    % firt from high to low state
    lowRef = [low(i0,1) - 20,low(i0,1)+min([low(i0,1) + 20,min_saddle(1)]);low(i0,2) - 20,low(i0,2)+min([low(i0,2) + 20,min_saddle(2)])]; %the border of lower state
%     highRef = [high(i0,1) - 20,high(i0,1)+min([high(i0,1) + 20,min_saddle(1)]);high(i0,2) - 20,high(i0,2)+high([low(i0,2) + 20,min_saddle(2)])]; %the border of higher state
    n = 0;
    flag = 0;
    while(n<N && flag<MAXTRY)
        ts = GlpSwitchExtrinsic(high(i0,:)',param,lowRef,OMIGA);
        if(~isnan(ts))
            n = n+1;
            t_swi_h2l(n,i0) = ts;
        end
            flag = flag + 1;
    end
end
saveFileName = ['swiTime-GLP-','tau',num2str(param(6)),'-','se',num2str(param(7)),'-',date,'.mat'];
save(saveFileName,'t_swi_h2l','a1list')
end