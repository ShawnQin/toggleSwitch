%**********************************************************************
% this grogram repeat the work done be Micheal Assaf in 2013, published on
% PRL. This is a very simple auto activation gene circuit
% The authors discussed the effect of extrinsic noise in detail 
% Here, I simluate this model using modiied Gillespie algorithm
% to show how extrinsic nose may affect the distribution, bistable region
% and even te auto correlation function 
% last revised on Dec. 19,2015
%************************************************************************
function selfactivationGlp()
 
%parameters
a0 = 0.05;
OMIGA = 300;  %system size
x0list = [0.2:0.05:0.4,0.42:0.02:0.56,0.6,0.7];
x0 = 0.42;
%set the initial value
tspan = [0 1000];
y0 = 10;
[t,y] = ode15s(@odeSelfActi,tspan,y0,[],a0,OMIGA,x0);
inV = round(y(end));
totTime = 20000;
tauc = 1000;
sigmae = 2*(0.02*inV)^2/tauc;
[time,Y]=YGlpExtri(inV,OMIGA,a0,x0,totTime,sigmae,tauc);
plot(time,Y);

end
function dydt = odeSelfActi(t,y,a0,OMIGA,x0)
        dydt = a0*OMIGA + (1-a0)*OMIGA*y^2/(y^2+ (OMIGA*x0)^2) - y;
end
function [time,Y]=YGlpExtri(inV,OMIGA,a0,x0,totTime,sigmae,tauc)
%ths is a modified Gillespie algorithm
%here the sigame is not the same as in the paper

%first, produce all the extrinsic noise\
dt = 0.05;
Next = round(1.1*totTime/dt);
mu = exp(-dt/tauc);
Dext = sigmae*sqrt((1-mu^2)*tauc/2);
Xtemp = 0;
Xe = zeros(Next,1);
for i1 = 1:Next;
    Xtemp = Xtemp*mu+ Dext*randn;
    Xe(i1) = Xtemp;
end

Y(1,1) = inV;
Y(1,2) = 0;  % the second column store extrinsic noise
S = [1,-1];
deg = 1;
t = 0;
te = dt;
i = 1;
j = 0;
while(t < totTime)
    a = [OMIGA*(a0 + (1-a0)*Y(i,1)^2/(Y(i,1)^2+(OMIGA*x0)^2)); deg*Y(i,1)];
    [tau,inx] = min(-1./a.*log(rand(2,1)));
    if(t+tau<te)
        Y(i+1,1) = Y(i,1) + S(:,inx);
        Y(i+1,2) = deg-1;
        t = t + tau;  
        time(i+1)= t;
        i = i+1;
    else
        deg = 1+Xe(j+1);
        if(deg<0)
            deg = 0;%negtive degradation is impossible
        end
        t = te;
        te = te+dt;   %updata the extrinsic time point
        j = j+1;
    end
    
end
time = [0;time(1:end-1)'];
end