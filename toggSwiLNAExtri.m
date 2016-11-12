function toggSwiLNAExtri()
% In this program, you can compare different source of extrinsic nosie
%

% default parameters
sigmae = 0.1;
tauc = 1000;
param = [15;10;2;2;0.03;tauc;sigmae];
OMIGA = 20;        %system size
Yini = [0,10];
% try different amplitude
paraInx = [2,5];
extrVal = [0.1,1000;0.1,1000];
% allTau = [1,10,100];


earlyWarningSum = toggSwiLNA(OMIGA,param,Yini,paraInx,extrVal);

%plot the figure
plotLNA(earlyWarningSum);

end

function earlyWarningSum = toggSwiLNA(OMIGA,param,Yini,paraInx,extrVal)
%
% this program use linear noise approximation to get the early warning
% signals of toggle switch model subject to both intrinsic and extrinsic
% noise, here I try to integrate previous fucntions into one so that
% different source, amplitude and also correlation time scale of extrinsic
% noise can be handled
% it returns a data struct containing the early warning signals: mean
% value, variance, correlation time scale
% specifying the positions of extrinsic noise and its corresponding tauc
% and amplitude
% last revised on 11/11/2016
 if length(paraInx) ~=  size(extrVal,1);
     error(message('number if extrinsic noise and input values has to be the same'))
 end
a1_s = 4;
a1_e = 16;
a0 = param(5);
N_a1 = 40;
runTime = 1000;


% ss= [1e-3,1e-2,1e-1,0.25];
% sigmae = 2*ss.^2*param(6);
sigmae = 2*extrVal(:,1).^2./extrVal(:,2);
corrcoef_ext_cle_theory = zeros(N_a1,1);
variance_ext_cle_theory = zeros(N_a1,2);
lagOneCLEext_theory = zeros(N_a1,2);
variance_ext_c_theory = zeros(N_a1,2);
corrcoef_ext_c_theory = zeros(N_a1,1);
lagOneC_theory = zeros(N_a1,2);

for j1 = 1:N_a1;
    a1_list(j1) = a1_s + (a1_e - a1_s)*j1/N_a1;
    param(1)= a1_list(j1);
    [td,Yd]=odeRun(param,runTime,Yini);
    steadyODE(j1,:)=Yd(end,:);
    for j2 = 1:length(sigmae)

        [Ac,Acle] = jacMatrx(param,steadyODE(j1,:),OMIGA,paraInx,extrVal);
        Dc = [a1_list(j1)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,2)^param(3)))+steadyODE(j1,1)*OMIGA,0;...
                0,param(2)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,1)^param(4)))+steadyODE(j1,2)*OMIGA];
        if(size(extrVal,1) == 1)
            Dcle = [a1_list(j1)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,2)^param(3)))+steadyODE(j1,1)*OMIGA,0,0;...
                0,param(2)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,1)^param(4)))+steadyODE(j1,2)*OMIGA,0;0,0,sigmae];
        elseif(size(extrVal,1) == 2)
            Dcle = [a1_list(j1)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,2)^param(3)))+steadyODE(j1,1)*OMIGA,0,0,0;...
                0,param(2)*OMIGA*(a0 + (1-a0)/(1+steadyODE(j1,1)^param(4)))+steadyODE(j1,2)*OMIGA,0,0;0,0,sigmae(1),0;0,0,0,sigmae(2)];
        else
            error('only two sources if extrinsic')
        end
        
        C3 = lyap(Ac,Dc);            
        corrcoef_ext_c_theory(j1) = C3(1,2)/(sqrt(C3(1,1)*C3(2,2)));
        variance_ext_c_theory(j1,2) = C3(2,2);
        variance_ext_c_theory(j1,1) = C3(1,1);
        G3 = C3*expm(Ac');
        lagOneC_theory(j1,2) = G3(2,2)/C3(2,2);
        lagOneC_theory(j1,1) = G3(1,1)/C3(1,1);
        
        C4 = lyap(Acle,Dcle);
        corrcoef_ext_cle_theory(j1) = C4(1,2)/(sqrt(C4(1,1)*C4(2,2)));
        
        variance_ext_cle_theory(j1,2) = C4(2,2);
        variance_ext_cle_theory(j1,1) = C4(1,1);
        G4 = C4*expm(Acle');
        lagOneCLEext_theory(j1,2) = G4(2,2)/C4(2,2);
        lagOneCLEext_theory(j1,1) = G4(1,1)/C4(1,1);    

    end
end
    earlyWarningSum = struct('a1',[],'meanVal',[],'var0',[],'corr0',[],'lag0',[],...
        'var',[],'corr',[],'lag',[]);
    earlyWarningSum.a1 = a1_list';
    earlyWarningSum.meanVal = steadyODE;
    earlyWarningSum.var0 = variance_ext_c_theory;
    earlyWarningSum.corr0 = corrcoef_ext_c_theory;
    earlyWarningSum.lag0 = lagOneC_theory;
    earlyWarningSum.var = variance_ext_cle_theory;
    earlyWarningSum.corr = corrcoef_ext_cle_theory;
    earlyWarningSum.lag = lagOneCLEext_theory;
end

function [Ac,Ace] = jacMatrx(param,vec,OMIGA,paraInx,extrVal)
% x y E are variables
% A is the jacobian matrix for the deterministic equation
% Ae is the jacobian matrix for the model with extrinsic noise
% modified the reaction with basal synthesis rate
syms x y E E1 E2 a1 a2 a0 n1 n2 b1 b2 N
    v1 = N*a1*(a0 + (1-a0)/(1+(y/N)^n1));
    v2 =  -x;
    v3 = N*a2*(a0 + (1-a0)/(1+(x/N)^n2));
    v4 = -y;

if (length(paraInx)==1)
    switch paraInx
        case {1,2}
            v1 = N*a1*(a0 + (1-a0)*(1+E1)/(1+(y/N)^n1));
            v3 = N*a2*(a0 + (1-a0)*(1+E1)/(1+(x/N)^n2));
            beta1 = 1/param(6);
        case {3,4}
            v1 = N*a1*(a0 + (1-a0)/(1+E1+(y/N)^n1));
            v3 = N*a2*(a0 + (1-a0)/(1+E1+(x/N)^n2));
            beta1 = 1/param(6);
        case {5,6}
            v2 = -(1 + E1)*x;
            v4 = -(1 + E1)*y;
            beta1 = 1/extrVal(2);
        otherwise,
            error('the range of index has to be 1 ~ 6 !')
        
    end
    J4 = jacobian([v1 + v2;v3+v4;-b1*E1],[x y E1]);
    Ace = double(subs(J4,{x,y,E1,a1 a2 n1 n2 a0 b1 N},...
         {vec(1)*OMIGA,vec(2)*OMIGA,0,param(1),param(2),param(3),param(4),param(5),beta1,OMIGA}));

else
    
    flag = 1; %since we only consider two extrinsic noise situation
    beta1 = 1/param(6);
    beta2 = 1/param(6);
    if(any(paraInx==1))
        v1 = N*a1*(a0 + (1-a0)*(1+E)/(1+(y/N)^n1));
        flag = 0;          
    end 
    
    if(any(paraInx==2))
         if flag
            v3 = N*a2*(a0 + (1-a0)*(1+E1)/(1+(x/N)^n2));
            flag = 0;
        else
            v3 = N*a2*(a0 + (1-a0)*(1+E2)/(1+(x/N)^n2));
         end
    end
    
    if(any(paraInx==3))
         if flag
            v1 = N*a1*(a0 + (1-a0)/(1+E2+(y/N)^n1));
            flag = 0;
         else
            v1 = N*a1*(a0 + (1-a0)/(1+E2+(y/N)^n1));
         end
    end
    
    
     if(any(paraInx==3))
         if flag
            v1 = N*a1*(a0 + (1-a0)/(1+E2+(y/N)^n1));
            flag = 0;
         else
            v1 = N*a1*(a0 + (1-a0)/(1+E2+(y/N)^n1));
         end
     end
     
     
     if(any(paraInx==4))
         if flag
            v3 = N*a2*(a0 + (1-a0)/(1+E1+(x/N)^n2));
            flag = 0;
         else
            v3 = N*a2*(a0 + (1-a0)/(1+E2+(x/N)^n2));
         end 
     end
     
     
     if(any(paraInx==5))
         if flag
            v2 = -(1 + E1)*x;
            beta1 = 1/extrVal(1,2);
            flag = 0;
         else
            v2 = -(1 + E2)*x;
            beta1 = 1/extrVal(2,2);
         end
     end
     
     
     if(any(paraInx==6))
         if flag
            v4 = -(1 + E1)*y;
            beta2 = 1/extrVal(1,2);
            flag = 0;
         else
            v4 = -(1 + E2)*y;
            beta2 = 1/extrVal(2,2);
         end           
     end
    J4 = jacobian([v1 + v2;v3+v4;-b1*E1;-b2*E2],[x y E1 E2]);
    Ace = double(subs(J4,{x,y,E1,E2,a1 a2 n1 n2 a0 b1 b2 N},...
        {vec(1)*OMIGA,vec(2)*OMIGA,0,0,param(1),param(2),param(3),param(4),param(5),beta1,beta2,OMIGA}));
end
J3 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - y],[x y]);
Ac = double(subs(J3,{x,y,a1 a2 n1 n2 a0 b1 N},{vec(1)*OMIGA,vec(2)*OMIGA,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));
% J4 = jacobian([a1*N/(1+E+(y/N)^n1) - x;a2*N/(1+E+(x/N)^n2) - y;-b*E],[x y E]);
% J4 = jacobian([N*a1*(a0 + (1-a0)/(1+(y/N)^n1)) - (1+E)*x;N*a2*(a0 + (1-a0)/(1+(x/N)^n2)) - (1+E)*y;-b*E],[x y E]);
% J4 = jacobian([N*a1*(a0 + (1+E)*(1-a0)/(1+(y/N)^n1)) - x;N*a2*(a0 + (1+E)*(1-a0)/(1+(x/N)^n2)) - y;-b*E],[x y E]);
% J4 = jacobian([N*a1*(a0 + (1-a0)/(1+E+(y/N)^n1)) - x;N*a2*(a0 + (1-a0)/(1+E+(x/N)^n2)) - y;-b*E],[x y E]);
% Ace = double(subs(J4,{x,y,E,a1 a2 n1 n2 a0 b N},{vec(1)*OMIGA,vec(2)*OMIGA,0,param(1),param(2),param(3),param(4),param(5),param(6),OMIGA}));

end

function [t,y]=odeRun(param,totTime,Yini)

tspan = [0 totTime];
[t,y] = ode45(@toggleSwitch,tspan,Yini,[],param);
% plot(t,y,'LineWidth',2)
% legend('gene1','gene2')

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

function plotLNA(earlyWarningSum)
% input is a data struct that containing all the information of signals
a1_list = earlyWarningSum.a1;
figure(1)
hold on
plot(a1_list,earlyWarningSum.corr0,'k-','Linewidth',3)
plot(a1_list,earlyWarningSum.corr,'Linewidth',3,'LineStyle','--')
hold off


figure(2)
hold on
plot(a1_list,earlyWarningSum.var0,'k-','Linewidth',3)
plot(a1_list,earlyWarningSum.var,'Linewidth',3,'LineStyle','--')
hold on

figure(3)
hold on
plot(a1_list,earlyWarningSum.lag0,'k-','Linewidth',3)
plot(a1_list,earlyWarningSum.lag,'Linewidth',3,'LineStyle','--')    
hold off

end