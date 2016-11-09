% stochastic integration of toggle swich model using Euler-Maruyama method
% stochasitc equations:
%----------------------------------------------------------------------
%
% du = a1/1+v^beta - u + sigma1*dW1
% dv = a2/(1+ u^gamma) - v + sigma2*dW2
%
% with paramers a2 = 10, beta = 3, gamma = 3, sigma1 = sigma2 = 0.5;
clear

a1 = 5.8;
a2 = 10;
beta = 2;
gamma = 2;
sigma1 = 0.2;
sigma2 = 0.2;
T = 400;         %total simulation time
dt = 0.05;
N = T/dt;       %total simulation steps
Xzero = [30;1]; %inital values
X = zeros(N,2);
% dW = sqrt(dt)*randn(N,2); % a N by 2 matrix with normal distribution elements
% W = cumsum(dW);            %cummulative add of dW
Xtemp = Xzero;
dT = [dt;dt];
SIGMA = [sigma1,0;0,sigma2];
for i1 = 1:N;
    a1(i1) = 2 + 0.02*i1*dt;
    TEMP = [a1(i1)/(1+Xtemp(2)^beta)-Xtemp(1);a2/(1+Xtemp(1)^gamma)-Xtemp(2)];
    Xtemp = Xtemp + TEMP.*dT + SIGMA*(sqrt(dT).*randn(2,1));
    X(i1,:) = Xtemp;
end
[variance, lagAutocorrelation ,corrFluc] =covAutoCorr( X(:,1),X(:,2) );
%plot the stochastic trajactories
plot((0:dt:T)',[Xzero';X],'LineWidth',2)
legend('gene1','gene2')
xlabel('time','FontSize',24,'FontName','calibri')
ylabel('gene expression','FontSize',24,'FontName','calibri')
hold on 

% using deterministic method to calculate 
tspan = [0 T];
iniVal = Xzero;
parameter = [a1, a2, beta, gamma];
[t,y] = ode45(@toggleSwitch,tspan,iniVal,[],parameter);
plot(t,y,'LineWidth',2)
legend('gene1','gene2')
%usign Euler method to calculate
XtempEu = Xzero;
for i2 = 1:N;
    stepTemp = [a1/(1+XtempEu(2)^beta)-XtempEu(1);a2/(1+XtempEu(1)^gamma)-XtempEu(2)];
    XtempEu = XtempEu + stepTemp.*dT;
    Xeuler(i2,:) = XtempEu;
end
plot((0:dt:T)',[Xzero';Xeuler])
hold off
