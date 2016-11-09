a2 = 10;
beta = 2;
gamma = 2;
sigma1 = 0.02;
sigma2 = 0.02;
T = 400;         %total simulation time
dt = 0.05;
N = T/dt;       %total simulation steps
numP = 101;      %number of parameters to scan
NUMPATHS = 25;

Xzero = [1;10]; %inital values
X = zeros(N,numP);
Y = zeros(N,numP);

aVAR1 = zeros(numP,1);
aVAR2 = zeros(numP,1);
aCORR = zeros(numP,1);
par = zeros(1,4);
par(2) = a2;
par(3) = beta;
par(4) = gamma;
a1 = zeros(numP,1);

for i1 = 1:numP;
    tspan = [0 T];
    iniVal = Xzero;
    par(1) = 1 + 0.5*(i1 - 1);
    a1(i1) = par(1);
    parameter = [par(1), a2, beta, gamma];
    [t,y] = ode45(@toggleSwitch,tspan,iniVal,[],parameter);
    uss = y(length(y),1);
    vss = y(length(y),2);
    D = [sigma1*sigma1,0;0,sigma2*sigma2];
    covariance = analyCov(uss,vss,par,D);
    aVAR1(i1) = covariance(1,1);
    aVAR2(i1) = covariance(1,1);
    aCORR(i1) = covariance(1,2)/sqrt(covariance(1,1))/sqrt(covariance(2,2));
end