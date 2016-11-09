% stochastic integration of toggle swich model using Euler-Maruyama method
% this program scaning parameter of protein sythesis rate a1 or a1 and plot
% the correlation of fluctuatioin with a1(or a2) to check wether it has
% critiality at bifurcation point(tipping point)
% stochasitc equations:
%----------------------------------------------------------------------
%
% du = a1/1+v^beta - u + sigma1*dW1
% dv = a2/(1+ u^gamma) - v + sigma2*dW2
%
% with paramers a2 = 10, beta = 3, gamma = 3, sigma1 = sigma2 = 0.5 fixed;

clear
%prameter assignment
a2 = 10;
beta = 2;
gamma = 2;
sigma1 = 0.04;
sigma2 = 0.04;
T = 400;         %total simulation time
dt = 0.04;
N = T/dt;       %total simulation steps
numP = 101;      %number of parameters to scan
NUMPATHS = 25;

Xzero = [1;10]; %inital values
X = zeros(N,NUMPATHS);
Y = zeros(N,NUMPATHS);
varGene1 = zeros(numP,NUMPATHS);
varGene2 = zeros(numP,NUMPATHS);
lagGene1 = zeros(numP,NUMPATHS);
lagGene2 = zeros(numP,NUMPATHS);
CORRFLUC = zeros(numP,NUMPATHS);


dT = [dt;dt];
SIGMA = [sigma1;sigma2];
% for j0 = 1:NUMPATHS;        % total paths simulated
%     for i0 = 1:numP;
%         a1(i0) = 1 + 0.5*(i0 - 1);
% %     Xtemp1 = Xzero(1);
% %     Xtemp2 = Xzero(2);
%         Xtemp = Xzero;
%         for i1 = 1:N;
%             TEMP = [a1(i0)/(1+Xtemp(2)^beta)-Xtemp(1);a2/(1+Xtemp(1)^gamma)-Xtemp(2)];
%             Xtemp = Xtemp + TEMP.*dT + SIGMA.*(sqrt(dT).*randn(2,1));
%             X(i1,i0) = Xtemp(1);
%             Y(i1,i0) = Xtemp(2);
%         end
%     end
% [variance, lagAutocorrelation ,corrFluc] =covAutoCorr( X , Y );
% varGene1(j0,:) = variance(1,:);
% varGene2(j0,:) = variance(2,:);
% lagGene1(j0,:) = lagAutocorrelation(1,:);
% lagGene2(j0,:) = lagAutocorrelation(2,:);
% CORRFLUC(j0,:) = corrFluc;
% 
% end
Xstor = zeros(N,1);
Ystor = zeros(N,1);
for i0 = 1:numP;
    a1(i0) = 1 + 0.5*(i0 - 1);
    INDEX = 0;
    TRY = 0;
    DET = determine(a1(i0),Xzero,T);
    
    while(INDEX<NUMPATHS && TRY<=200);                 % at most 100 trajectories
        Xtemp = Xzero;
        for i1 = 1:N;
            TEMP = [a1(i0)/(1+Xtemp(2)^beta)-Xtemp(1);a2/(1+Xtemp(1)^gamma)-Xtemp(2)];
            Xtemp = Xtemp + TEMP.*dT + SIGMA.*(sqrt(dT).*randn(2,1));
            Xstor(i1,1) = Xtemp(1);
            Ystor(i1,1) = Xtemp(2);
        end
        STO = mean([Xstor(N-1001:N,1),Ystor(N-1001:N,1)]);  % the last 1000 points are used
  
        if (max(abs(STO-DET)./DET) <=0.5 || max(abs(STO-DET)) <=0.5 )
            X(:,INDEX+1) = Xstor;
            Y(:,INDEX+1) = Ystor;
            INDEX = 1+INDEX;
        end
        TRY = TRY + 1;
    end
% numerical results
[variance, lagAutocorrelation ,corrFluc] =covAutoCorr( X , Y );
varGene1(i0,:) = variance(1,:);
varGene2(i0,:) = variance(2,:);
lagGene1(i0,:) = lagAutocorrelation(1,:);
lagGene2(i0,:) = lagAutocorrelation(2,:);
CORRFLUC(i0,:) = corrFluc;

end
% mean value and standard deviation
VAR1 = mean(varGene1,2);  %modified
ERR1 = std(varGene1,0,2);
VAR2 = mean(varGene2,2);
ERR2 = std(varGene1,0,2);
LAG1 = mean(lagGene1,2);
ERR3 = std(lagGene1,0,2);
LAG2 = mean(lagGene2,2);
ERR4 = std(lagGene2,0,2);
CORRFN = mean(CORRFLUC,2);
ERR5 = std(CORRFLUC,0,2);

%-------------------------------------------------------------------------
% analytical results
aVAR1 = zeros(numP,1);
aVAR2 = zeros(numP,1);
aCORR = zeros(numP,1);
lagAutoCor1 = zeros(numP,1);
lagAutoCor2 = zeros(numP,1);
par = zeros(4,1);
par(2) = a2;
par(3) = beta;
par(4) = gamma;
a1 = zeros(numP,1);
for i1 = 1:numP;
    tspan = [0 T];
    iniVal = Xzero;
    par(1) = 1 + 0.5*(i1 - 1);
    a1(i1) = par(1);
    parameter = [a1(i1), a2, beta, gamma];
    [t,y] = ode45(@toggleSwitch,tspan,iniVal,[],parameter);
    uss = y(length(y),1);
    vss = y(length(y),2);
    D = [sigma1*sigma1,0;0,sigma2*sigma2];
    [covariance,lagCov] = analyCov(uss,vss,par,D,'saddlenode');
    aVAR1(i1) = covariance(1,1);
    aVAR2(i1) = covariance(2,2);
    aCORR(i1) = covariance(1,2)/sqrt(covariance(1,1))/sqrt(covariance(2,2));
    lagAutoCor1(i1) = lagCov(1,1);
    lagAutoCor2(i1) = lagCov(2,2);
end
%----------------------------------------------------------------------------
%plot variance and correlation of fluctuation over a1
figure(1)
errorbar(a1,VAR1,ERR1,'g-o')
hold on 
errorbar(a1,VAR2,ERR2,'g-*')
legend('gene1','gene2')
plot(a1,aVAR1,'r+',a1,aVAR2,'c+')
xlabel('a1(sythesis rate)','FontSize',24,'FontName','calibri')
ylabel('variance','FontSize',24,'FontName','calibri')
hold off

figure(2)
errorbar(a1,CORRFN,ERR5,'b-o')
hold on
plot(a1,aCORR,'r+')
xlabel('a1(sythesis rate)','FontSize',24,'FontName','calibri')
ylabel('correlation of fluctuation','FontSize',24,'FontName','calibri')
hold off

figure(3)
errorbar(a1,LAG1,ERR3,'r-o')
hold on 
errorbar(a1,LAG2,ERR4,'g-*')
legend('gene1','gene2')
plot(a1,lagAutoCor1,'r<',a1,lagAutoCor2,'g+')
xlabel('a1(sythesis rate)','FontSize',24,'FontName','calibri')
ylabel('lag 1 autocorrealtion','FontSize',24,'FontName','calibri')
hold off

% using deterministic method to calculate 

