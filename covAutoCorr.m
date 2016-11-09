function [variance, lagAutocorrelation, corrFluc] =covAutoCorr( matrx1,matrx2,dt,Ns)
% this fucntion take two matrx input and return covariance and lag 1
% autocorrelation of the variables X and Y
% if Ns is not set, then  use all the data
[ROW1,COL1] = size(matrx1);
[ROW2,COL2] = size(matrx2);
%--------------------------time series analysis, using only one trajectory
% Ns = 6000;              % total of points selected from a time seris

if nargin<3 || isempty(dt)
    fprintf('you should secify the time step!\n')
end

if nargin<4 || isempty(Ns)
    Data1 = matrx1;
    Data2 = matrx2;
    Ns = length(matrx1);

elseif (ROW1 <= Ns)
    fprintf('the data points of this stochastic trajectory is too short, better input a longer one!\n')
else
    Data1 = matrx1(ROW1-Ns + 1:ROW1,:);        %data used for statistics, assuming that the final 1000 points is stationary process
    Data2 = matrx2(ROW2-Ns + 1:ROW2,:);        %ROW1 and ROW1 should be the same
end
% define 1 time steps length 
% dt = 0.05;
NUMSTEP = round(1/dt);  % time interval is 1/dt length,so it measures lag 1 autocorrelation
FIRST1 = Data1(1:NUMSTEP,:);
FIRST2 = Data2(1:NUMSTEP,:);
END1 = Data1(Ns + 1-NUMSTEP:Ns,:);
END2 = Data2(Ns + 1-NUMSTEP:Ns,:);
TEMP1 = Data1;
TEMP2 = Data2;

TEMP1(1:NUMSTEP,:) = [];
TEMP1(Ns-2*NUMSTEP+1:Ns-NUMSTEP,:) = [];
TEMP2(1:NUMSTEP,:) = [];
TEMP2(Ns-2*NUMSTEP+1:Ns-NUMSTEP,:) = [];
Vector1 = [FIRST1;TEMP1];
Vector2 = [TEMP1;END1];
Vector3 = [FIRST2;TEMP2];
Vector4 = [TEMP2;END2];

for i1 = 1:COL1;
    CORR1 = corrcoef(Vector1(:,i1),Vector2(:,i1));
    lagAutocorrelation(1,i1) = CORR1(1,2);         %colum represent different trajactories
    CORR2 = corrcoef(Vector3(:,i1),Vector4(:,i1));
    lagAutocorrelation(2,i1) = CORR2(1,2);
    
    CORR3 = corrcoef(Data1(:,i1),Data2(:,i1));
    corrFluc(i1) = CORR3(1,2);
end
    variance(1,:) = var(Data1,0,1);
    variance(2,:) = var(Data2,0,1);
%---------------------esembles averaging, with many trajectories-----------
% only carring out when more than 15 trajectoies are computed
% if (COL1 >= 15)
%     for j1 = 1:COL1;
%         CORR4 = corrcoef(Data1(1,:),Data2(2,:));
%         corrFluc(j1) = CORR4(1,2);           %correlation of fluctuation between gen1 and two 
%         variance(j1,1) = var(Data1(j1,:));   %variance of gene1 at each time point
%         variance(j1,2) = var(Data2(j1,:));   %variance of gene2 at each time point
%     end
% end

end

