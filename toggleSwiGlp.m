%this program simulates the toggle switch with
%Gillespie method
%model description
%**************************************************************
% A |--|B   %mutual inhibiton 

%**************************************************************
function [time,Y] = toggleSwiGlp(totTime,inV,param,OMIGA)
lamda1 = 0.3;
lamda2 = 1/2;
lamda3 = 2;
lamda4 = 1/2;
beta1 = 1/10;
beta2 = 1/50;
beta3 = 1/10;
beta4 = 1/50;
n = 2;                      %hill coefficient set as the same
% OMIGA = 100;                %system size
NUM = 100000;                %total number point
% time = zeros(NUM+1,1);
% Y = zeros(NUM+1,4);
% Y(1,:)= [10,80,80,20];       %initial value
Y(1,:) = inV;
% Y1(1) = 100;
% Y2 = zeros(NUM);
% Y2(1) = 100;
t = 0;
i = 1;

while(t < totTime && i<= NUM+1)
    a = [lamda1/(1+(Y(i,4)/OMIGA)^n); beta1*Y(i,1); lamda2*Y(i,1);...
         beta2*Y(i,2); lamda3/(1+(Y(i,2)/OMIGA).^n); beta3*Y(i,3);...
         lamda4*Y(i,3); beta4*Y(i,4)];
    a0 = sum(a);   
    t = t + 1/a0*log(1/rand);
    time(i)= t;
    reference = rand*a0;
    if(reference<a(1))
        Y(i+1,1) = Y(i,1) + 1;
        Y(i+1,2:4) = Y(i,2:4);
        i = i +1;
    elseif(reference > a(1) && (reference <= sum(a(1:2))))
        Y(i+1,1) = Y(i,1) -1;
        Y(i+1,2:4) = Y(i,2:4);
        i = i+1;
    elseif(reference > sum(a(1:2)) && (reference <= sum(a(1:3))))
        Y(i+1,2) = Y(i,2) + 1;
        Y(i+1,1) = Y(i,1);
        Y(i+1,3:4) = Y(i,3:4);
        i =i +1;
    elseif(reference > sum(a(1:3)) && (reference <= sum(a(1:4))))
        Y(i+1,2) = Y(i,2) - 1;
        Y(i+1,1) = Y(i,1);
        Y(i+1,3:4) = Y(i,3:4);
        i =i +1;
    elseif(reference > sum(a(1:4)) && (reference <= sum(a(1:5))))
        Y(i+1,3) = Y(i,3) + 1;
        Y(i+1,1:2) = Y(i,1:2);
        Y(i+1,4) = Y(i,4);
        i =i +1; 
    elseif(reference > sum(a(1:5)) && (reference <= sum(a(1:6))))
        Y(i+1,3) = Y(i,3) - 1;
        Y(i+1,1:2) = Y(i,1:2);
        Y(i+1,4) = Y(i,4);
        i =i +1; 
    elseif(reference > sum(a(1:6)) && (reference <= sum(a(1:7))))
        Y(i+1,4) = Y(i,4) + 1;
        Y(i+1,1:3) = Y(i,1:3);
        i = i+1;
    elseif(reference > sum(a(1:7)))
        Y(i+1,4) = Y(i,4) - 1;
        Y(i+1,1:3) = Y(i,1:3);
        i = i+1;
        
    end
end
time = [0;time'];
plot(time,Y)
% hold on
%--------------------------------------------------------------------------
%ode
% tspan = [0 time(length(time))];
% [t,y] = ode45(@GapGeneNetwork,tspan,Y(1,:)',[],Km,tao,n,OMIGA);
% plot(t,y)
% xlabel('Time(a.u)','FontSize',30,'FontName','Times New Roman')
% ylabel('Gene Expression(a.u)','FontSize',30,'FontName','Times New Roman')
% set(gca,'FontSize',20)
end