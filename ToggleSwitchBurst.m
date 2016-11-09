cd('/Users/ningwang/Documents/work/correlation of fluctuation/burst toggle switch');
clear;
close all
%set an error
options=odeset('RelTol',1e-8,'AbsTol',1e-11);
color=[0 0 0;1 0 0;0 1 0;0 0 1];
%timespan

TIME=1*10^6;
tspan = [0,TIME];
flag=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial value
global namada1 namada2 namada3 namada4 beta1 beta2 beta3 beta4 alpha1 alpha2 k gama;
namada1=1/50;
namada3=1/50;
beta1=1/500;
beta3=1/500;
namada2=1/100;
namada4=1/100;
beta2=1/10000;
beta4=1/10000;
alpha1=2.5;
alpha2=3;
k=500;
gama=1;

Y01=[0,1000,0,0];
Y02=[0,0,1000,0];
M=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ser=[0.4 0.5 0.6 0.7 0.71 0.72 0.721 0.722 0.723 0.724 0.725 0.73 0.75 0.9 1 1.1 1.15 1.16 1.17 1.18 1.182 1.1826 1.183 1.186 1.19 1.2 1.3 1.4 1.5]
%ser=[0.4 0.5]; 
for m=1:length(ser)
    gama=ser(m);
    flag(m)=gama;
    [t,y] = ode15s(@BTS,tspan,Y01,options);
    
    ys1=[];
    ys2=[];
    ys3=[];
    ys4=[];
    Ts=[];
    
    for j=1:M
y1(1)=round(y(end,1));
y2(1)=round(y(end,2));
y3(1)=round(y(end,3));
y4(1)=round(y(end,4)); 
 
N=1*10^6; % Initial values and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T(1)=0;


for i=1:N % N simulation runs
rates=[namada1/(1+(y4(i)/k)^alpha1) beta1*y1(i) namada2*y1(i) beta2*y2(i) namada3/(1+(y2(i)/(k*gama)).^alpha2) beta3*y3(i) namada4*y3(i) beta4*y4(i)];
kT=sum(rates);
u_time=rand;
u_event=rand;
t(i)=-log(1-u_time)/kT; % Exponentially distributed number
T(i+1)=T(i)+t(i);

 if  u_event < rates(1)/kT
       y1(i+1)=y1(i)+1;   % Birth chosen
 elseif  rates(1)/kT < u_event && u_event < (rates(1)+rates(2))/kT
       y1(i+1)=y1(i)-1;   % Death chosen
 else
     y1(i+1)=y1(i);
 end
 
 
 if  (rates(1)+rates(2))/kT < u_event && u_event < (rates(1)+rates(2)+rates(3))/kT
       y2(i+1)=y2(i)+1;   % Birth chosen
 elseif  (rates(1)+rates(2)+rates(3))/kT < u_event && u_event < (rates(1)+rates(2)+rates(3)+rates(4))/kT
       y2(i+1)=y2(i)-1;   % Death chosen
 else
     y2(i+1)=y2(i);
 end

 
 if  (rates(1)+rates(2)+rates(3)+rates(4))/kT < u_event && u_event < (rates(1)+rates(2)+rates(3)+rates(4)+rates(5))/kT
       y3(i+1)=y3(i)+1;   % Birth chosen
 elseif  (rates(1)+rates(2)+rates(3)+rates(4)+rates(5))/kT < u_event && u_event < (rates(1)+rates(2)+rates(3)+rates(4)+rates(5)+rates(6))/kT
       y3(i+1)=y3(i)-1;   % Death chosen
 else
     y3(i+1)=y3(i);
 end

  if  (rates(1)+rates(2)+rates(3)+rates(4)+rates(5)+rates(6))/kT< u_event && u_event < (rates(1)+rates(2)+rates(3)+rates(4)+rates(5)+rates(6)+rates(7))/kT
       y4(i+1)=y4(i)+1;   % Birth chosen
 elseif  (rates(1)+rates(2)+rates(3)+rates(4)+rates(5)+rates(6)+rates(7))/kT< u_event && u_event < 1
       y4(i+1)=y4(i)-1;   % Death chosen
 else
     y4(i+1)=y4(i);
 end
 
end

ys1(j,:)=y1;
ys2(j,:)=y2;
ys3(j,:)=y3;
ys4(j,:)=y4;
TS(j,:)=T;
%%%%%%%%%%%%%%
        Mean2step(j,m)=mean(ys2(j,:));
        Var2step(j,m)=var(ys2(j,:));
        
        r=corrcoef(ys2(j,:),ys4(j,:));
        Corr2step(j,m)=r(1,2);
        
        a=autocorr(ys2(j,:),1);
        Auto2step(j,m)=a(2);
end

       
  waitbar(m/length(ser),'Simulation in Process:');

 end
save(strcat('ToggleSwitchburst',num2str(M),'repeats',num2str(N),'steps'));

