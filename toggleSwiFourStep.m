function dydt = toggleSwiFourStep(t,y,n,lamda,beta,OMIGA)
% dydt = zeros(y);

dydt = [lamda(1)/(1+(y(4)/OMIGA)^n)- beta(1)*y(1); lamda(2)*y(1)-...
         beta(2)*y(2); lamda(3)/(1+(y(2)/OMIGA).^n)- beta(3)*y(3);...
         lamda(4)*y(3)- beta(4)*y(4)];
end