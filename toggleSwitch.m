function dydt = toggleSwitch(t,y,parameter)

% sythestic toggle switch model adopted form James J.Collins' 2000 Nature
% paper

a1 = parameter(1);          %Gene1 sythesis rate
a2 = parameter(2);          %Gene2 sythesis rate
beta = parameter(3);        %hill coefficient
gamma = parameter(4);       %hill coefficient

dydt = zeros(size(y));
u = y(1);
v = y(2);
dudt = a1/(1+v^beta) - u;
dvdt = a2/(1+u^gamma) - v;
dydt = [dudt;dvdt];
end



