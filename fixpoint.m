function answer = fixpoint( parameter )
% equations to solve:
% -------------------------------------------------
% du = a1/1+v^beta - u
% dv = a2/(1+ u^gamma) - v
% -------------------------------------------------
a1 = parameter(1);
a2 = parameter(2);
beta = parameter(3);
gamma = parameter(4);

syms x
 ANS = solve(x*a1^beta/(1+x^beta)^gamma + x - a2 == 0, 'Real', true);
 answer(:,1) = sort(double(ANS),1);
 answer(:,2) = a1./(1+answer(:,1).^beta);
end