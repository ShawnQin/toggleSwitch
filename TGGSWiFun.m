function F = TGGSWiFun(x,a1)
a2 = 10;
beta = 3;
gamma = 2;
F = [a1/1+x(2)^beta - x(1);a2/1+x(1)^gamma - x(2)];