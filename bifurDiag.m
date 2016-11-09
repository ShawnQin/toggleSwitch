% bigurcation diagram of toggle switch model
% bifurcation parameter is a1
clear
FIX_u = zeros(100,3);        %preallocation 
FIX_v = zeros(100,3);
PRAM = zeros(100,4); 
PRAM(:,2)  = 10;                % a2 = 10;
PRAM(:,3)  = 2;                 % beta = 3;
PRAM(:,4)  = 2;                 % gamma = 2;
for i1 = 1:100;
    PRAM(i1,1) = 0.1*i1;
    TEMP = fixpoint( PRAM(i1,:));
    FIX_v(i1,:) = TEMP(:,1)';
%     FIX_u(i1,:) = i1./(1+FIX_v(i1,:).^PRAM(i1,3));
    FIX_u(i1,:) = TEMP(:,2)';
end

plot((1:1:100)',FIX_v,'ro')
legend('gene2')
xlabel('a1(gene1 sythesis rate)','FontSize',24,'FontName','calibri')
ylabel('gene1 expression','FontSize',24,'FontName','calibri')
hold on 
plot((1:1:100)',FIX_u,'g<')
legend('gene1')