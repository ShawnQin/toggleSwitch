% toggle switch ode model
% NAME: tggSwiRun
clear
% fprintf('Please input the parameter array,1*4,:\n')
parameter = input('Please input the parameter array,1*4,:\n');

tspan = [0 100];
upp = max(parameter(1),parameter(1));
for i1 = 1:20;
    iniVal = [i1*upp/20,i1*upp/20+0.1];
% iniVal = [1;1];
    [t,y] = ode45(@toggleSwitch,tspan,iniVal,[],parameter);
    plot(t,y,'LineWidth',1), hold on
end
xlabel('time','FontSize',24,'FontName','Calibri');
ylabel('gene expression','FontSize',24,'FontName','Calibri');
hold off