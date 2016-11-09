
clear
clc


fileName = 'simu_varCoefLag_tau1e_se0.25-1231.mat';
load(fileName)

dt = 0.1;
t_s = (0:0.1:1e3)';
for i= 1:length(a1list);
    for j = 1:size(dataSelect,2);
        if(~isempty(dataSelect{i,j}))   %revised on Jan 1,2016
            detrendData1 = dataDetrend(dataSelect{i,j}(:,1),t_s,10);  %detrending methods
            detrendData2 = dataDetrend(dataSelect{i,j}(:,2),t_s,10);  %detrending methods
            meanVal1(i,j) = mean(detrendData1);
            meanVal2(i,j) = mean(detrendData2);            
            variance1(i,j) = var(detrendData1);
            variance2(i,j) = var(detrendData2);
            C1 = corrcoef(detrendData1,detrendData2);
            corrCLE(i,j) = C1(1,2);
            
            C3 = corrcoef(detrendData1(1:end-round(1/dt),1),detrendData1(round(1/dt)+1:end,1));
            lagAuto1(i,j) = C3(1,2);
            C4 = corrcoef(detrendData2(1:end-round(1/dt)),detrendData2(round(1/dt)+1:end));
            lagAuto2(i,j) = C4(1,2);

        end
    end
end

    notNanInx = ~isnan(mean(meanVal2,2));
    figure(1)
%     plot(dataB(:,1),OMIGA*dataB(:,2))
    hold on
    errorbar(a1list(notNanInx)',mean(meanVal2(notNanInx,:),2),std(meanVal2(notNanInx,:),0,2))
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('mean expression level','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    hold off
    
    figure(2)
    hold on
    errorbar(a1list(notNanInx)',mean(variance2(notNanInx,:),2),std(variance2(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(variance2(notNanInx,:),2),[std(variance2(notNanInx,:),0,2),std(variance2(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('variance','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    
    figure(3)
    hold on
    errorbar(a1list(notNanInx)',mean(corrCLE(notNanInx,:),2),std(corrCLE(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(corrCLE(notNanInx,:),2),[std(corrCLE(notNanInx,:),0,2),std(corrCLE(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('correlation coefficient of noise','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')
    
    figure(4)
    hold on
    errorbar(a1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),std(lagAuto2(notNanInx,:),0,2))
%     boundedline(a1list(notNanInx)',mean(lagAuto2(notNanInx,:),2),[std(lagAuto2(notNanInx,:),0,2),std(lagAuto2(notNanInx,:),0,2)],'alpha')
    xlabel('a1','FontSize',24,'FontWeight','Bold')
    ylabel('lag one auto correlation','FontSize',24,'FontWeight','Bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','Bold')