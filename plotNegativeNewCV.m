%this function replot the coefficient of variance
load('/media/M_fM__VM_0M_eM__JM__M_eM__MM_7/MATLAB/criticalityData/simu_varCoefLag_tau1e2_se0.1-0101.mat');
notNanInx = ~isnan(mean(meanVal2,2));
variance2New = variance2(notNanInx,:);
for i = 1:size(variance2New,1)
    list = [];
    for j = 1:size(variance2New,2)
        if(variance2New(i,j)<5e3)
            list  = [list,variance2New(i,j)/meanVal2(i,j)];
        end
    end
    cvNew(i,:) = [mean(list),std(list)];
end
hold on
errorbar(a1list(notNanInx)',cvNew(:,1),cvNew(:,2))