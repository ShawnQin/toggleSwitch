%% replot figure 3
fh1 = figure(1);
set(fh1,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
set(gca,'ylim',[0,10])
print('-dpdf','negative_compare_var_tau10_degrad_1117.pdf')


fh2 = figure(2);
set(fh2,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
set(gca,'ylim',[-0.3,1])
print('-dpdf','negative_compare_corr_tau10_degrad_1117.pdf')


fh3 = figure(3);
set(fh3,'Units','inches','Position',[0 0 8 6],...
'PaperPositionMode','auto');
set(gca,'ylim',[-0.2,1])
print('-dpdf','negative_compare_lag_tau10_degrad_1117.pdf')
