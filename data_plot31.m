%% 
if exist('fig_stat_reconfig', 'var') && fig_stat_reconfig.isvalid
    figure(fig_stat_reconfig);
else
    fig_stat_reconfig = figure('Name', 'Satistics of Reconfiguration');
end
idx = 51:NUM_EVENT;
x = 1:length(options.theta);
avgrat_reconfig = mean(results.reconfig.rat_reconfig(idx))*ones(length(x),1);
avgrat_fastconfig = zeros(length(x),1);
for i = x
    avgrat_fastconfig(i) = mean(results.fastconfig(i).rat_reconfig(idx));
end
plot(x, avgrat_reconfig, '-o', x, avgrat_fastconfig, '-x');
ylim([0, 0.4]);
legend({'Reconfiguration', 'ReF'});
xlabel('Experiment Time');
ylabel('Avrage reconfigration ratio');