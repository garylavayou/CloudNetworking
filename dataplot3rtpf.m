%% Performance Evaluation of Dual-ADMM Method with the Normal Convex Optimization Method
% Dual-ADMM: distributed parallel implementation;
% Normal: interior-point method provided by the bulit-in <fmincon> function of Matlab.

%% Running time of ADMM and the normal method.
load('./Results/singles/EXP6002.mat', 'runtime')
ll = length(runtime.varsize);
mean_time_normal = zeros(ll,1);
for i = 1:ll
    mean_time_normal(i) = mean(runtime.varsize(i).normal);
end
load('./Results/singles/EXP6002e1w10p4.mat', 'runtime', 'numberflow')
mean_time_admm = zeros(ll,1);
for i = 1:ll
    mean_time_admm(i) = mean(runtime.varsize(i).admm)/numberflow(i)*6;
end
load('./Results/singles/EXP6002e1w10p4li.mat', 'runtime', 'numberflow')
mean_time_admm_limit = zeros(ll,1);
for i = 1:ll
    mean_time_admm_limit(i) = mean(runtime.varsize(i).admm_limit)/numberflow(i)*6;
end
if exist('fig_runtime', 'var') && fig_runtime.isvalid
    figure(fig_runtime);
    fig_runtime.Children.delete;
else
    fig_runtime = figure('Name', 'Runtime');
end
fig_runtime.OuterPosition(3:4) = [360 375];  % [496 476];
hl = plot(numberflow, mean_time_normal, '-^', ...
    numberflow, mean_time_admm, '-o',...
    numberflow, mean_time_admm_limit, '--x');
color_set = [Color.MildBlue; Color.MildGreen; Color.Red];
for k=1:length(hl)
    hl(k).LineWidth = 1;
    hl(k).Color = color_set(k).RGB;
end
legend({'Interior-Point', 'Dual-ADMM', 'Dual-ADMM-EB'}, 'Location', 'northwest');
xlabel('Number of flows');
ylabel('Run time (seconds)');
xlim([100,1000]);
export_fig(fig_runtime, './Figures/runtime-varsize', '-pdf', '-transparent');
%% Reconfigration Ratio
load('./Results/singles/EXP6002.mat', 'results')
ll = length(results.Fastconfig);
mean_rr_normal = zeros(ll,1);
for i = 1:ll
    mean_rr_normal(i) = mean(results.Fastconfig{i}.NumberReconfigVariables...
        ./results.Fastconfig{i}.NumberVariables);
end
load('./Results/singles/EXP6002e1w10p4.mat', 'results')
mean_rr_admm = zeros(ll,1);
for i = 1:ll
    mean_rr_admm(i) = mean(results.Fastconfig{i}.NumberReconfigVariables...
        ./results.Fastconfig{i}.NumberVariables);
end
load('./Results/singles/EXP6002e1w10p4li.mat', 'results')
mean_rr_admm_limit = zeros(ll,1);
for i = 1:ll
    mean_rr_admm_limit(i) = mean(results.Fastconfig{i}.NumberReconfigVariables...
        ./results.Fastconfig{i}.NumberVariables);
end
if exist('fig_reconfigratio', 'var') && fig_reconfigratio.isvalid
    figure(fig_reconfigratio);
    fig_reconfigratio.Children.delete;
else
    fig_reconfigratio = figure('Name', 'Reconfiguration Ratio');
end
fig_reconfigratio.OuterPosition(3:4) = [360 375];  % [496 476];
hl = plot(numberflow, mean_rr_normal, '-^', ...
    numberflow, mean_rr_admm, '-o',...
    numberflow, mean_rr_admm_limit, '--x');
color_set = [Color.MildBlue; Color.MildGreen; Color.Red];
for k=1:length(hl)
    hl(k).LineWidth = 1;
    hl(k).Color = color_set(k).RGB;
end
legend({'Interior-Point', 'Dual-ADMM', 'Dual-ADMM-EB'}, 'Location', 'northeast');
xlabel('Number of flows');
ylabel('Reconfiguration ratio');
xlim([100,1000]);
export_fig(fig_reconfigratio, './Figures/reconfigratio-distr-varsize', '-pdf', '-transparent');
