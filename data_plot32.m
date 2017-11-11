%%
idx = 2:NUM_EVENT;
[~,i] = min(abs(options.Theta-1));      % choose theta close to 1.

%% Number of Reconfiguration
if exist('fig_num_reconfig', 'var') && fig_num_reconfig.isvalid
    figure(fig_num_reconfig);
else
    fig_num_reconfig = figure('Name', 'Number of Reconfiguration');
end
xlabel('Experiment Time');
% fig_num_reconfig.OuterPosition = [100 400 400 380];
yyaxis left;
hl = plot(idx, results.Reconfig.NumberReconfigVariables(idx), 'g-.', ...
    idx, results.Fastconfig(i).NumberReconfigVariables(idx), 'r-',...
    idx, results.Fastconfig2(i).NumberReconfigVariables(idx), 'b--');
hl(1).Parent.YLim = [0 900];
hl(1).Parent.XLim = [idx(1) idx(end)];
ylabel('Number of Reconfiguration');
yyaxis right
% hr = plot(idx, g_results.slice.num_flows(idx), 'r--');
hr = plot(idx, results.Reconfig.NumberVariables(idx), 'k:');
ylabel('Number of Variables');
%hr.Parent.YLim = [2500 4500];
hr.Parent.YLim = [0 900];
hr.Parent.YColor = 'k';
% title(sprintf('\\theta=%.2f', options.Theta(i)));
legend({'Reconfig', 'ReF', 'ReF2', 'Variables'}, 'Location', 'northwest');

%% Ratio of Reconfiguration
if exist('fig_ratio_reconfig', 'var') && fig_ratio_reconfig.isvalid
    figure(fig_ratio_reconfig);
else
    fig_ratio_reconfig = figure('Name', 'Ratio of Reconfiguration');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
ratio_reconfig = results.Reconfig(i).NumberReconfigVariables(idx)./...
    results.Reconfig(i).NumberVariables(idx);
ratio_fastconfig = results.Fastconfig(i).NumberReconfigVariables(idx)./...
    results.Fastconfig(i).NumberVariables(idx);
ratio_fastconfig2 = results.Fastconfig2(i).NumberReconfigVariables(idx)./...
    results.Fastconfig2(i).NumberVariables(idx);
hl = plot(idx, results.Reconfig.RatioReconfigVariables(idx), 'b--', ...
    idx, results.Fastconfig(i).RatioReconfigVariables(idx), 'g-',...
    idx, results.Fastconfig2(i).RatioReconfigVariables(idx), 'r-.');
ylabel('Ratio of Reconfiguration');
xlabel('Experiment Time');
hl(1).Parent.XLim = [idx(1) idx(end)];
% hl(1).Color = [0, 0.447, 0.471];
% hl(2).Color = [0, 0, 1];        % blue
legend({'Reconfig', 'ReF', 'ReF2'}, 'Location', 'best');
% title(sprintf('\\theta=%.2f', options.Theta(i)));

%% Cost of Reconfiguration
if exist('fig_cost_reconfig', 'var') && fig_cost_reconfig.isvalid
    figure(fig_cost_reconfig);
else
    fig_cost_reconfig = figure('Name', 'Cost of Reconfiguration');
end
h = plot(idx, results.Reconfig.Cost(idx)*options.Theta(i), 'g--', ...
    idx, results.Fastconfig(i).Cost(idx), 'b-', ...
    idx, results.Fastconfig2(i).Cost(idx), 'r-.');
h(1).Parent.XLim = [idx(1) idx(end)];
ylabel('Cost');
xlabel('Experiment Time');
legend({'Reconfig', 'ReF', 'ReF2'}, 'Location', 'northwest');
% title(sprintf('\\theta=%.2f', options.Theta(i)));

%% Profit
if exist('fig_profit_reconfig', 'var') && fig_profit_reconfig.isvalid
    figure(fig_profit_reconfig);
else
    fig_profit_reconfig = figure('Name', 'Profit with Reconfiguration');
end
% fig_profit_reconfig.OuterPosition = [100 400 400 380];
xlabel('Experiment Time');
yyaxis left;
profit = results.Reconfig.Profit + (1-options.Theta(i))*results.Reconfig.Cost;
hl = plot(idx, profit(idx), 'g-.', idx, results.Fastconfig(i).Profit(idx), 'b--',...
    idx, results.Fastconfig2(i).Profit(idx), 'r-');
hl(1).Parent.YLim = [3000 4000];
ylabel('Profit');
yyaxis right;
hr = plot(idx, results.Reconfig.NumberFlows(idx), 'k:');
hr.Parent.YLim = [0 100];
hl(1).Parent.XLim = [idx(1) idx(end)];
ylabel('Number of Flows');
legend({'Reconfig', 'ReF', 'ReF2', 'Flow'}, 'Location', 'southeast');
% title(sprintf('\\theta=%.2f', options.Theta(i)));