%%
idx = 51:NUM_EVENT;

%% Number of Reconfiguration
if exist('fig_num_reconfig', 'var') && fig_num_reconfig.isvalid
    figure(fig_num_reconfig);
else
    fig_num_reconfig = figure('Name', 'Number of Reconfiguration');
end
xlabel('Experiment Time');
% fig_num_reconfig.OuterPosition = [100 400 400 380];
for i = 1:length(options.theta)
    yyaxis left;
    hl = plot(idx, results.reconfig.num_reconfig(idx), 'k-.', ...
        idx, results.fastconfig(i).num_reconfig(idx), 'b-');
    ylabel('Number of Reconfiguration');
    yyaxis right
    % hr = plot(idx, g_results.slice.num_flows(idx), 'r--');
    hr = plot(idx, results.reconfig.num_vars(idx), 'r--');
    ylabel('Number of Variables');
    %hr.Parent.YLim = [2500 4500];
    hr.Parent.YLim(1) = 0;
    title(sprintf('\\theta=%.2f', options.theta(i)));
    legend({'Reconfig', 'ReF', 'Variables'}, 'Location', 'northwest');
    pause(1);
end

%% Ratio of Reconfiguration
if exist('fig_ratio_reconfig', 'var') && fig_ratio_reconfig.isvalid
    figure(fig_ratio_reconfig);
else
    fig_ratio_reconfig = figure('Name', 'Ratio of Reconfiguration');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
% yyaxis left;
for i = 1:length(options.theta)
    hl = plot(idx, results.reconfig.rat_reconfig(idx), '-.', ...
        idx, results.fastconfig(i).rat_reconfig(idx), '-');
    ylabel('Ratio of Reconfiguration');
    xlabel('Experiment Time');
    hl(1).Color = [0, 0.447, 0.471];
    hl(2).Color = [0, 0, 1];        % blue
    % yyaxis right
    % hr = plot(idx, g_results.slice.num_flows(idx), 'r:',...
    %     idx, num_nonzeros_reconfig(idx), '-.',...
    %     idx, num_nonzeros_fastconfig(idx), '--');
    % hr(2).Color = [0.871, 0.49, 0];
    % hr(3).Color = [0.749, 0, 0.749];
    % ylabel('Number of Reconfiguration');
    % legend({'Ratio-Reconfig', 'Ratio-FastConfig', 'Flows', 'Reconfig', 'Fast Reconfig'}, ...
    %     'Location', 'northwest');
    legend({'Reconfig', 'ReF'}, 'Location', 'northwest');
    title(sprintf('\\theta=%.2f', options.theta(i)));
    pause(1);
end
%% Cost of Reconfiguration
if exist('fig_cost_reconfig', 'var') && fig_cost_reconfig.isvalid
    figure(fig_cost_reconfig);
else
    fig_cost_reconfig = figure('Name', 'Cost of Reconfiguration');
end
for i = 1:length(options.theta)
    h = plot(idx, results.reconfig.cost(idx)*options.theta(i)/options.theta(1), '-.', ...
        idx, results.fastconfig(i).cost(idx));
    ylabel('Cost');
    xlabel('Experiment Time');
    legend({'Reconfig', 'ReF'}, 'Location', 'northwest');
    title(sprintf('\\theta=%.2f', options.theta(i)));
    pause(0.5);
end
%% Profit
if exist('fig_profit_reconfig', 'var') && fig_profit_reconfig.isvalid
    figure(fig_profit_reconfig);
else
    fig_profit_reconfig = figure('Name', 'Profit with Reconfiguration');
end
% fig_profit_reconfig.OuterPosition = [100 400 400 380];
xlabel('Experiment Time');
for i = 1:length(options.theta)
    yyaxis left;
    hl = plot(idx, results.reconfig.profit(idx), 'b-', ...
        idx, results.fastconfig(i).profit(idx), 'k-.');
    hl(1).Parent.YLim(1) = 0;
    ylabel('Profit');
    yyaxis right;
    hr = plot(idx, results.reconfig.num_flows(idx), 'r--');
    hr.Parent.YLim(1) = 0;
    ylabel('Number of Flows');
    legend({'Reconfig', 'ReF', 'Flow'}, 'Location', 'northwest');
    title(sprintf('\\theta=%.2f', options.theta(i)));
    pause(0.5);
end