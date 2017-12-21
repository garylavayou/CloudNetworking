%% single plot and accumulate quantity
% The independent variables is time.
line_width = [1 1 1 1 1 1];
marker_size = [6 7 6 8 6 6 6];
marker = {'none', '+', 'none', 'x', 's', 'd', 'o'};
line_style = {'-.', 'none', '-', 'none', '--', '--', '--'};
color_set = [Color.MildGreen; Color.Purple; Color.MildBlue; Color.Red; Color.Black; Color.Gray];
legend_label = {'RFFV', 'RFEV', 'Hybrid-1', 'Hybrid-2', 'Baseline', '#Variables', 'Hybrid-0'};
i = 1;
tx = 1:NUM_EVENT;
t = results.Reconfig.Time(tx);
if length(tx)>=10
    marker_index = round(linspace(1,length(tx), 10));
end
%% Number of Reconfiguration (with benchmark)
if exist('fig_num_reconfig', 'var') && fig_num_reconfig.isvalid
    figure(fig_num_reconfig);
else
    fig_num_reconfig = figure('Name', 'Number of Reconfiguration (With benchmark)');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
yyaxis('left');
hl = plot(t, cumsum(results.Fastconfig{i}{tx,'NumberReconfigVariables'}), '-.',...
    t, cumsum(results.Fastconfig2{i}{tx,'NumberReconfigVariables'}), '+',...
    t, cumsum(results.Dimconfig{i}{tx,'NumberReconfigVariables'}),'-',...
    t, cumsum(results.Dimconfig2{i}{tx,'NumberReconfigVariables'}),'x',...
    t, cumsum(results.Reconfig{tx,'NumberReconfigVariables'}), '--s');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).MarkerIndices = marker_index;
    hl(k).MarkerSize = marker_size(k);
end
ylabel('Number of Reconfiguration');
ylimmax = hl(1).Parent.YLim(2);
yyaxis right
hr = plot(t, cumsum(results.Reconfig{tx,'NumberVariables'}), '--d');
hr(1).Color = Color.Gray.RGB;
hr(1).LineWidth = 1;
hr(1).MarkerIndices = marker_index;
ylabel('Number of Variables');
hr.Parent.YColor = Color.Gray.RGB;
% title(sprintf('\\eta=%.2f', thetas(i)));
legend(legend_label, 'Location', 'northwest');
ylimmax = max(ylimmax, hr.Parent.YLim(2));
yyaxis left;
ylim(hl(1).Parent, [0, ylimmax]);
yyaxis right;
ylim(hr.Parent, [0, ylimmax]);
xlabel('Time');
xlim([t(1), t(end)]);
%% Number of Reconfiguration (No benchmark)
if exist('fig_num_reconfig_nb', 'var') && fig_num_reconfig_nb.isvalid
    figure(fig_num_reconfig_nb);
else
    fig_num_reconfig_nb = figure('Name', 'Number of Reconfiguration (No benchmark)');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
hl = plot(t, cumsum(results.Fastconfig{i}{tx,'NumberReconfigVariables'}), '-.',...
    t, cumsum(results.Fastconfig2{i}{tx,'NumberReconfigVariables'}), '--+',...
    t, cumsum(results.Dimconfig{i}{tx,'NumberReconfigVariables'}),'-',...
    t, cumsum(results.Dimconfig2{i}{tx,'NumberReconfigVariables'}),'--x');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = 1;
    hl(k).MarkerIndices = marker_index;
    hl(k).MarkerSize = marker_size(k);
end
xlabel('Time');
ylabel('Number of Reconfiguration');
legend(legend_label(1:4), 'Location', 'northwest');
xlim([t(1), t(end)]);

%% Cost of Reconfiguration (with benchmark)
if exist('fig_cost_reconfig', 'var') && fig_cost_reconfig.isvalid
    figure(fig_cost_reconfig);
else
    fig_cost_reconfig = figure('Name', 'Cost of Reconfiguration (With benchmark)');
end
cost = [results.Fastconfig{i}{tx,'Cost'}, results.Fastconfig2{i}{tx,'Cost'},...
    results.Dimconfig{i}{tx,'Cost'}, results.Dimconfig2{i}{tx,'Cost'},...
    results.Reconfig{tx,'Cost'}*thetas(i)];
reconfig_type = [results.Fastconfig{i}{tx,'ReconfigType'}, results.Fastconfig2{i}{tx,'ReconfigType'},...
    results.Dimconfig{i}{tx,'ReconfigType'}, results.Dimconfig2{i}{tx,'ReconfigType'},...
    results.Reconfig{tx,'ReconfigType'}];
cost = cost/slice_opt.Flow.ArrivalRate;
for k = 1:size(cost,2)
    idx = reconfig_type(:,k)==ReconfigType.Dimensioning;
    cost(idx,k) = cost(idx,k)*slice_opt.EventInterval;
end
hl = plot(t, cumsum(cost));
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).LineStyle = line_style{k};
    hl(k).Marker = marker{k};
    hl(k).MarkerIndices = marker_index;
end
ylabel('Cummulated Cost');
xlabel('Time');
xlim([t(1),t(end)]);
legend(legend_label(1:5), 'Location', 'northwest');
%% Cost of Reconfiguration (no benchmark)
if exist('fig_cost_reconfig_nb', 'var') && fig_cost_reconfig_nb.isvalid
    figure(fig_cost_reconfig_nb);
else
    fig_cost_reconfig_nb = figure('Name', 'Cost of Reconfiguration (No benchmark)');
end
cost = [results.Fastconfig{i}{tx,'Cost'}, results.Fastconfig2{i}{tx,'Cost'},...
    results.Dimconfig{i}{tx,'Cost'}, results.Dimconfig2{i}{tx,'Cost'}];
reconfig_type = [results.Fastconfig{i}{tx,'ReconfigType'}, results.Fastconfig2{i}{tx,'ReconfigType'},...
    results.Dimconfig{i}{tx,'ReconfigType'}, results.Dimconfig2{i}{tx,'ReconfigType'}];
cost = cost/slice_opt.Flow.ArrivalRate;
for k = 1:size(cost,2)
    idx = reconfig_type(:,k)==ReconfigType.Dimensioning;
    cost(idx,k) = cost(idx,k)*slice_opt.EventInterval;
end
hl = plot(t, cumsum(cost));
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).LineStyle = line_style{k};
    hl(k).Marker = marker{k};
    hl(k).MarkerIndices = marker_index;
end
ylabel('Cummulated Cost');
xlabel('Time');
xlim([t(1),t(end)]);
legend(legend_label(1:4), 'Location', 'northwest');

%% Profit
if exist('fig_profit_reconfig', 'var') && fig_profit_reconfig.isvalid
    figure(fig_profit_reconfig);
else
    fig_profit_reconfig = figure('Name', 'Profit with Reconfiguration');
end
% fig_profit_reconfig.OuterPosition = [100 400 400 380];
cost = [results.Fastconfig{i}{tx,'Cost'}, results.Fastconfig2{i}{tx,'Cost'},...
    results.Dimconfig{i}{tx,'Cost'}, results.Dimconfig2{i}{tx,'Cost'},...
    results.Reconfig{tx,'Cost'}*thetas(i),results.Dimconfig0{i}{tx,'Cost'}];
reconfig_type = [results.Fastconfig{i}{tx,'ReconfigType'}, results.Fastconfig2{i}{tx,'ReconfigType'},...
    results.Dimconfig{i}{tx,'ReconfigType'}, results.Dimconfig2{i}{tx,'ReconfigType'},...
    results.Reconfig{tx,'ReconfigType'}, results.Dimconfig0{i}{tx,'ReconfigType'}];
cost = cost/slice_opt.Flow.ArrivalRate;
for k = 1:size(cost,2)
    idx = reconfig_type(:,k)==ReconfigType.Dimensioning;
    cost(idx,k) = cost(idx,k)*slice_opt.EventInterval;
end
profit = [results.Fastconfig{i}{tx,'Profit'}+results.Fastconfig{i}{tx,'Cost'},...
    results.Fastconfig2{i}{tx,'Profit'}+results.Fastconfig2{i}{tx,'Cost'},...
    results.Dimconfig{i}{tx,'Profit'}+results.Dimconfig{i}{tx,'Cost'},...
    results.Dimconfig2{i}{tx,'Profit'}+results.Dimconfig2{i}{tx,'Cost'},...
    results.Reconfig{tx,'Profit'}+results.Reconfig{tx,'Cost'},...
    results.Dimconfig0{i}{tx,'Profit'}+results.Dimconfig0{i}{tx,'Cost'}];   % recover profit without reconfiguration cost
t_diff = [0;diff(results.Reconfig{tx,'Time'})];
cum_profit = cumsum(profit.*t_diff) - cumsum(cost);
hl = plot(t, cum_profit);
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).LineStyle = line_style{k};
    hl(k).Marker = marker{k};
    hl(k).MarkerIndices = marker_index;
end
ylabel('Cummulated Profit');
xlabel('Time');
xlim([t(50), t(end)]);
legend(legend_label([1:5, 7]), 'Location', 'northwest');
