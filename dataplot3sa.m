%% Single plot of accumulative quantity
% The independent variables is time.
line_width = [1 1 1 1 1 1];
marker_size = [6 7 6 8 6 6 6];
marker = {'none', '+', 'none', 'x', 's', 'd', 'o'};
line_style = {'-.', 'none', '-', 'none', '--', '--', '--'};
color_set = [Color.MildGreen; Color.Purple; Color.MildBlue; Color.Red; Color.Black; Color.Gray];
legend_label = {'FSR0', 'FSR', 'HSR0', 'HSR', 'Baseline', '#Variables'};
i = 4;
if strcmpi(mode, 'var-eta')
    j = 1;
else
    j = i;
end
tx = 50:NUM_EVENT;
t = results.DimBaseline{j}.Time(tx);
if length(tx)>=10
    marker_index = round(linspace(1,length(tx), 10));
end
%% Number of Reconfiguration (with benchmark)
% TODO: exclude those with zero reconfigureation cost (linear).
if exist('fig_num_reconfig', 'var') && fig_num_reconfig.isvalid
    figure(fig_num_reconfig);
else
    fig_num_reconfig = figure('Name', 'Number of Reconfiguration (With benchmark)');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
yyaxis('left');
hl = plot(t, cumsum(results.Fastconfig{i}{tx,'NumberReconfigVariables'}), '-.',...
    t, cumsum(results.FastconfigReserve{i}{tx,'NumberReconfigVariables'}), '+',...
    t, cumsum(results.Dimconfig{i}{tx,'NumberReconfigVariables'}),'-',...
    t, cumsum(results.DimconfigReserve{i}{tx,'NumberReconfigVariables'}),'x',...
    t, cumsum(results.DimBaseline{j}{tx,'NumberReconfigVariables'}), '--s');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).MarkerIndices = marker_index;
    hl(k).MarkerSize = marker_size(k);
end
ylabel('Number of Reconfiguration');
ylimmax = hl(1).Parent.YLim(2);
yyaxis right
hr = plot(t, cumsum(results.DimBaseline{j}{tx,'NumberVariables'}), '--d');
hr(1).Color = Color.Gray.RGB;
hr(1).LineWidth = 1;
hr(1).MarkerIndices = marker_index;
ylabel('Number of Variables');
hr.Parent.YColor = Color.Gray.RGB;
% title(sprintf('\\eta=%.2f', etas(i)));
legend(legend_label(1:6), 'Location', 'northwest');
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
    t, cumsum(results.FastconfigReserve{i}{tx,'NumberReconfigVariables'}), '--+',...
    t, cumsum(results.Dimconfig{i}{tx,'NumberReconfigVariables'}),'-',...
    t, cumsum(results.DimconfigReserve{i}{tx,'NumberReconfigVariables'}),'--x');
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
% If dimensioning is triggered, the reconfiguration cost has been
% distibuted to the whole dimensioning interval. To calculate the
% reconfiguration cost conveniently, we recollect the cost to a single
% configure interval.
if exist('fig_cost_reconfig', 'var') && fig_cost_reconfig.isvalid
    figure(fig_cost_reconfig);
else
    fig_cost_reconfig = figure('Name', 'Cost of Reconfiguration (With benchmark)');
end
cost = [results.Fastconfig{i}{tx,'Cost'}.*results.Fastconfig{i}{tx,'Interval'}, ...
    results.FastconfigReserve{i}{tx,'Cost'}.*results.FastconfigReserve{i}{tx,'Interval'},...
    results.Dimconfig{i}{tx,'Cost'}.*results.Dimconfig{i}{tx,'Interval'}, ...
    results.DimconfigReserve{i}{tx,'Cost'}.*results.DimconfigReserve{i}{tx,'Interval'}];
if strcmpi(mode, 'var-eta')
    cost = [cost, results.DimBaseline{1}{tx,'Cost'}*(etas(i)/etas(1))...
        .*results.DimBaseline{1}{tx,'Interval'}];
else
    cost = [cost, results.DimBaseline{i}{tx,'Cost'}.*results.DimBaseline{i}{tx,'Interval'}];
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

%% Profit
% Reconfiguration cost is counted as a instantanuous value, while the profit excluding
% reconfiguration cost is distrbuted in the whole configuration interval.
if exist('fig_profit_reconfig', 'var') && fig_profit_reconfig.isvalid
    figure(fig_profit_reconfig);
else
    fig_profit_reconfig = figure('Name', 'Profit with Reconfiguration');
end
% fig_profit_reconfig.OuterPosition = [100 400 400 380];
profit = [results.Fastconfig{i}{tx,'Profit'}+results.Fastconfig{i}{tx,'Cost'},...
    results.FastconfigReserve{i}{tx,'Profit'}+results.FastconfigReserve{i}{tx,'Cost'},...
    results.Dimconfig{i}{tx,'Profit'}+results.Dimconfig{i}{tx,'Cost'},...
    results.DimconfigReserve{i}{tx,'Profit'}+results.DimconfigReserve{i}{tx,'Cost'},...
    results.DimBaseline{j}{tx,'Profit'}+results.DimBaseline{j}{tx,'Cost'}];   % recover profit without reconfiguration cost
%     results.Dimconfig0{i}{tx,'Profit'}+results.Dimconfig0{i}{tx,'Cost'}
t_diff = diff(results.DimBaseline{j}{tx,'Time'});
cum_profit = cumsum(profit(end-1,:).*t_diff) - cumsum(cost(end-1,:),1);
hl = plot(t(1:end-1), cum_profit);
mi = marker_index;
mi(end) = mi(end) - 1;
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).LineStyle = line_style{k};
    hl(k).Marker = marker{k};
    hl(k).MarkerIndices = mi;
end
ylabel('Cummulated Profit');
xlabel('Time');
xlim([t(1), t(end-1)]);
legend(legend_label(1:5), 'Location', 'northwest');

%% Cost of Reconfiguration (no benchmark)
if exist('fig_cost_reconfig_nb', 'var') && fig_cost_reconfig_nb.isvalid
    figure(fig_cost_reconfig_nb);
else
    fig_cost_reconfig_nb = figure('Name', 'Cost of Reconfiguration (No benchmark)');
end
hl = plot(t, cumsum(cost(:,1:4)));  % the first four methods' cost
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

%% Number of Flow Reconfiguration (with benchmark)
if exist('fig_num_reconfig_flow', 'var') && fig_num_reconfig_flow.isvalid
    figure(fig_num_reconfig_flow);
else
    fig_num_reconfig_flow = figure('Name', 'Number of Reconfigured Flows (With benchmark)');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
yyaxis('left');
hl = plot(t, cumsum(results.Fastconfig{i}{tx,'NumberReconfigFlows'}), '-.',...
    t, cumsum(results.FastconfigReserve{i}{tx,'NumberReconfigFlows'}), '+',...
    t, cumsum(results.Dimconfig{i}{tx,'NumberReconfigFlows'}),'-',...
    t, cumsum(results.DimconfigReserve{i}{tx,'NumberReconfigFlows'}),'x',...
    t, cumsum(results.DimBaseline{j}{tx,'NumberReconfigFlows'}), '--s');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).MarkerIndices = marker_index;
    hl(k).MarkerSize = marker_size(k);
end
ylabel('Number of Reconfiguration');
ylimmax = hl(1).Parent.YLim(2);
yyaxis right
hr = plot(t, cumsum(results.DimBaseline{j}{tx,'NumberFlows'}), '--d');
hr(1).Color = Color.Gray.RGB;
hr(1).LineWidth = 1;
hr(1).MarkerIndices = marker_index;
ylabel('Number of Flows');
hr.Parent.YColor = Color.Gray.RGB;
% title(sprintf('\\eta=%.2f', etas(i)));
legend(legend_label(1:6), 'Location', 'northwest');
ylimmax = max(ylimmax, hr.Parent.YLim(2));
yyaxis left;
ylim(hl(1).Parent, [0, ylimmax]);
yyaxis right;
ylim(hr.Parent, [0, ylimmax]);
xlabel('Time');
xlim([t(1), t(end)]);
%% Number of Flow Reconfiguration (No benchmark)
if exist('fig_num_reconfig_flow_nb', 'var') && fig_num_reconfig_flow_nb.isvalid
    figure(fig_num_reconfig_flow_nb);
else
    fig_num_reconfig_flow_nb = figure('Name', 'Number of Reconfigured Flows (No benchmark)');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
hl = plot(t, cumsum(results.Fastconfig{i}{tx,'NumberReconfigFlows'}), '-.',...
    t, cumsum(results.FastconfigReserve{i}{tx,'NumberReconfigFlows'}), '--+',...
    t, cumsum(results.Dimconfig{i}{tx,'NumberReconfigFlows'}),'-',...
    t, cumsum(results.DimconfigReserve{i}{tx,'NumberReconfigFlows'}),'--x');
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