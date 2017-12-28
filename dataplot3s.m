%%
line_width = [1 1 1 1 0.5];
color_set = [Color.Black; Color.MildGreen; Color.Purple; Color.MildBlue; Color.Red];
legend_label = {'Benchmark', 'RFFV', 'RFEV', 'Hybrid-1', 'Hybrid-2'};
i = 1;
tx = 1:NUM_EVENT;

%% Number of Reconfiguration
if exist('fig_num_reconfig', 'var') && fig_num_reconfig.isvalid
    figure(fig_num_reconfig);
else
    fig_num_reconfig = figure('Name', 'Number of Reconfiguration');
end
xlabel('Experiment Time');
% fig_num_reconfig.OuterPosition = [100 400 400 380];
yyaxis('left');
hl = plot(tx, results.Reconfig{tx,'NumberReconfigVariables'}, '-.', ...
    tx, results.Fastconfig{i}{tx,'NumberReconfigVariables'}, ':',...
    tx, results.Fastconfig2{i}{tx,'NumberReconfigVariables'}, '.',...
    tx, results.Dimconfig{i}{tx,'NumberReconfigVariables'},'-',...
    tx, results.Dimconfig2{i}{tx,'NumberReconfigVariables'},'+');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
end
ylabel('Number of Reconfiguration');
ylimmax = hl(1).Parent.YLim(2);
yyaxis right
hr = plot(tx, results.Reconfig{tx,'NumberVariables'}, '--');
hr(1).Color = Color.Red.RGB;
hr(1).LineWidth = 1;
ylabel('Number of Variables');
hr.Parent.YColor = Color.Red.RGB;
% title(sprintf('\\eta=%.2f', etas(i)));
legend([legend_label, {'Variable#'}], 'Location', 'northwest');
ylimmax = max(ylimmax, hr.Parent.YLim(2));
yyaxis left;
ylim(hl(1).Parent, [0, ylimmax]);
yyaxis right;
ylim(hr.Parent, [0, ylimmax]);

%% Ratio of Reconfiguration
if exist('fig_ratio_reconfig', 'var') && fig_ratio_reconfig.isvalid
    figure(fig_ratio_reconfig);
else
    fig_ratio_reconfig = figure('Name', 'Ratio of Reconfiguration');
end
% fig_num_reconfig.OuterPosition = [100 400 400 380];
% yyaxis left;
ratio_reconfig = results.Reconfig{tx,'NumberReconfigVariables'}./...
    results.Reconfig{tx,'NumberVariables'};
ratio_fastconfig = results.Fastconfig{i}{tx,'NumberReconfigVariables'}./...
    results.Fastconfig{i}{tx,'NumberVariables'};
ratio_fastconfig2 = results.Fastconfig2{i}{tx,'NumberReconfigVariables'}./...
    results.Fastconfig2{i}{tx,'NumberVariables'};
ratio_dimconfig = results.Dimconfig{i}{tx,'NumberReconfigVariables'}./...
    results.Dimconfig{i}{tx,'NumberVariables'};
ratio_dimconfig2 = results.Dimconfig2{i}{tx,'NumberReconfigVariables'}./...
    results.Dimconfig2{i}{tx,'NumberVariables'};
hl = plot(tx, ratio_reconfig, '-.', ...
    tx, ratio_fastconfig, ':',...
    tx, ratio_fastconfig2, '.',...
    tx, ratio_dimconfig, '-',...
    tx, ratio_dimconfig2, '+');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
end
ylabel('Ratio of Reconfiguration');
xlabel('Experiment Time');
% yyaxis right
% hr = plot(tx, g_results.slice.num_flows(tx), 'r:',...
%     tx, num_nonzeros_reconfig(tx), '-.',...
%     tx, num_nonzeros_fastconfig(tx), '--');
% hr(2).Color = [0.871, 0.49, 0];
% hr(3).Color = [0.749, 0, 0.749];
% ylabel('Number of Reconfiguration');
% legend({'Ratio-Reconfig', 'Ratio-FastConfig', 'Flows', 'Reconfig', 'Fast Reconfig'}, ...
%     'Location', 'northwest');
ylim([0,1]);
legend(legend_label, 'Location', 'best', 'Orientation', 'vertical');
% title(sprintf('\\eta=%.2f', etas(i)));

%% Cost of Reconfiguration
if exist('fig_cost_reconfig', 'var') && fig_cost_reconfig.isvalid
    figure(fig_cost_reconfig);
else
    fig_cost_reconfig = figure('Name', 'Cost of Reconfiguration');
end
hl = plot(tx, results.Reconfig{tx,'Cost'}*etas(i), '-.', ...
    tx, results.Fastconfig{i}{tx,'Cost'}, ':', ...
    tx, results.Fastconfig2{i}{tx,'Cost'}, '.', ...
    tx, results.Dimconfig{i}{tx,'Cost'}, '-', ...
    tx, results.Dimconfig2{i}{tx,'Cost'}, '+');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
end
ylabel('Cost');
xlabel('Experiment Time');
% hold on;
% plot(tx, results.Reconfig{tx,'ResourceCost'}, '--');
% hold off;
% hr.Color = Color.Orange.RGB;
% hr.Parent.YColor = Color.Orange.RGB;
legend(legend_label, 'Location', 'northwest');
% title(sprintf('\\eta=%.2f', etas(i)));

%% Number of Reconfigured Flows
if exist('fig_flow_reconfig', 'var') && fig_flow_reconfig.isvalid
    figure(fig_flow_reconfig);
else
    fig_flow_reconfig = figure('Name', 'Number of Reconfigured Flows');
end
hl = plot(tx, results.Reconfig{tx,'NumberReconfigFlows'}, '-.', ...
    tx, results.Fastconfig{i}{tx,'NumberReconfigFlows'}, ':', ...
    tx, results.Fastconfig2{i}{tx,'NumberReconfigFlows'}, '.', ...
    tx, results.Dimconfig{i}{tx,'NumberReconfigFlows'}, '-', ...
    tx, results.Dimconfig2{i}{tx,'NumberReconfigFlows'}, '+');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
end
ylabel('Number of Reconfigured Flows');
ylimmax(1) = hl(1).Parent.YLim(2);
yyaxis right;
hr = plot(tx, results.Reconfig{tx,'NumberFlows'});
hr(1).Color = Color.Red.RGB;
hr(1).LineWidth = 1;
ylimmax(2) = hr(1).Parent.YLim(2);
ylim([0, max(ylimmax)+2])
ylabel('Number of Flows');
hr.Parent.YColor = Color.Red.RGB;
yyaxis left;
ylim([0, max(ylimmax)+2]);
xlabel('Experiment Time');
legend([legend_label,{'Flows #'}], 'Location', 'northwest');
% title(sprintf('\\eta=%.2f', etas(i)));

%% Profit
if exist('fig_profit_reconfig', 'var') && fig_profit_reconfig.isvalid
    figure(fig_profit_reconfig);
else
    fig_profit_reconfig = figure('Name', 'Profit with Reconfiguration');
end
% fig_profit_reconfig.OuterPosition = [100 400 400 380];
xlabel('Experiment Time');
yyaxis left;
profit = results.Reconfig{tx,'Profit'} + (1-etas(i))*results.Reconfig{tx,'Cost'};
hl = plot(tx, profit, '-.', ...
    tx, results.Fastconfig{i}{tx,'Profit'}, ':',...
    tx, results.Fastconfig2{i}{tx,'Profit'}, '.',...
    tx, results.Dimconfig{i}{tx,'Profit'}, '-',...
    tx, results.Dimconfig2{i}{tx,'Profit'}, '+');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
end
%     hl(1).Parent.YLim(1) = 0;
ylabel('Profit');
yyaxis right;
hr = plot(tx, results.Reconfig{tx,'NumberFlows'}, '--');
hr.Color = Color.Orange.RGB;
hr.Parent.YLim(1) = 0;
hr.Parent.YLim(2) = hr.Parent.YLim(2)+2;
ylabel('Number of Flows');
legend([legend_label, {'Flow #'}], 'Location', 'northwest');
% title(sprintf('\\eta=%.2f', etas(i)));
