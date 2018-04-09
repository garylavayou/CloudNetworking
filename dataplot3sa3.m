%% FSR and Baseline
% Single experiment.
%{
description = sprintf('%s\n%s\n%s\n%s\n%s',...
    sprintf('Experiment 502%d: Fast slice reconfiguration scheme.', type.Permanent),...
    'Topology=Sample-2.',...
    'Experiment without warm-up phase.',...
    sprintf('Slice Type %d (disable ad-hoc mode, enable dimension-trigger).', type.Index(type.Permanent))...
    );
output_name = sprintf('Results/singles/EXP502%d_fast.mat', type.Permanent);
save(output_name, 'description', 'results', 'NUM_EVENT', 'etas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}
%% declaration
line_width = [1 1 1 1 1 1];
marker_size = [6 7 6 8 6 6 6];
marker = {'none', 'none', '+', 'x', 's', 'd', 'o'};
line_style = {'-.', '-', 'none', 'none', '--', '--', '--'};
color_set = [Color.Red; Color.MildGreen; Color.MildBlue; Color.Purple; Color.Black; Color.Gray];
legend_label = {'FSR', 'Baseline'};
%%
load('Results/EXP5024_varweight_fast.mat')
ex_id = 3;
tx = 1:NUM_EVENT;
t = results.DimBaseline{1}.Time(tx);
if length(tx)>=10
    marker_index = round(linspace(1,length(tx), 10));
end
%% Number of Reconfiguration Variables
if exist('fig_num_reconfig', 'var') && fig_num_reconfig.isvalid
    figure(fig_num_reconfig);
    fig_num_reconfig.Children.delete;
else
    fig_num_reconfig = figure('Name', 'Number of Reconfiguration');
end
out_pos = fig_num_reconfig.OuterPosition;
%     out_pos(3:4) = [496 476];
out_pos(3:4) = [360 380];
fig_num_reconfig.OuterPosition = out_pos;
hl = plot(t, cumsum(results.Fastconfig{ex_id}{tx,'NumberReconfigVariables'}),'-',...
    t, cumsum(results.DimBaseline{ex_id}{tx,'NumberReconfigVariables'}), '--s');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).MarkerIndices = marker_index;
    hl(k).MarkerSize = marker_size(k);
end
ax = hl(1).Parent;
ylabel('# of Reconfigured Variables');
xlabel('Time');
% ax.OuterPosition = [-0.045 -0.01 1.14 1.04];
ax.OuterPosition = [0 0 1.08 1.01];
axis tight;
legend(legend_label, 'Location', 'northwest');
% axes('OuterPosition', [0.58, 0.20, 0.41, 0.39]);
export_fig(fig_num_reconfig, 'Figures/numac-reconfig-vartime-fast', '-pdf', '-transparent');
%% Number of Reconfigured Flows
if exist('fig_num_reconfig_flow', 'var') && fig_num_reconfig_flow.isvalid
    figure(fig_num_reconfig_flow);
    fig_num_reconfig_flow.Children.delete;
else
    fig_num_reconfig_flow = figure('Name', 'Number of Reconfigured Flows');
end
if exist('textBox', 'var')
    textBox.delete;
end
out_pos = fig_num_reconfig_flow.OuterPosition;
%     out_pos(3:4) = [496 476];
out_pos(3:4) = [360 380];
fig_num_reconfig_flow.OuterPosition = out_pos;
% fig_num_reconfig.OuterPosition = [100 400 400 380];
hl = plot(t, cumsum(results.Fastconfig{ex_id}{tx,'NumberReconfigFlows'}),'-',...
    t, cumsum(results.DimBaseline{ex_id}{tx,'NumberReconfigFlows'}), '--s');
for k=1:length(hl)
    hl(k).Color = color_set(k).RGB;
    hl(k).LineWidth = line_width(k);
    hl(k).MarkerIndices = marker_index;
    hl(k).MarkerSize = marker_size(k);
end
legend(legend_label, 'Location', 'northwest');
ylabel('# of Reconfigured Flows');
xlabel('Time');
ax = hl(1).Parent;
% ax.OuterPosition = [-0.045 -0.01 1.14 1.04];
ax.OuterPosition = [-0.02 0 1.07 1];
axis tight
scale = 10^3; scale_string = sprintf('\\times10^{%d}',log10(scale));
ax.YTickLabel = num2str(str2double(ax.YTickLabel)/scale);
textBox = text(0, 1.04, scale_string, 'Interpreter', 'tex', 'Units', 'normalized', ...
    'FontSize', 9); 
export_fig(fig_num_reconfig_flow, 'Figures/numac-reconfigflow-vartime-fast', '-pdf', '-transparent');
% posAxes = get(gca,'position');
% textBox = annotation('textbox','linestyle','none','string',['x10^{' sprintf('%d',log10(1./scale)) '}']);
% posAn = get(textBox,'position');
% set(textBox,'position',[posAxes(1) posAxes(2)+posAxes(4)-0.02 posAn(3) posAn(4)],'VerticalAlignment','cap');
%% Profit with Reconfiguration
if exist('fig_profit_reconfig', 'var') && fig_profit_reconfig.isvalid
    figure(fig_profit_reconfig);
    fig_profit_reconfig.Children.delete;
else
    fig_profit_reconfig = figure('Name', 'Profit with Reconfiguration');
end
out_pos = fig_profit_reconfig.OuterPosition;
% out_pos(3:4) = [496 476];
out_pos(3:4) = [360 380];
fig_profit_reconfig.OuterPosition = out_pos;
cost = [results.Fastconfig{ex_id}{tx,'Cost'}.*results.Fastconfig{ex_id}{tx,'Interval'},...
    results.DimBaseline{1}{tx,'Cost'}.*results.DimBaseline{1}{tx,'Interval'}];
profit = [results.Fastconfig{ex_id}{tx,'Profit'}+results.Fastconfig{ex_id}{tx,'Cost'},...
    results.DimBaseline{ex_id}{tx,'Profit'}+results.DimBaseline{ex_id}{tx,'Cost'}];   % recover profit without reconfiguration cost
%     results.Dimconfig0{i}{tx,'Profit'}+results.Dimconfig0{i}{tx,'Cost'}
t_diff = diff(results.DimBaseline{ex_id}{tx,'Time'});
cum_profit = cumsum(profit(1:(end-1),:).*t_diff) - cumsum(cost(1:(end-1),:),1);
hl = plot(t(1:end-1), cum_profit);
mi = marker_index;
mi(end) = mi(end) - 1;
for k=1:length(hl)
    kr = length(hl) - k + 1;
    hl(k).Color = color_set(kr).RGB;
    hl(k).LineWidth = line_width(kr);
    hl(k).LineStyle = line_style{kr};
    hl(k).Marker = marker{kr};
    hl(k).MarkerIndices = mi;
end
ylabel('Cummulated Profit');
xlabel('Time');
xlim([t(1), t(end-1)]);
legend(legend_label, 'Location', 'northwest');
ax = hl(1).Parent;
% ax.OuterPosition = [-0.035 -0.01 1.13 1.04];
ax.OuterPosition = [0 0 1.08 1.01];
% axes('OuterPosition', [0.55, 0.125, 0.43, 0.4]);
% hl = plot(t(1:end-1), cum_profit);
% mi = round(linspace(1,length(tx), 40));
% for k=1:length(hl)
%     kr = length(hl) - k + 1;
%     hl(k).Color = color_set(kr).RGB;
%     hl(k).LineWidth = line_width(kr);
%     hl(k).Marker = marker{kr};
%     hl(k).MarkerIndices = mi;
%     hl(k).LineStyle = line_style{kr};
% end
% axis tight
% xlim([300, 360]);
export_fig(fig_profit_reconfig, 'Figures/profit-vartime-fast', '-pdf', '-transparent');
%% Cost
if exist('fig_cost', 'var') && fig_cost.isvalid
    figure(fig_cost);
    fig_cost.Children.delete;
else
    fig_cost = figure('Name', 'Reconfiguration Cost');
end
out_pos = fig_cost.OuterPosition;
% out_pos(3:4) = [496 476];
out_pos(3:4) = [360 380];
fig_cost.OuterPosition = out_pos;
hl = plot(t, cumsum(cost(:,1)),'-',...
    t, cumsum(cost(:,2)), '--s');
for k=1:length(hl)
    kr = length(hl) - k + 1;
    hl(k).Color = color_set(kr).RGB;
    hl(k).LineWidth = line_width(kr);
    hl(k).MarkerIndices = marker_index;
    hl(k).MarkerSize = marker_size(kr);
end
ylabel('Reconfiguration Cost');
xlabel('Time');
xlim([t(1), t(end)]);
legend(legend_label, 'Location', 'northwest');
ax = hl(1).Parent;
% ax.OuterPosition = [-0.04 -0.01 1.13 1.06];
ax.OuterPosition = [-0.02 0 1.09 1.01];
scale = 10^3; scale_string = sprintf('\\times10^{%d}',log10(scale));
ax.YTickLabel = num2str(str2double(ax.YTickLabel)/scale);
textBox = text(0, 1.04, scale_string, 'Interpreter', 'tex', 'Units', 'normalized', ...
    'FontSize', 9); 
export_fig(fig_cost, 'Figures/cost-vartime-fast', '-pdf', '-transparent');
