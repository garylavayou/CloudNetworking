%%
legend_label_full = {'Price-SPP', 'Dynamic-Slicing','Dynamic-Slicing2',...
    'Static-Slicing', 'Upperbound-SPP'};
legend_label_compact = {'PS', 'DS', 'SS'};
line_width = 1.5;
font_name = 'Calibri';
font_size = 12;
x = slice_config.Count;
num_config = length(slice_config.RatioSet);

%% net social welfare
if exist('fig_net_social_welfare', 'var') && fig_net_social_welfare.isvalid
    figure(fig_net_social_welfare);
else
    fig_net_social_welfare = figure('Name', 'Net Social Welfare');
end
switch num_config
    case 1
        fig_net_social_welfare.OuterPosition = [100 400 400 380];
    case 4
        fig_net_social_welfare.OuterPosition = [100 400 1600 380];
end

for i = 1:num_config
    subplot(1,4,i);
%         x, output_results(i).Price1.Stat.AccurateWelfare, '--^',...
    h = plot(x, output_results(i).Optimal.Stat.AccurateWelfare(:,2), 'o',...
        x, output_results(i).Price2.Stat.AccurateWelfare, '--+',...
        x, output_results(i).Static.Stat.AccurateWelfare, '-.x',...
        x, output_results(i).Optimal.Stat.AccurateWelfare(:,1), '--');
%     legend(legend_label_full, 'Location', 'northwest');
    legend(legend_label_full([1 2 4 5]), 'Location', 'northwest');
    xlabel('Total number of slices');
    ylabel('Net social welfare');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    ylim(welfare_limit);
    xlim([slice_config.Count(1), slice_config.Count(end)]);
end

%% node utilization
if exist('fig_node_utilization', 'var') && fig_node_utilization.isvalid
    figure(fig_node_utilization);
    clf(fig_node_utilization)
else
    fig_node_utilization = figure('Name', 'Node Utilization');
end
switch num_config
    case 1
        fig_node_utilization.OuterPosition = [100 400 300 380];
    case 4
        fig_node_utilization.OuterPosition = [100 400 1600 380];
end

for i = 1:num_config
    subplot(1,4,i);
%         x, output_results(i).Price1.Stat.NodeUtilization(:,1), '--^',...
    h = plot(x, output_results(i).Optimal.Stat.NodeUtilization(:,1), '-o',...
        x, output_results(i).Price2.Stat.NodeUtilization(:,1), '--+',...
        x, output_results(i).Static.Stat.NodeUtilization(:,1), '-.x');
%     hold on;
%     h(1) = errorbar(x, u_node_optimal.Average, u_node_optimal.StdDev/4, '-');
%     h(2) = errorbar(x, u_node_price.Average, u_node_price.StdDev/4, '--');
%     h(3) = errorbar(x, u_node_static.Average, u_node_static.StdDev/4, '-.');
%     legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
    legend(legend_label_full([1 2 4]), 'Location', 'southeast');
    xlabel('Total number of slices');
    ylabel('Node Utilization');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    ylim([0, 1]);
    xlim([slice_config.Count(1), slice_config.Count(end)]);
end

%% link utilization
% Since the links are bi-directional, some links might not be utilized. Therefore, the
% link utilization cannot approach 1.
if exist('fig_link_utilization', 'var') && fig_link_utilization.isvalid
    figure(fig_link_utilization);
else
    fig_link_utilization = figure('Name', 'Link Utilization');
end
switch num_config
    case 1
        fig_link_utilization.OuterPosition = [100 400 400 380];
    case 4
        fig_link_utilization.OuterPosition = [100 400 1600 380];
end

for i = 1:num_config
    subplot(1,4,i);
%         x, output_results(i).Price1.Stat.LinkUtilization(:,1), '--^',...
    h = plot(x, output_results(i).Optimal.Stat.LinkUtilization(:,1), '-o',...
        x, output_results(i).Price2.Stat.LinkUtilization(:,1), '--+',...
        x, output_results(i).Static.Stat.LinkUtilization(:,1), '-.x');
%     legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
    legend(legend_label_full([1 2 4]), 'Location', 'northwest');
    xlabel('Total number of slices');
    ylabel('Link Utilization');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    axis([slice_config.Count(1), slice_config.Count(end), 0, 1]);
end
%% Total network operational cost
if exist('fig_network_cost', 'var') && fig_network_cost.isvalid
    figure(fig_network_cost);
else
    fig_network_cost = figure('Name', 'Network Operation Cost');
end
switch num_config
    case 1
        fig_network_cost.OuterPosition = [100 400 400 380];
    case 4
        fig_network_cost.OuterPosition = [100 400 1600 380];
end

for i = 1:num_config
    subplot(1,4,i);
%         x, output_results(i).Price1.Stat.Cost(:,1), '--^',...
    h = plot(x, output_results(i).Optimal.Stat.Cost(:,1), '-o',...
        x, output_results(i).Price2.Stat.Cost(:,1), '--+',...
        x, output_results(i).Static.Stat.Cost(:,1), '-.x');
%     legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
    legend(legend_label_full([1 2 4]), 'Location', 'southeast');
    xlabel('Total number of slices');
    ylabel('Network Operation Cost');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    xlim([slice_config.Count(1), slice_config.Count(end)]);
end

%% SP's profit
if exist('fig_sp_profit', 'var') && fig_sp_profit.isvalid
    figure(fig_sp_profit);
else
    fig_sp_profit = figure('Name', 'SP Profit');
end
switch num_config
    case 1
        fig_sp_profit.OuterPosition = [100 400 400 380];
    case 4
        fig_sp_profit.OuterPosition = [100 400 1600 380];
end

for i = 1:num_config
    subplot(1,num_config,i);
%         x, output_results(i).Price1.Stat.Profit(:,1), '--^',...
    h = plot(x, output_results(i).Optimal.Stat.Profit(:,1), '-o',...
        x, output_results(i).Price2.Stat.Profit(:,1), '--+',...
        x, output_results(i).Static.Stat.Profit(:,1), '-.x');
%     legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
    legend(legend_label_full([1 2 4]), 'Location', 'northwest');
    xlabel('Total number of slices');
    ylabel('Profit of Slice Provider');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    axis([slice_config.Count(1), slice_config.Count(end), profit_limit]);
%     oldLabels = str2double(get(gca,'YTickLabel'));
%     scale = 10^-3;newLabels = num2str(oldLabels*scale);
%     set(gca,'YTickLabel',newLabels,'units','normalized');
%     posAxes = get(gca,'position');
%     textBox = annotation('textbox','linestyle','none','string',['x 10\it^{' sprintf('%d',log10(1./scale)) '}']);
%     posAn = get(textBox,'position');
%     set(textBox,'position',[posAxes(1) posAxes(2)+posAxes(4) posAn(3) posAn(4)],'VerticalAlignment','cap');

end

%% SC's profit
% if exist('fig_sc_profit', 'var') && fig_sc_profit.isvalid
%     figure(fig_sc_profit);
% else
%     fig_sc_profit = figure('Name', 'SC profit');
% end
