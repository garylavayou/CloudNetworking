%%
legend_label_full = {'Optimal-SPP', 'Price-SPP', 'Dynamic-Slicing','Dynamic-Slicing2', ...
    'Static-Slicing'};
legend_label_compact = {'OS', 'DS', 'SS'};
line_width = 1.5;
font_name = 'Calibri';
font_size = 12;
x = slice_config.Count;

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
    subplot(1,num_config,i);
    h = plot(x, output_results{i}{1}.optimal(:,1), '-',...
        x, output_results{i}{1}.optimal(:,2), 'o',...
        x, output_results{i}{1}.price1, '--^',...
        x, output_results{i}{1}.price2, '--+',...
        x,output_results{i}{1}.static, '-.x');
    legend(legend_label_full, 'Location', 'northwest');
    xlabel('Total number of slices');
    ylabel('Net social welfare');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    ylim([0, h(1).Parent.YLim(2)]*1.2);
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
        fig_node_utilization.OuterPosition = [100 400 400 380];
    case 4
        fig_node_utilization.OuterPosition = [100 400 1600 380];
end

for i = 1:num_config
    subplot(1,num_config,i);
    h = plot(x, output_results{i}{9}.optimal.Average, '-o',...
        x, output_results{i}{9}.price1.Average, '--^',...
        x, output_results{i}{9}.price2.Average, '--+',...
        x, output_results{i}{9}.static.Average, '-.x');
%     hold on;
%     h(1) = errorbar(x, u_node_optimal.Average, u_node_optimal.StdDev/4, '-');
%     h(2) = errorbar(x, u_node_price.Average, u_node_price.StdDev/4, '--');
%     h(3) = errorbar(x, u_node_static.Average, u_node_static.StdDev/4, '-.');
    legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
    xlabel('Total number of slices');
    ylabel('Node Utilization');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    axis([slice_config.Count(1), slice_config.Count(end), 0, 1.2]);
end

%% link utilization
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
    subplot(1,num_config,i);
    % , u_link_optimal.StdDev, , u_link_price.StdDev, u_link_static.StdDev, 
    h = plot(x, output_results{i}{10}.optimal.Average, '-o',...
        x, output_results{i}{10}.price1.Average, '--^',...
        x, output_results{i}{10}.price2.Average, '--+',...
        x, output_results{i}{10}.static.Average, '-.x');
    legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
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
    subplot(1,num_config,i);
    h = plot(x, output_results{i}{11}.optimal, '-o',...
        x, output_results{i}{11}.price1, '--^',...
        x, output_results{i}{11}.price2, '--+',...
        x, output_results{i}{11}.static, '-.x');
    legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
    xlabel('Total number of slices');
    ylabel('Network Operation Cost');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    axis([slice_config.Count(1), slice_config.Count(end), cost_limit]);
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
    h = plot(x, output_results{i}{3}.optimal{end}, '-o',...
        x, output_results{i}{3}.price1{end}, '--^',...
        x, output_results{i}{3}.price2{end}, '--+',...
        x,output_results{i}{3}.static{end}, '-.x');
    legend(legend_label_full([1 3 4 5]), 'Location', 'northwest');
    xlabel('Total number of slices');
    ylabel('Profit of Slice Provider');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size;
    for j = 1:length(h)
        h(j).LineWidth = line_width;
    end
    axis([slice_config.Count(1), slice_config.Count(end), profit_limit]);
end

%% SC's profit
% if exist('fig_sc_profit', 'var') && fig_sc_profit.isvalid
%     figure(fig_sc_profit);
% else
%     fig_sc_profit = figure('Name', 'SC profit');
% end
