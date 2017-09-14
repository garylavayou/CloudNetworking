%%
% if b_single_optimal
%     %     plot(x, profit_approx.optimal(x), x, profit_accurate.optimal(x));
% end
% if b_price_adjust
% %     plot(x, profit_approx.price(x), x, profit_accurate.price(x));
% end
%%
color.blue = 'blue';
color.red = 'red';
color.black = 'black';
color.green = [0 0.69 0.314];
color.midblue = [0 0.45 0.74];
color.orange = [0.85 0.33 0.1];
color.darkpurple = [0.49 0.18 0.56];
color.lightgreen = [0.47 0.67 0.19];
color.lightblue = [0.3 0.75 0.93];
color.purple = [1 0 1];
color.midgreen = [0.13 0.69 0.3];
legend_label_full = {'Price-SPP', 'Dynamic-Slicing','Static-Slicing'};
legend_label_compact = {'PS', 'DS', 'SS'};
line_width = 1;
font_name = 'Calibri';
font_size = 12;
%%
if exist('fig_number_slice', 'var') && fig_number_slice.isvalid
    figure(fig_number_slice);
else
    fig_number_slice = figure('Name', 'Number of Slices');
%     scrsz = get(groot,'ScreenSize');
%     fig_1 = figure('OuterPosition',[scrsz(1)-5 scrsz(2)+34 scrsz(3)+10 scrsz(4)-34], ...
%         'Name', 'Numerical Results of Slices and Network');
end
% offset = 0;
% x = (1+offset):(num_events-offset);
% x = 1:i-1;
%     x, sum(stat.number_reject(x,:),2), '-m',...
%     x, sum(stat.number_part_reject(x,:),2), '--b'
% , 'Reject', 'Part-Reject'
% h = plot(stat.times(x_tick), stat_optimal.NumberSlices(x_tick,1), '--', ...
%     stat.times(x_tick), stat_optimal.NumberSlices(x_tick,2), ':', ...
%     stat.times(x_tick), stat_optimal.NumberSlices(x_tick,3), '-.', ...
%     stat.times(x_tick), sum(stat_optimal.NumberSlices(x_tick,:),2), '-',...
%     stat.times(x_tick), sum(stat_static.NumberSlicesStatic(x_tick,:),2), 'x');
h = plot(stat.times(x_tick), number_slices(x_tick,1), '--', ...
    stat.times(x_tick), number_slices(x_tick,2), ':', ...
    stat.times(x_tick), number_slices(x_tick,3), '-.', ...
    stat.times(x_tick), sum(number_slices(x_tick,:),2), '-',...
    stat.times(x_tick), sum(stat_static.NumberSlicesStatic(x_tick,:),2), 'x');
h(1).Parent.FontName = font_name;
h(1).Parent.FontSize = font_size;
h(1).LineWidth = line_width;
h(2).LineWidth = line_width+1;
h(3).LineWidth = line_width;
h(4).LineWidth = line_width;
h(5).LineWidth = line_width;
h(1).Color = color.midblue;
h(2).Color = color.red;
h(3).Color = color.green;
h(4).Color = color.orange;
h(5).Color = color.blue;
legend({'Type-1', 'Type-2', 'Type-3', 'Total', 'Total Static'}, ...
    'Location', 'northwest')
xlabel('Time (hour)');
ylabel('Number of Slices');
xlim([stat.times(1), stat.times(x_tick(end))]);
%%
if exist('fig_netsocialwelfare', 'var') && fig_netsocialwelfare.isvalid
    figure(fig_netsocialwelfare);
else
    fig_netsocialwelfare = figure('Name', 'Network Social Welfare');
end
h = plot(stat.times(x_tick), stat_optimal.AccurateWelfare(x_tick,2),'--r', ...
    stat.times(x_tick), stat_price2.AccurateWelfare(x_tick), '-b', ...
    stat.times(x_tick), stat_static.AccurateWelfare(x_tick),'-.g');
%     x, profit_accurate.price2(x), '-c', ...
%     x, profit_accurate.part(x), '-g',...
%     x, profit_accurate.partprice(x),'-k'
h(1).Parent.FontName = font_name;
h(1).Parent.FontSize = font_size;
h(1).LineWidth = line_width;
h(2).LineWidth = line_width;
h(3).LineWidth = line_width;
h(3).Color = color.green;
legend(legend_label_full, 'Location', 'northwest');  % 'price2', 'partition', 
xlabel('Time (hour)');
ylabel('Net Social Welfare');
xlim([stat.times(1), stat.times(x_tick(end))]);
%%
if exist('fig_others', 'var') && fig_others.isvalid
    figure(fig_others);
else
    scrsz = get(groot,'ScreenSize');
    fig_others = figure('OuterPosition',[scrsz(1)-5 scrsz(2)+34 scrsz(3)+10 scrsz(4)-34], ...
        'Name', 'Others');
end
subplot(2,3,1);
devi1 = profit_accurate.optimal(x_tick) - profit_accurate.price(x_tick);
devi2 = profit_accurate.optimal(x_tick) - profit_accurate.static(x_tick);
% devi2 = profit_accurate.optimal(x) - profit_accurate.part(x);
% devi4 = profit_accurate.optimal(x) - profit_accurate.price2(x);
% devi5 = profit_accurate.optimal(x) - profit_accurate.partprice(x);
plot(stat.times(x_tick), devi1,'-r', stat.times(x_tick), devi2, '-b'); % , x, devi3, '-g', x, devi4, '-c', x, devi5, '-k'
clear devi1 devi2;
xlim(xlimits);
legend({'diff\_price', 'diff\_static'}, 'Location', 'northwest') % 'diff\_part', 'diff\_price2'
xlabel('Time (hour)');
ylabel('Difference of optimal objective value');
title('Distance to optimal net social welfare');
subplot(2,3,2);
plot(stat.times(x_tick), utilization.optimal(x_tick), '-r', ...
    stat.times(x_tick), utilization.price(x_tick), '-b', ...
    stat.times(x_tick), utilization.static(x_tick), '-m');
%     x, utilization.price2(x), '-c', ...
%     x, utilization.part(x), '-g', ...
%     x, utilization.partprice(x), '-k'
xlim(xlimits);
clear h lim;
legend(legend_label_full,'Location', 'southeast') % 'price2', 'partition', 
xlabel('Time (hour)');
ylabel('Network Utilization');
title('Network utilization from different methods');
subplot(2,3,3)
plot(stat.times(x_tick), stat.number_flows(x_tick,1), '-b', ...
    stat.times(x_tick), stat.number_flows_static(x_tick,1), 'x');
legend({'Dynamic', 'Static'},'Location', 'northwest')
xlabel('Time (hour)');
ylabel('Number of Flows');
title('Total number of flows in the network');
xlim(xlimits);
subplot(2,3,4)
plot(stat.times(x_tick), runtime.optimal(x_tick,1), '-b', ...
    stat.times(x_tick), runtime.price(x_tick,1), '-r',...
    stat.times(x_tick), runtime.static(x_tick), '-m');
%     x, runtime.price2(x)./sum(stat.number_slices(x,:),2), '-c',...
%     x, runtime.part(x)./sum(stat.number_slices(x,:),2), '-g',...
xlabel('Time (hour)');
ylabel('Run Time (seconds)');
title('Run time of the optimization process');
legend(legend_label_full, 'Location', 'northwest'); % 'price2', 'part', 
xlim(xlimits);
%clear x;

%% statistic of the profit of different type of slice.
% if exist('fig_2', 'var') && fig_2.isvalid
%     figure(fig_2);
% else
%     fig_2 = figure('Name', 'Comparison of Slices');
% end
% x_tick = 1;
% num_column = 3;
% Y = zeros(num_type+1,num_column);
% for t = 1:num_type
%     Y(t,:) = [profit_stat.optimal{t}.Average(x_tick) ...
%         profit_stat.price{t}.Average(x_tick) ...
%         profit_stat.static{t}.Average(x_tick)];
% %         profit_stat.price2{i}.Average(x) ...
% %         profit_stat.part{i}.Average(x) ...
% end
% % Y(end,:) = [profit_stat.optimal{end}(x), profit_stat.price{end}(x), ...
% %     profit_stat.part{end}(x), profit_stat.static{end}(x)];
% Y(end,1) = profit_approx.optimal(x_tick) ...
%     - profit_stat.optimal{1}.Average(x_tick)*stat.number_slices(x_tick,1) ...
%     - profit_stat.optimal{2}.Average(x_tick)*stat.number_slices(x_tick,2) ...
%     - profit_stat.optimal{3}.Average(x_tick)*stat.number_slices(x_tick,3);
% Y(end,2) = profit_approx.price(x_tick) ...
%     - profit_stat.price{1}.Average(x_tick)*stat.number_slices(x_tick,1) ...
%     - profit_stat.price{2}.Average(x_tick)*stat.number_slices(x_tick,2) ...
%     - profit_stat.price{3}.Average(x_tick)*stat.number_slices(x_tick,3);
% Y(end,3) = profit_approx.static(x_tick) ...
%     - profit_stat.static{1}.Average(x_tick)*stat.number_slices_static(x_tick,1) ...
%     - profit_stat.static{2}.Average(x_tick)*stat.number_slices_static(x_tick,2) ...
%     - profit_stat.static{3}.Average(x_tick)*stat.number_slices_static(x_tick,3);
% % Y(end,3) = profit_approx.price2(x) ...
% %     - profit_stat.price2{1}.Average(x)*stat.number_slices(x,1) ...
% %     - profit_stat.price2{2}.Average(x)*stat.number_slices(x,2) ...
% %     - profit_stat.price2{3}.Average(x)*stat.number_slices(x,3);
% % Y(end,4) = profit_approx.part(x) ...
% %     - profit_stat.part{1}.Average(x)*stat.number_slices(x,1) ...
% %     - profit_stat.part{2}.Average(x)*stat.number_slices(x,2) ...
% %     - profit_stat.part{3}.Average(x)*stat.number_slices(x,3);
% b = bar(Y);
% b(1).Parent.XTickLabel = {'Type-1', 'Type-2', 'Type-3', 'Network'};
% if b(1).Parent.YLim(1) < -1000
%     b(1).Parent.YLim = [-1000, b(1).Parent.YLim(2)];
% end
% legend({'optimal', 'price', 'static'}, 'Location', 'northwest'); % 'price2', 'partition', 
% xlabel('Type of slices or network');
% ylabel('Avearge profit of each slice/network');
% title(sprintf('Comparison of profit (Event ID: %d)', x_tick));

%% statistic of the profit of different type of slice.
if exist('fig_sp_profit', 'var') && fig_sp_profit.isvalid
    figure(fig_sp_profit);
else
    fig_sp_profit = figure('Name', 'Comparison the Profit of Slice Provider');
end
% limits = [0, 40, 0, 25000;
%     0, 40, 0, 2000;
%     0, 40, 1000, 6000;
%     0, 40, 3000, 11000];
h = plot(stat.times(x_tick), stat_optimal.Profit(x_tick), 'r--',...
    stat.times(x_tick), stat_price2.Profit(x_tick), 'b-',...
    stat.times(x_tick), stat_static.Profit(x_tick), 'g-.');
h(1).Parent.FontName = font_name;
h(1).Parent.FontSize = font_size;
h(1).LineWidth = line_width;
h(2).LineWidth = line_width;
h(3).LineWidth = line_width;
h(3).Color = color.green;
legend(legend_label_full, 'Location', 'northwest');
ylabel('Profit');
xlabel('Time (hour)');
xlim([stat.times(1), stat.times(x_tick(end))]);
%%
if exist('fig_sc_profit', 'var') && fig_sc_profit.isvalid
    figure(fig_sc_profit);
else
    fig_sc_profit = figure('Name', 'Comparison the Profit of Slices');
%     fig_sc_profit.Resize
end
scrsz = get(groot,'ScreenSize');  
fig_sc_profit.OuterPosition = [100 400 1200 380];
% plot the profit of each type of slice
for t = 1:length(type.Index)
    subplot(1,3,t)
    h = plot(stat.times(x_tick), slice_stat_optimal{t}.Profit(x_tick,1), 'r--',...
        stat.times(x_tick), slice_stat_price2{t}.Profit(x_tick,1), 'b-',...
        stat.times(x_tick), slice_stat_static{t}.Profit(x_tick,1), 'k-.');
    h(1).Parent.FontName = font_name;
    h(1).Parent.FontSize = font_size+2;
    h(1).LineWidth = line_width+0.5;
    h(2).LineWidth = line_width+0.5;
    h(3).LineWidth = line_width+0.5;
    h(3).Color = color.green;
    legend(legend_label_compact);
    ylabel('Profit');
    xlabel('Time (hour)');
    axis([stat.times(1), stat.times(x_tick(end)) slice_profit_limit(t,:)]);
    oldLabels = str2double(get(gca,'YTickLabel'));
    scale = 10^-3;newLabels = num2str(oldLabels*scale);
    set(gca,'YTickLabel',newLabels,'units','normalized');
    posAxes = get(gca,'position');
    textBox = annotation('textbox','linestyle','none','string',['x 10\it^{' sprintf('%d',log10(1./scale)) '}']);
    posAn = get(textBox,'position');
    set(textBox,'position',[posAxes(1) posAxes(2)+posAxes(4) posAn(3) posAn(4)],'VerticalAlignment','cap');
end 
