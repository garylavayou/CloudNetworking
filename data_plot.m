%%
% if b_single_optimal
%     %     plot(x, profit_approx.optimal(x), x, profit_accurate.optimal(x));
% end
% if b_price_adjust
% %     plot(x, profit_approx.price(x), x, profit_accurate.price(x));
% end
%%
if exist('fig_1', 'var') && fig_1.isvalid
    figure(fig_1);
else
    scrsz = get(groot,'ScreenSize');
    fig_1 = figure('OuterPosition',[scrsz(1)-5 scrsz(2)+34 scrsz(3)+10 scrsz(4)-34], ...
        'Name', 'Numerical Results of Slices and Network');
end
offset = 0;
x = (1+offset):(num_events-offset);
subplot(2,3,1);
plot(x, profit_accurate.optimal(x),'-r', ...
    x, profit_accurate.price(x), '-b', ...
    x, profit_accurate.price2(x), '-c', ...
    x, profit_accurate.part(x), '-g',...
    x, profit_accurate.static(x),'-m');
legend({'optimal', 'price', 'price2', 'partition', 'static'}, 'Location', 'northwest');
xlabel('Network Events (slice arrival/departure)');
ylabel('Net Social Welfare');
title('Network social welfare from different methods');
h = subplot(2,3,2);
devi1 = profit_accurate.optimal(x) - profit_accurate.price(x);
devi2 = profit_accurate.optimal(x) - profit_accurate.part(x);
devi3 = profit_accurate.optimal(x) - profit_accurate.static(x);
devi4 = profit_accurate.optimal(x) - profit_accurate.price2(x);
plot(x, devi1,'-r', x, devi2, '-b', x, devi3, '-g', x, devi4, '-c');
clear devi1 devi2 devi3 devi4;
lim = h.YLim;
lim(1) = 0;
h.YLim = lim;
legend({'diff\_price', 'diff\_part', 'diff\_static', 'diff\_price2'}, 'Location', 'northwest')
xlabel('Network Events (slice arrival/departure)');
ylabel('Difference of optimal objective value');
title('Distance to optimal net social welfare');
h = subplot(2,3,3);
plot(x, utilization.optimal(x), '-r', x, utilization.price(x), '-b', ...
    x, utilization.price2(x), '-c', ...
    x, utilization.part(x), '-g', x, utilization.static(x), '-m');
lim = h.YLim;
lim(2) = 1;
h.YLim = lim;
clear h lim;
legend({'optimal', 'price', 'price2', 'partition', 'static'},'Location', 'northeast')
xlabel('Network Events (slice arrival/departure)');
ylabel('Network Utilization');
title('Network utilization from different methods');
subplot(2,3,4)
plot(x, stat.number_slices(x,1), '-g', x, stat.number_slices(x,2), '-b', ...
    x, stat.number_slices(x,3), '-r', x, sum(stat.number_slices(x,:),2), '-k',...
    x, sum(stat.number_slices_static(x,:),2), '-c',...
    x, sum(stat.number_reject(x,:),2), '-m',...
    x, sum(stat.number_part_reject(x,:),2), '--b');
legend({'Type-1', 'Type-2', 'Type-3', 'Total', 'Total Static', 'Reject', 'Part-Reject'}, ...
    'Location', 'northwest')
xlabel('Network Events (slice arrival/departure)');
ylabel('Number of Slices');
title('Number of Slices');
subplot(2,3,5)
plot(x, stat.number_flows(x,1), '-b', x, stat.number_flows_static(x,1), '-r');
legend({'Others', 'Static'},'Location', 'northwest')
xlabel('Network Events (slice arrival/departure)');
ylabel('Number of Flows');
title('Total number of flows in the network');
subplot(2,3,6)
plot(x, runtime.optimal(x), '-b', ...
    x, runtime.price(x)./sum(stat.number_slices(x,:),2), '-r',...
    x, runtime.price2(x)./sum(stat.number_slices(x,:),2), '-c',...
    x, runtime.part(x)./sum(stat.number_slices(x,:),2), '-g',...
    x, runtime.static(x), '-m');
xlabel('Network Events (slice arrival/departure)');
ylabel('Run Time (seconds)');
title('Run time of the optimization process');
legend({'optimal', 'price', 'price2', 'part', 'static'}, 'Location', 'northwest');
clear x;

%% statistic of the profit of different type of slice.
if exist('fig_2', 'var') && fig_2.isvalid
    figure(fig_2);
else
    fig_2 = figure('Name', 'Comparison of Slices');
end
x = 50;
Y = zeros(num_type+1,5);
for i = 1:num_type
    Y(i,:) = [profit_stat.optimal{i}.Average(x) profit_stat.price{i}.Average(x) ...
        profit_stat.price{i}.Average(x) ...
        profit_stat.part{i}.Average(x) profit_stat.static{i}.Average(x)];
end
% Y(end,:) = [profit_stat.optimal{end}(x), profit_stat.price{end}(x), ...
%     profit_stat.part{end}(x), profit_stat.static{end}(x)];
Y(end,1) = profit_approx.optimal(x) ...
    - profit_stat.optimal{1}.Average(x)*stat.number_slices(x,1) ...
    - profit_stat.optimal{2}.Average(x)*stat.number_slices(x,2) ...
    - profit_stat.optimal{3}.Average(x)*stat.number_slices(x,3);
Y(end,2) = profit_approx.price(x) ...
    - profit_stat.price{1}.Average(x)*stat.number_slices(x,1) ...
    - profit_stat.price{2}.Average(x)*stat.number_slices(x,2) ...
    - profit_stat.price{3}.Average(x)*stat.number_slices(x,3);
Y(end,3) = profit_approx.price2(x) ...
    - profit_stat.price2{1}.Average(x)*stat.number_slices(x,1) ...
    - profit_stat.price2{2}.Average(x)*stat.number_slices(x,2) ...
    - profit_stat.price2{3}.Average(x)*stat.number_slices(x,3);
Y(end,4) = profit_approx.part(x) ...
    - profit_stat.part{1}.Average(x)*stat.number_slices(x,1) ...
    - profit_stat.part{2}.Average(x)*stat.number_slices(x,2) ...
    - profit_stat.part{3}.Average(x)*stat.number_slices(x,3);
Y(end,5) = profit_approx.static(x) ...
    - profit_stat.static{1}.Average(x)*stat.number_slices_static(x,1) ...
    - profit_stat.static{2}.Average(x)*stat.number_slices_static(x,2) ...
    - profit_stat.static{3}.Average(x)*stat.number_slices_static(x,3);
b = bar(Y);
b(1).Parent.XTickLabel = {'Type-1', 'Type-2', 'Type-3', 'Network'};
if b(1).Parent.YLim(1) < -1000
    b(1).Parent.YLim = [-1000, b(1).Parent.YLim(2)];
end
legend({'optimal', 'price', 'price2', 'partition', 'static'}, 'Location', 'northwest');
xlabel('Type of slices or network');
ylabel('Avearge profit of each slice/network');
title(sprintf('Comparison of profit (Event ID: %d)', x));
