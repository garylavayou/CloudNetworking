%% Compare the L1-Approximation and the original reconfiguration cost
% data = results.Dimconfig{1};
% data = results.Fastconfig{1};
% data = results.FastconfigReserve;
idx = 1;
data = results.Dimconfig{idx};
x = 2:height(data);
T = data.Time;
if exist('fig_l1approx', 'var') && fig_l1approx.isvalid
    figure(fig_l1approx);
else
    fig_l1approx = figure('Name', 'L1 Approximation Cost');
end
fig_l1approx.OuterPosition = [100   500   600   350];
ax = subplot(121);
h = plot(T(x), data.Cost(x), T(x), data.LinearCost(x), '--');
xlim([T(1), T(end)]);
xlabel('Time(s)');
ylabel('Reconfiguration Cost');
ax.OuterPosition = [0.01 0 0.48 0.9948];
legend({'Origin', 'L1-Appox'});
ax = subplot(122);
r = data.Cost(x)./data.LinearCost(x);
r(isnan(r)) = 1;
plot(T(x), r);
xlim([T(1), T(end)]);
ylim([0, 8]);
xlabel('Time(s)');
ylabel('Approximation Ratio');
ax.OuterPosition = [0.5100 0 0.51 0.9948];

%% Split view
idx = 4;
x = 2:height(results.Dimconfig{idx});
T = results.Dimconfig{idx}.Time;
if exist('fig_l1approx1', 'var') && fig_l1approx1.isvalid
    figure(fig_l1approx1);
    fig_l1approx1.Children.delete;
else
    fig_l1approx1 = figure('Name', 'L1 Approximation Cost');
end
fig_l1approx1.OuterPosition = [100   700   300   340];
ax = axes;
ax.OuterPosition = [0,0,1.07,1.04];
plot(T(x), results.Dimconfig{1}.Cost(x), T(x), results.Dimconfig{1}.LinearCost(x), '--');
xlim([T(1), T(end)]);
xlabel('Time(s)');
ylabel('Reconfiguration Cost');
legend({'Origin', 'L1-Appox'});
export_fig(fig_l1approx1, 'Figures/approx-reconfig-cost', '-pdf', '-transparent');
if exist('fig_l1approx2', 'var') && fig_l1approx2.isvalid
    figure(fig_l1approx2);
    fig_l1approx2.Children.delete;
else
    fig_l1approx2 = figure('Name', 'L1 Approximation Cost');
end
fig_l1approx2.OuterPosition = [450   700   300   340];
ax = axes;
ax.OuterPosition = [0,0,1.07,1.04];
ax.FontName = font_name;
r = results.Dimconfig{idx}.Cost(x)./results.Dimconfig{idx}.LinearCost(x);
r(isnan(r)) = 1;
plot(T(x), r);
xlim([T(1), T(end)]);
ylim([0,10]);
xlabel('Time(s)');
ylabel('Approximation Ratio');
export_fig(fig_l1approx2, 'Figures/approx-reconfig-cost-ratio', '-pdf', '-transparent');