%%
line_width = [1 1 1 1 0.5];
color_set = [Color.Black; Color.MildGreen; Color.Purple; Color.MildBlue; Color.Red];
legend_label = {'Benchmark', 'RFFV', 'RFEV', 'Hybrid-1', 'Hybrid-2'};

%% Number of Reconfiguration
if exist('fig_num_reconfig', 'var') && fig_num_reconfig.isvalid
    figure(fig_num_reconfig);
else
    fig_num_reconfig = figure('Name', 'Number of Reconfiguration');
end
xlabel('Experiment Time');
% fig_num_reconfig.OuterPosition = [100 400 400 380];
for i = 1:length(etas)
    yyaxis('left');
    hl = plot(idx, results.Reconfig{idx,'ReVariables'}, '-.', ...
        idx, results.Fastconfig{i}{idx,'ReVariables'}, ':',...
        idx, results.Fastconfig2{i}{idx,'ReVariables'}, '.');
    if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
        hold on;
        plot(idx, results.Dimconfig{i}{idx,'ReVariables'},'-',...
            idx, results.Dimconfig{i}{idx,'ReVariables'},'+');
        hold off;
    end
    hl = hl(1).Parent.Children;
    nl = length(hl);
    for k=1:nl
        hl(k).Color = color_set(nl-k+1).RGB;
        hl(k).LineWidth = line_width(nl-k+1);
    end
    ylabel('Number of Reconfiguration');
    ylimmax = hl(1).Parent.YLim(2);
    yyaxis right
    hr = plot(idx, results.Reconfig{idx,'NumVariables'}, '--');
    hr(1).Color = Color.Orange.RGB;
    hr(1).LineWidth = 1;
    ylabel('Number of Variables');
    hr.Parent.YColor = Color.Orange.RGB;
    title(sprintf('\\eta=%.2f', etas(i)));
    legend([legend_label, {'Variable#'}], 'Location', 'northwest');
    ylimmax = max(ylimmax, hr.Parent.YLim(2));
    yyaxis left;
    ylim(hl(1).Parent, [0, ylimmax]);
    yyaxis right;
    ylim(hr.Parent, [0, ylimmax]);
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
for i = 1:length(etas)
    ratio_reconfig = results.Reconfig{idx,'ReVariables'}./...
        results.Reconfig{idx,'NumVariables'};
    ratio_fastconfig = results.Fastconfig{i}{idx,'ReVariables'}./...
        results.Fastconfig{i}{idx,'NumVariables'};
    ratio_fastconfig2 = results.Fastconfig2{i}{idx,'ReVariables'}./...
        results.Fastconfig2{i}{idx,'NumVariables'};
    hl = plot(idx, ratio_reconfig, '-.', ...
        idx, ratio_fastconfig, ':',...
        idx, ratio_fastconfig2, '.');
    if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
        ratio_dimconfig = results.Dimconfig{i}{idx,'ReVariables'}./...
            results.Dimconfig{i}{idx,'NumVariables'};
        ratio_dimconfig2 = results.Dimconfig2{i}{idx,'ReVariables'}./...
            results.Dimconfig2{i}{idx,'NumVariables'};
        hold on;
        plot(idx, ratio_dimconfig, '-',...
            idx, ratio_dimconfig2, '+');
        hold off;
    end
    hl = hl(1).Parent.Children;
    nl = length(hl);
    for k=1:nl
        hl(k).Color = color_set(nl-k+1).RGB;
        hl(k).LineWidth = line_width(nl-k+1);
    end
    ylabel('Ratio of Reconfiguration');
    xlabel('Experiment Time');
    % yyaxis right
    % hr = plot(idx, g_results.slice.num_flows(idx), 'r:',...
    %     idx, num_nonzeros_reconfig(idx), '-.',...
    %     idx, num_nonzeros_fastconfig(idx), '--');
    % hr(2).Color = [0.871, 0.49, 0];
    % hr(3).Color = [0.749, 0, 0.749];
    % ylabel('Number of Reconfiguration');
    % legend({'Ratio-Reconfig', 'Ratio-FastConfig', 'Flows', 'Reconfig', 'Fast Reconfig'}, ...
    %     'Location', 'northwest');
    ylim([0,1]);
    legend(legend_label, 'Location', 'best', 'Orientation', 'vertical');
    title(sprintf('\\eta=%.2f', etas(i)));
    pause(1);
end
%% Cost of Reconfiguration
if exist('fig_cost_reconfig', 'var') && fig_cost_reconfig.isvalid
    figure(fig_cost_reconfig);
else
    fig_cost_reconfig = figure('Name', 'Cost of Reconfiguration');
end
for i = 1:length(etas)
    hl = plot(idx, results.Reconfig{idx,'Cost'}*etas(i), '-.', ...
        idx, results.Fastconfig{i}{idx,'Cost'}, ':', ...
        idx, results.Fastconfig2{i}{idx,'Cost'}, '.');
    if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
        hold on;
        plot(idx, results.Dimconfig{i}{idx,'Cost'}, '-', ...
            idx, results.Dimconfig2{i}{idx,'Cost'}, '+');
        hold off;
    end
    hl = hl(1).Parent.Children;
    nl = length(hl);
    for k=1:nl
        hl(k).Color = color_set(nl-k+1).RGB;
        hl(k).LineWidth = line_width(nl-k+1);
    end
    ylabel('Cost');
    xlabel('Experiment Time');
    legend(legend_label, 'Location', 'northwest');
    title(sprintf('\\eta=%.2f', etas(i)));
    pause(1);
end

%% Number of Reconfigured Flows
if exist('fig_flow_reconfig', 'var') && fig_flow_reconfig.isvalid
    figure(fig_flow_reconfig);
else
    fig_flow_reconfig = figure('Name', 'Number of Reconfigured Flows');
end
for i = 1:length(etas)
    yyaxis left;
    hl = plot(idx, results.Reconfig{idx,'ReFlows'}, '-.', ...
        idx, results.Fastconfig{i}{idx,'ReFlows'}, ':', ...
        idx, results.Fastconfig2{i}{idx,'ReFlows'}, '.');
    if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
        hold on;
        plot(idx, results.Dimconfig{i}{idx,'ReFlows'}, '-', ...
            idx, results.Dimconfig2{i}{idx,'ReFlows'}, '+');
        hold off;
    end
    hl = hl(1).Parent.Children;
    nl = length(hl);
    for k=1:nl
        hl(k).Color = color_set(nl-k+1).RGB;
        hl(k).LineWidth = line_width(nl-k+1);
    end
    ylabel('Number of Reconfigured Flows');
    ylimmax(1) = hl(1).Parent.YLim(2);
    yyaxis right;
    hr = plot(idx, results.Reconfig{idx,'Flows'});
    hr(1).Color = Color.Orange.RGB;
    hr(1).LineWidth = 1;
    ylimmax(2) = hr(1).Parent.YLim(2);
    ylim([0, max(ylimmax)+2])
    ylabel('Number of Flows');
    hr.Parent.YColor = Color.Orange.RGB;
    yyaxis left;
    ylim([0, max(ylimmax)+2]);
    xlabel('Experiment Time');
    legend([legend_label,{'Flows #'}], 'Location', 'northwest');
    title(sprintf('\\eta=%.2f', etas(i)));
    pause(1);
end

%% Profit
if exist('fig_profit_reconfig', 'var') && fig_profit_reconfig.isvalid
    figure(fig_profit_reconfig);
else
    fig_profit_reconfig = figure('Name', 'Profit with Reconfiguration');
end
% fig_profit_reconfig.OuterPosition = [100 400 400 380];
xlabel('Experiment Time');
for i = 1:length(etas)
    yyaxis left;
    profit = results.Reconfig{idx,'Profit'} + (1-etas(i))*results.Reconfig{idx,'Cost'};
    hl = plot(idx, profit, '-.', ...
        idx, results.Fastconfig{i}{idx,'Profit'}, ':',...
        idx, results.Fastconfig2{i}{idx,'Profit'}, '.');
    if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
        hold on;
        plot(idx, results.Dimconfig{i}{idx,'Profit'}, '-',...
            idx, results.Dimconfig2{i}{idx,'Profit'}, '+');
        hold off;
    end
    hl = hl(1).Parent.Children;
    nl = length(hl);
    for k=1:nl
        hl(k).Color = color_set(nl-k+1).RGB;
        hl(k).LineWidth = line_width(nl-k+1);
    end
    %     hl(1).Parent.YLim(1) = 0;
    ylabel('Profit');
    yyaxis right;
    hr = plot(idx, results.Reconfig{idx,'Flows'}, 'r--');
    hr.Parent.YLim(1) = 0;
    ylabel('Number of Flows');
    legend([legend_label, {'Flow #'}], 'Location', 'northwest');
    title(sprintf('\\eta=%.2f', etas(i)));
    pause(1);
    % pause;
end