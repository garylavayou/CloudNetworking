%% Reconfiguration Ratio (Number of Variables)
plot_lines = {'Fastconfig', 'DimBaseline'};  % 'Fastconfig','FastconfigReserve', 
plot_metrics = {'ReconfigRatio'}; % 'FairIndex', 'Utilization', 'Profit'
legend_label = {'FSR', 'Baseline'};
marker = {'x', 'o', 's', 'd'};
line_style = {'-',  '--', '-.'};
color_set = [Color.MildBlue; Color.Red; Color.MildGreen; Color.Purple; Color.Black; Color.Gray];
idx_offset = 2;

if ~exist('figs', 'var') || ~isstruct(figs)
    figs = struct;
end
switch mode
    case 'var-eta'
        vars = etas;
        x_label = '\eta';
    case 'var-weight'
        vars = weight;
        x_label = 'flow weight';
    case 'var-number'
        vars = numberflow;
        x_label = 'number of flows';
end
y_label = {'Reconfigure Ratio', 'Fairness Index', 'Utilization Ratio', 'Profit'};
for k = 1:length(plot_metrics)
    metric = plot_metrics{k};
    if isfield(figs, metric) && figs.(metric).isvalid
        figs.(metric).Children.delete;
        figure(figs.(metric));
    else
        figs.(metric) = figure('Name', metric);
        figs.(metric).OuterPosition = [100+10*k   500-10*k   350   350];
    end
    for j = 1:length(plot_lines)
        name = plot_lines{j};
        if (strcmpi(name, 'Baseline') || strcmpi(name, 'dimbaseline')) && strcmpi(mode, 'var-eta')
            if strcmpi(metric, 'reconfigratio')
                mean_value = mean(results.(name){1}{idx,'NumberReconfigVariables'}./...
                        results.(name){1}{idx,'NumberVariables'})*ones(length(vars),1);
            else
                mean_value = mean(results.(name){1}{idx,metric})*ones(length(vars),1);
            end
        else
            mean_value = zeros(length(results.(name)),1);
            for i = 1:length(results.(name))
                idx = idx_offset:height(results.(name){i});
                if strcmpi(metric, 'reconfigratio')
                    mean_value(i) = mean(results.(name){i}{idx,'NumberReconfigVariables'}./...
                        results.(name){i}{idx,'NumberVariables'});
                else
                    mean_value(i) = mean(results.(name){i}{idx,metric});
                end
            end
        end
        if strcmpi(mode, 'var-eta') 
            semilogx(vars, mean_value, line_style{length(plot_lines)-j+1});
        else
            plot(vars, mean_value, line_style{length(plot_lines)-j+1});
        end
        hold on;
    end
    h = figs.(metric).Children(end).Children;
    for j=1:length(h)
        h(j).Color = color_set(length(h)-j+1).RGB;
        h(j).Marker = marker{length(h)-j+1};
        h(j).LineWidth = 1;
    end
    ax = h(1).Parent;
    if strcmpi(plot_metrics, 'Utilization')
        ax.YLim(2) = 1.1;
    end
    if strcmpi(mode, 'var-eta')
        ax.XLim = [vars(1) vars(end)];
        ax.XTick = vars;
        %         h(1).Parent.XTickLabel = split(num2str(vars));
        ax.TickLabelInterpreter = 'latex';
        ax.XTickLabel = {'$^1\!/\!_{32}$', '$^1\!/\!_{16}$', '$^1\!/\!_{8}$', ...
            '$^1\!/\!_{4}$', '$^1\!/\!_{2}$', '1', '2', '4', '8'};
    else
        xlim([vars(1), vars(end)])
    end
    switch mode
        case 'var-eta'
            ax.OuterPosition = [0 0.01 1.08 1.04];
        case 'var-weight'
            ax.OuterPosition = [0 0.01 1.05 1.04];
        case 'var-number'
            ax.OuterPosition = [0 0.01 1.08 1.04];
    end
    legend(legend_label, 'Location', 'best');
    if strcmpi(metric, 'Profit')
        if strcmpi(mode, 'var-weight')
            ax.OuterPosition = [0 0.01 1.05 1.00];
        end
    end
    xlabel(x_label);
    ylabel(y_label{k});
    hold off;
end
switch mode
    case 'var-eta'
        filename = 'Figures/reconfig-ratio-vareta-fastcomp';
    case 'var-number'
        filename = 'Figures/reconfig-ratio-varnum-fastcomp';
    case 'var-weight'
        filename = 'figures/reconfig-ratio-varweight-fastcomp';
end
export_fig(figs.ReconfigRatio, filename, '-pdf', '-transparent');
%%
% Average profit over time, derived from accumulate profit devided by time.
function mean_profit(profit, reconfig_cost)

end
%%
% # Save Results
% A group of experiments with varying parameters
%{
varname = split(mode, '-'); varname = varname{2};
description = sprintf('%s\n%s\n%s',...
    sprintf('Experiment 502%d: verify the influence of network settings to FSR (without warm-up phase).', type.Permanent),...
    sprintf('Topology=Sample-2, SliceType = %d (disable ad-hoc mode, enable dimensioning).', type.Index(type.Permanent)),...
    sprintf('variables = %s', varname)...
    );
output_name = sprintf('Results/EXP502%d_var%s_fast.mat', type.Permanent, varname);
save(output_name, 'description', 'results', 'EXPNAME', 'mode', 'etas', 'numberflow', 'weight', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'NUM_EVENT');
%}
