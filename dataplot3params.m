%% Performance Metrics Varying with Parameters
% mode = 'vareta'|'varweight'|'varnumber'

%% For HSR/HSR-RSV/Baseline (Default)
%{
lines.Sources = {'Dimconfig','DimconfigReserve', 'DimBaseline'};  % 'Fastconfig','FastconfigReserve', 
lines.Lables = {'HSR', 'HSR-RSV', 'Baseline'};
metrics.Names = {'ReconfigRatio'}; % 'FairIndex', 'Utilization', 'Profit'
metrics.Labels = {'Reconfigure Ratio'}; % 'Fairness Index', 'Utilization Ratio', 'Profit';
options.Suffix = '';		
%}
%% For FSR/Baseline
% Modify following Arguments
%{
lines.Sources = {'Fastconfig', 'DimBaseline'};  % 'Fastconfig','FastconfigReserve', 
lines.Labels = {'FSR', 'Baseline'};
options.Suffix = 'fast'; % 'fastcomp' for no normalized version.
%}
function dataplot3params(mode, lines, metrics, idx_offset, options)
global figs;
if ~exist('figs', 'var') || isempty(figs) || ~isstruct(figs)
    figs = struct;
end
%%
if nargin <= 0 || isempty(mode)
	mode = 'vareta';
end
if nargin <= 1 || isempty(lines)
	lines = struct;
end
if ~isfield(lines, 'Sources')
	lines.Sources = {'Dimconfig','DimconfigReserve', 'DimBaseline'};
end
if ~isfield(lines, 'Labels')
	lines.Labels = {'HSR', 'HSR-RSV', 'Baseline'};
end
if ~isfield(lines, 'Markers')
	lines.Markers = {'x', 'o', 's'};
end
if ~isfield(lines, 'Styles')
	lines.Styles = {'-',  '--', '-.'};
end
if ~isfield(lines, 'Colors')
	lines.Colors = [Color.MildBlue; Color.Red; Color.MildGreen];
end
if nargin <= 2 || isempty(metrics)
	metrics = struct;
end
if ~isfield(metrics, 'Names')
	metrics.Names = {'ReconfigRatio'};
end
if ~isfield(metrics, 'Labels')
	metrics.Labels = {'Reconfigure Ratio'};
end
if nargin <= 3 || isempty(idx_offset)
	idx_offset = 1;
end
if nargin <= 4
	options = struct;
end
if ~isfield(options, 'Suffix')
	options.Suffix = '';		
end
options.bNormalizeRatio = true;
if strcmpi(options.Suffix, '-fastcomp')
	options.bNormalizeRatio = false;
end
if ~isfield(options, 'bSavePlot')
	options.bSavePlot = false;
end

%%
switch mode
	case {'vareta', 'var-eta'}
		if length(lines.Sources) == 2
			load('Results/EXP5024_vareta_fast.mat', 'results', 'etas');
		else
			load('Results/EXP5024_vareta.mat', 'results', 'etas');
		end
		mode = 'vareta';
		vars = etas;
		x_label = '\eta';
	case {'varweight', 'var-weight'}
		if length(lines.Sources) == 2
			load('Results/EXP5024_varweight_fast.mat', 'results', 'weight');
		else
			load('Results/EXP5024_varweight.mat', 'results', 'weight');
		end
		mode = 'varweight';
		vars = weight;
		x_label = 'flow weight';
	case {'varnumber', 'var-number'}
		if length(lines.Sources) == 2
			load('Results/EXP5024_varnumber_fast.mat', 'results', 'numberflow');
		else
			load('Results/singles/EXP4_OUTPUT241s0100.mat', 'results', 'numberflow');
% 			load('Results/EXP5024_varnumber.mat', 'results', 'numberflow');
		end
		mode = 'varnumber';
		vars = numberflow;
		x_label = 'number of flows';
end
% load(['Results/5024_', mode, '.mat'], 'results');		% TODO, update the file name.
for k = 1:length(metrics.Names)
    metric = metrics.Names{k};
    if isfield(figs, metric) && figs.(metric).isvalid
        figs.(metric).Children.delete;
        figure(figs.(metric));
    else
        figs.(metric) = figure('Name', metric);
        figs.(metric).OuterPosition = [100+10*k   500-10*k   350   350];
    end
    for j = 1:length(lines.Sources)
        name = lines.Sources{j};
				if (strcmpi(name, 'Baseline') || strcmpi(name, 'dimbaseline')) ...
					&& strcmpi(metric, 'ReconfigRatio') ...
					&& options.bNormalizeRatio
							continue;
				end
        if (strcmpi(name, 'Baseline') || strcmpi(name, 'dimbaseline')) ...
						&& strcmpi(mode, 'vareta')
					if strcmpi(metric, 'reconfigratio')
						mean_value = mean(results.(name){1}{idx,'ReVariables'}./...
							results.(name){1}{idx,'NumVariables'})*ones(length(vars),1);
					else
						mean_value = mean(results.(name){1}{idx,metric})*ones(length(vars),1);
					end
        else
            mean_value = zeros(length(results.(name)),1);
            for i = 1:length(results.(name))
                idx = (idx_offset+1):height(results.(name){i});
                if strcmpi(metric, 'reconfigratio')
									if options.bNormalizeRatio
										%% Normalize the reconfiguration ratio based on the Baseline
										mean_value(i) = mean(results.(name){i}{idx,'ReVariables'}./...
											results.(name){i}{idx,'NumVariables'})/...
											mean(results.DimBaseline{1}{idx,'ReVariables'}./...
											results.DimBaseline{1}{idx,'NumVariables'});
									else
										%% Not Normalize
										mean_value(i) = mean(results.(name){i}{idx,'ReVariables'}./...
                        results.(name){i}{idx,'NumVariables'});
									end
                else
                    mean_value(i) = mean(results.(name){i}{idx,metric});
                end
            end
        end
        if strcmpi(mode, 'vareta') 
            semilogx(vars, mean_value, lines.Styles{length(lines.Sources)-j+1});
        else
            plot(vars, mean_value, lines.Styles{length(lines.Sources)-j+1});
        end
        hold on;
    end
    h = figs.(metric).Children(end).Children;
    for j=1:length(h)
        h(j).Color = lines.Colors(length(h)-j+1).RGB;
        h(j).Marker = lines.Markers{length(h)-j+1};
        h(j).LineWidth = 1;
    end
    ax = h(1).Parent;
    if strcmpi(metric, 'Utilization')
        ax.YLim(2) = 1.1;
    end
    if strcmpi(mode, 'vareta')
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
        case 'vareta'
            ax.OuterPosition = [0 0.01 1.08 1.04];
        case 'varweight'
            ax.OuterPosition = [0 0.01 1.05 1.04];
        case 'varnumber'
            ax.OuterPosition = [0 0.01 1.08 1.04];
    end
    if strcmpi(metric, 'ReconfigRatio')
			if ~options.bNormalizeRatio
				legend(lines.Labels, 'Location', 'best');
			elseif strcmpi(mode, 'varweight')
				legend(lines.Labels, 'Location', 'northwest');
			else
				legend(lines.Labels, 'Location', 'northeast');
			end
    else
        legend(lines.Labels, 'Location', 'best');
    end
    if strcmpi(metric, 'Profit')
        if strcmpi(mode, 'varweight')
            ax.OuterPosition = [0 0.01 1.05 1.00];
        end
    end
    xlabel(x_label);
    ylabel(metrics.Labels{k});
    hold off;
		
		if strcmpi(metric, 'ReconfigRatio')
			switch mode
				case {'vareta','var-eta'}
					filename = 'Figures/reconfig-ratio-vareta';
				case {'var-number','varnumber'}
					filename = 'Figures/reconfig-ratio-varnum';
				case {'var-weight', 'varweight'}
					filename = 'figures/reconfig-ratio-varweight';
			end
			filename = [filename, options.Suffix]; %#ok<AGROW>
			if options.bSavePlot
				export_fig(figs.ReconfigRatio, filename, '-pdf', '-transparent');
			end
		end
end

end
%%
% Average profit over time, derived from accumulate profit devided by time.
function mean_profit(profit, reconfig_cost) %#ok<INUSD,DEFNU>

end
