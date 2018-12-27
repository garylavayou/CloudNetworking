%% Performance Evaluation of HSR/HSR-RSV/FSR
% results: providing one group of output for each scheme.
%
%% HSR/HSR-RSV/Baseline (default)

%% FSR and Baseline
%{
lines.Labels = {'FSR', 'Baseline'};

%}
function dataplot3sa2(results, lines, metrics, idx_offset, options)
global figs;
if ~exist('figs', 'var') || isempty(figs) || ~isstruct(figs)
    figs = struct;
end
%% Arguments
if nargin <= 0 || isempty(results)
	load('Results/EXP042_e100.mat', 'results');
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
	lines.Markers = {'none', 'none', 's', 'x', '+', 'd', 'o'};
end
if ~isfield(lines, 'MarkerSizes')
	lines.MarkerSizes = [6 7 6 8 6 6 6];
end
if ~isfield(lines, 'Styles')
	lines.Styles = {'-',  '-.', '--'};
end
if ~isfield(lines, 'Colors')
	lines.Colors = [Color.Red; Color.MildGreen; Color.MildBlue; Color.Purple; Color.Black; Color.Gray];
end
if ~isfield(lines, 'Width')
	lines.Width = [1 1 1 1 1 1];
end
if nargin <= 2 || isempty(metrics)
	metrics = struct;
end
if ~isfield(metrics, 'Names')
	metrics.Names = {'ReVariables', 'ReFlows',...
		'Profit', 'Cost', 'Utilization', 'FairIndex'};
end
if ~isfield(metrics, 'Labels')
	metrics.Labels = {'# of Reconfigured Variables',...
		'# of Reconfigured Flows',...
		'Cummulated Profit',...
		'Reconfiguration Cost',...
		'Utilization Ratio',...
		'Fairness Index'};
end
if ~isfield(metrics, 'Titles')
	metrics.Titles= {'Number of Reconfigured Variables',...
		'Number of Reconfigured Flows',...
		'Profit of Slice',...
		'Reconfiguration Cost',...
		'Utilization of Slice Resources',...
		'Fairness Index'};
end
if nargin <= 3 || isempty(idx_offset)
	idx_offset = 50;
end
if nargin <= 4
	options = struct;
end
if ~isfield(options, 'bSavePlot')
	options.bSavePlot = false;
end
if options.bSavePlot
	if ~isfield(options, 'Filenames')
	options.filenames = {'numac-reconfig-vartime',...
		'numac-reconfigflow-vartime', ...
		'profit-vartime',...
		'cost-vartime',...
		'utilization-vartime',...
		'fairness-vartime'};
	end
	if ~isfield(options, 'Suffix')
		options.Suffix = '';
	end
end

%%
tx = (idx_offset+1):height(results.Dimconfig{1});
t = results.(lines.Sources{1}){1}.Time(tx);
if length(tx)>=10
    marker_index = round(linspace(1,length(tx), 10));
end
for i = 1:length(metrics.Names)
	metric = metrics.Names{i};
	title = metrics.Titles{i};
	if isfield(figs, metric) && isvalid(figs.(metric))
		figs.(metric).Children.delete;
		figure(figs.(metric));
	else
		figs.(metric) = figure('Name', title);
		switch metric
			case 'Utilization'
				figs.(metric).OuterPosition(3:4) = [362 367];
			otherwise
				figs.(metric).OuterPosition(3:4) = [360 380];  % [496 476];
		end
	end
	hold on;
	switch metric
		case {'ReVariables', 'ReFlows'}  % accumulated
			%% Number of Reconfiguration Variables and Flows
			num_lines = length(lines.Sources);
			data = cell(num_lines, 1);
			for j = 1:num_lines
				name = lines.Sources{j};
				data{j} = cumsum(results.(name){1}{tx,metric});
				plot(t, data{j}, lines.Styles{j});
			end
			axis tight;
		case {'Utilization', 'FairIndex'}
			for j = 1:num_lines
				plot(t, results.(name){1}{tx,metric}, lines.Styles{j});
			end
			xlim([t(1), t(end-1)]);
			if strcmpi(metric, 'FairIndex')
				ylim([0.55, 0.9]);
			else
				ylim([0.8,1]);
			end
		case {'Profit', 'Cost'}
			cum_cost = cell(num_lines, 1);
			cum_profit = cell(num_lines, 1);
			tx1 = idx_offset:height(results.(lines.Sources{1}){1});
			t_diff = diff(results.(lines.Sources{1}){1}{tx1,'Time'});
			for j = 1:num_lines
				name = lines.Sources{j};
				cost = results.(name){1}{tx1,'Cost'}.*results.(name){1}{tx1,'Interval'};
				profit = results.(name){1}{tx1,'Profit'}+results.(name){1}{tx1,'Cost'};
				% recover profit without reconfiguration cost
				cum_cost{j} = cumsum(cost(1:(end-1)));
				cum_profit{j} = cumsum(profit(1:(end-1)).*t_diff) - cum_cost{j};
				if strcmpi(metric, 'Profit')
					plot(t, cum_profit{j}, lines.Styles{j});
				else
					plot(t, cum_cost{j}, lines.Styles{j});
				end
			end
			xlim([t(1), t(end-1)]);
		otherwise
	end
	ax = gca;
	for j = 1:num_lines
		k = num_lines-j+1;
		ax.Children(j).Color = lines.Colors(k).RGB;
		ax.Children(j).LineWidth = lines.Width(k);
		ax.Children(j).MarkerIndices = marker_index;
		ax.Children(j).Marker = lines.Markers{k};
		ax.Children(j).MarkerSize = lines.MarkerSizes(k);
	end
	switch metric
		case 'Cost'
			ax.OuterPosition = [-0.02 0 1.09 1.01];	% ax.OuterPosition = [-0.04 -0.01 1.13 1.06];
		case 'Utilization'
			ax.OuterPosition = [0 0 1.09 1.04];
		case 'FairIndex'
			ax.OuterPosition = [0 0 1.08 1.04];
		otherwise
			ax.OuterPosition = [0 0 1.08 1.01];			% ax.OuterPosition = [-0.045 -0.01 1.14 1.04];
	end
	ylabel(metrics.Labels{i});
	xlabel('Time');
	switch metric
		case {'Utilization', 'FairIndex'}
			legend(lines.Labels, 'Location', 'southwest');
		otherwise	
			legend(lines.Labels, 'Location', 'northwest');
	end
	%% small figure
	switch metric
		case {'ReVariables','ReFlows'}
			if num_lines >= 3
				if strcmpi(metric, 'ReVariables')
					% axes('OuterPosition', [0.58, 0.20, 0.41, 0.39]);
					ax = axes('OuterPosition', [0.58, 0.21, 0.42, 0.41]);
				else
					% ax = axes('OuterPosition', [0.61, 0.24, 0.38, 0.35]);
					ax = axes('OuterPosition', [0.61, 0.24, 0.38, 0.35]);
					if exist('textBox', 'var')
						textBox.delete;
					end
				end
				hold on;
				for j = 1:num_lines-1
					plot(t, data{j}, lines.Styles{j});
				end
				ax.XTickLabel = {};
				for j = 1:num_lines-1
					k = num_lines-j;
					ax.Children(j).Color = lines.Colors(k).RGB;
					ax.Children(j).LineWidth = lines.Width(k);
				end
				axis tight
				if strcmpi(metric, 'ReFlows')
					scale = 10^3; scale_string = sprintf('\\times10^{%d}',log10(scale));
					ax.YTickLabel = num2str(str2double(ax.YTickLabel)/scale);
					textBox = text(0, 1.125, scale_string, 'Interpreter', 'tex', ...
						'Units', 'normalized', 'FontSize', 9);
					% posAxes = get(gca,'position');
					% textBox = annotation('textbox','linestyle','none','string',['x10^{' sprintf('%d',log10(1./scale)) '}']);
					% posAn = get(textBox,'position');
					% set(textBox,'position',[posAxes(1) posAxes(2)+posAxes(4)-0.02 posAn(3) posAn(4)],'VerticalAlignment','cap');
				end
			end
		case 'Profit'
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
	end
	if options.bSavePlot
		export_fig(fig_num_reconfig, ['Figures/', options.Filenames{i}, options.Suffix], ...
			'-pdf', '-transparent');
	end
end
