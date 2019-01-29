%% Visualize the output of exp022
%   plotexp022('EXP022_Sample1', {'Dual', 'Resource', 'DistrScale', 'FixPricing'}, {'Welfare', 'Profit', 'Utilization'}, struct('SliceMetrics', {{'Profit', 'Rate'}}, 'ConfigId', 4, 'bSave', true))
function plotexp022(filename, plot_methods, metrics, options)
if nargin <= 0 || isempty(filename)
	filename = 'Results/EXP022_Sample1.mat';
else
	if ~startsWith(filename, 'Results/')
		filename = ['Results/' filename];
	end
	if ~endsWith(filename, '.mat')
		filename = [filename '.mat'];
	end
end
load(filename, 'results', 'slice_results', 'output_results', 'methods', 'type');
if nargin <= 1 || isempty(plot_methods)
	plot_methods = methods.Name;
end
if nargin <= 2 || (iscell(metrics) && isempty(metrics))
	metrics = {'Welfare', 'Profit', 'Cost'};
elseif isempty(metrics)
	metrics = {};
elseif ischar(metrics)
	metrics = {metrics};
end
defaultopts = struct('SliceMetrics', {{'Profit'}}, 'ConfigId', 1, 'bSave', false);
if nargin <= 3
	options = defaultopts;
else
	options = setdefault(options, defaultopts);
end
if options.bSave
	fig_filename = split(filename, {'/','.'});
	fig_filename = fig_filename{end-1};
end
%% Generate labels
% legend_label_full = {'Price-SPP', 'Dynamic-Slicing','Dynamic-Slicing2',...
%     'Static-Slicing', 'Upperbound-SPP', 'Optimal-SPP'};
% legend_label_compact = {'PS', 'DS', 'SS', 'OS'};
num_methods = length(plot_methods);
labels = cell(num_methods,1);
for i = 1:num_methods
	method = plot_methods{i};
	switch method
		case 'Dual'
			labels{i} = 'DualPricing';
		case 'DistrScale'
			labels{i} = 'Searching';
		case 'Resource'
			labels{i} = 'UsageBased';
		case 'StaticScale'
			labels{i} = 'OnlineScaling';
		case 'StaticDual'
			labels{i} = 'OnlineDual';
		case 'FixPricing'
			labels{i} = 'FixPrice';
		otherwise
			error('error: unrecognized method ''%s''.', method);
	end
end
ylabels = cell(num_methods,1);
for i = 1:length(metrics)
	switch metrics{i}
		case 'Welfare'
			ylabels{i} = 'Social Welfare';
		case 'NodeUtilization'
			ylabels{i} = 'Node Utilization';
		case 'LinkUtilization'
			% Since the links are bi-directional, some links might not be utilized. Therefore, the
			% link utilization cannot approach 1.
			ylabels{i} = 'Link Utilization';
		otherwise
			ylabels{i} = metrics{i};
	end
end
	
%% Figure Properties
line_width = 1.5;
% font_name = 'Calibri';
% font_size = 12;
%'+' | 'o' | '*' | '.' | 'x' |'square' | 'diamond' | 'v' | '^' | '>' | '<' | 'pentagram' |
%'hexagram' |'none'. 
line_markers = {'x', 'o', 's', '+', '^', 'd'}; 
line_styles = {'-',  '--', '-.', ':', '-'};  
line_colors = [Color.MildBlue; Color.Red; Color.MildGreen; Color.Purple; Color.Orange];

x_range = 1:length(type.Count);
x = type.Count(x_range);
num_config = length(type.RatioSet);

%% Initialize Figures
global figs;
if isempty(figs)
	figs = struct;
	for t = 1:length(metrics)
		metric = metrics{t};
		figs.(metric) = matlab.ui.Figure.empty(num_config,0);
		for i = 1:num_config
			figs.(metric)(i) = figure('Name', [metric ' Config ' num2str(i)], 'Visible', 'off');
		end
	end
else
	for t = 1:length(metrics)
		metric = metrics{t};
		if ~isfield(figs, metric)
			figs.(metric) = matlab.ui.Figure.empty(num_config,0);
			for i = 1:num_config
				figs.(metric)(i) = figure('Name', [metric ' Config ' num2str(i)], 'Visible', 'off');
			end
		else
			for i = 1:num_config
				if length(figs.(metric)) < i || ~figs.(metric)(i).isvalid
					figs.(metric)(i) = figure('Name', [metric ' Config ' num2str(i)], 'Visible', 'off');
				else
					figs.(metric)(i).Children.delete();
				end
			end
		end
	end
end

%% Plot
for t = 1:length(metrics)
	metric = metrics{t};
	for i = 1:num_config
		fig = figs.(metric)(i);
		fig.OuterPosition = [600+20*i 400+20*i 380 380];
		figure(fig);
		hold on;
		h = matlab.graphics.chart.primitive.Line.empty(num_methods, 0);
		for j = 1:num_methods
			method = plot_methods{j};
			switch metric
				case 'Utilization'
					data = results(i).(method){x_range, metric}(:,2);  % overall utilization
				case {'NodeUtilization', 'LinkUtilization'}
					data = zeros(length(x_range),1);
					for xi = 1:length(data)
						data(xi) = output_results(i,xi).(method).(metric).('Overall'); % overall utilization
					end
				case 'Runtime'
					data = results(i).(method){x_range, metric}(:,1);
				otherwise
					data = results(i).(method){x_range, metric};
			end
			h(j) = plot(x, data);
			h(j).Marker = line_markers{j};
			h(j).LineStyle = line_styles{j};
			h(j).Color = line_colors(j).RGB;
			h(j).LineWidth = line_width;
		end
		ax = h(1).Parent;
		xlim([type.Count(1), type.Count(end)]);
		if contains(metric, 'Utilization')
			ylim([0, 1]);
		else
			ax.YLim(1) = 0;
		end
		legend(labels, 'Location', 'best');
		xlabel('Number of slices');
		ylabel(ylabels{t});
		ax.YAxis.Exponent = floor(log10(ax.YLim(end)));
		ax.XTick = type.Count;
		% ax.XTickLabel = string(type.Count);
		% 	h(1).Parent.FontName = font_name;
		% 	h(1).Parent.FontSize = font_size;
		hold off;
		pause(1);
		figname = ['Figures/', fig_filename, '_config', num2str(i), '_', metric];
		export_fig(fig, figname, '-pdf', '-transparent');
	end
end

%% Visualize individual slice metrics
if ~isempty(options.SliceMetrics)
	if ~isfield(figs, 'Slice')
		figs.Slice = struct;
	end
	for t = 1:length(options.SliceMetrics)
		metric = options.SliceMetrics{t};
		if ~isfield(figs.Slice, metric)
			for i = 1:length(type.Index)
				if type.RatioSet{options.ConfigId}(i) == 0
					figs.Slice.(metric)(i) = matlab.graphics.GraphicsPlaceholder;
				else
					figs.Slice.(metric)(i) = figure('Name', ...
						['Slice ' num2str(type.Index(i)) ' ' metric ' Config ' num2str(options.ConfigId)], ...
						'Visible', 'off');
				end
			end
		else
			for i = 1:length(type.Index)
				if type.RatioSet{options.ConfigId}(i) ~= 0
					if length(figs.Slice.(metric)) < i || ~figs.Slice.(metric)(i).isvalid
						figs.Slice.(metric)(i) = figure(...
							'Name', ['Slice ' num2str(type.Index(i)) ' ' metric ' Config ' num2str(options.ConfigId)], ...
							'Visible', 'off');
					else
						figs.Slice.(metric)(i).Children.delete();
					end
				end
			end
		end
	end
	%% Plot
	for t = 1:length(options.SliceMetrics)
		metric = options.SliceMetrics{t};
		for i = 1:length(type.Index)
			if type.RatioSet{options.ConfigId}(i) ~= 0
				fig = figs.Slice.(metric)(i);
				fig.OuterPosition = [600+20*i 400+20*i 400 380];
				figure(fig);
				h = matlab.graphics.chart.primitive.Line.empty(num_methods, 0);
				hold on;
				for j = 1:num_methods
					method = plot_methods{j};
					switch metric
						case {'Profit', 'Rate'}
							data = slice_results(options.ConfigId, i).(method){x_range,metric}(:,1);
						otherwise
							error('error: not support.');
					end
					h(j) = plot(x, data);
					h(j).Marker = line_markers{j};
					h(j).LineStyle = line_styles{j};
					h(j).Color = line_colors(j).RGB;
					h(j).LineWidth = line_width;
				end
				ax = h(1).Parent;
				xlim([type.Count(1), type.Count(end)]);
				legend(labels, 'Location', 'best');
				xlabel('Number of slices');
				ylabel(metric);
				if (ax.YLim(2) > 1000)
					ax.YAxis.Exponent = floor(log10(ax.YLim(end)));
				end
				ax.XTick = type.Count;
				hold off
				pause(1);
				figname = ['Figures/', fig_filename, '_config', num2str(options.ConfigId), ...
					'_type', num2str(type.Index(i)), '_', metric];
				export_fig(fig, figname, '-pdf', '-transparent');
			end
		end
	end
end

%		ax.YAxis.Exponent = floor(log10(ax.YLim(end)));
% Legacy code
%     oldLabels = str2double(get(gca,'YTickLabel'));
%     scale = 10^-3;newLabels = num2str(oldLabels*scale);
%     set(gca,'YTickLabel',newLabels,'units','normalized');
%     posAxes = get(gca,'position');
%     textBox = annotation('textbox','linestyle','none','string',['x 10\it^{' sprintf('%d',log10(1./scale)) '}']);
%     posAn = get(textBox,'position');
%     set(textBox,'position',[posAxes(1) posAxes(2)+posAxes(4) posAn(3) posAn(4)],'VerticalAlignment','cap');
