%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The net social welfare;
% # The profit of each slice and the substrate network .
% # Flow rate of all flows in the network.
%%
% |argin|: existing field for output.
%		[stat, slice_stat, output] = calculateOutput(this, slices, fields, options)
function [stat, slice_stat, output] = calculateOutput(this, slices, fields, options)
if nargout == 0
	return;
end
if nargout >= 1
	stat = table;
end
if nargout >= 2
	slice_stat = struct;
end
if nargout >= 3
	output = struct;
end

if isempty(slices)
	slices = this.slices;
end
if isempty(fields)
	fields = {'all'};
elseif ischar(fields)
	fields = {fields};
end
if isa(options, 'Dictionary')
	options = options.data;
end
assert(isfield(options, 'PricingPolicy'), 'error: <PricingPolicy> must be specified.');
options.Stage = 'final';
options.bFinal = true;
options.bCompact = false;

Ns = length(slices);
for i = 1:length(fields)
	name = fields{i};
	if nargout >= 3
		if ismember(name, {'all', 'LinkPrice'})
			output.LinkPrice = this.readLink('Price');
		end
		if ismember(name, {'all', 'NodePrice'})
			output.NodePrice = this.readDataCenter('Price');
		end
		if ismember(name, {'all', 'LinkLoad'})
			output.LinkLoad = this.readLink('Load');
		end
		if ismember(name, {'all', 'NodeLoad'})
			output.NodeLoad = this.readDataCenter('Load');
		end
		if ismember(name, {'all', 'FlowRate'})
			output.FlowRate = [];
			for s = 1:Ns
				output.FlowRate = [output.FlowRate; slices(s).FlowTable.Rate];
			end
		end
		%% Net profit with offered price
		% For slices:       net_profit = utility - payment(price).
		% The price is calculated accoding to optimization procedure.
		if ismember(name, {'all', 'SliceProfit'})
			output.SliceProfit = zeros(Ns, 1);
			for s = 1:Ns
				output.SliceProfit(s) = slices(s).Optimizer.getProfit(options);
			end
		end
	end
	if nargout>= 1
		if ismember(name, {'all', 'Welfare'})
			stat.Welfare = 0;
			for s = 1:Ns
				% Calculate the net social welfare: the total utility less the total network cost.
				sl = slices(s);
				stat.Welfare = stat.Welfare + sl.Weight*sum(fcnUtility(sl.FlowTable.Rate));
			end
			stat.Welfare = stat.Welfare - this.totalCost();
		end
		if ismember(name, {'all', 'Profit'})
			% Calculate the profit of physical network
			stat.Profit = this.getSliceProviderProfit(slices, [], options);
		end
		if ismember(name, {'all', 'Cost'})
			stat.Cost = this.totalCost();
		end
		if ismember(name, {'all', 'Flows'})
			stat.Flows = this.NumberFlows;
		end
		if ismember(name, {'all', 'RejectFlows', 'ViolateSlices'})
			num_rejects = zeros(Ns,1);
			for s = 1:Ns
				num_rejects(s) = nnz(this.slices(s).FlowTable.Rate <= this.op.options.NonzeroTolerance);
			end
			if ismember(name, {'all', 'RejectFlows'})
				stat.RejectFlows = sum(num_rejects);
			end
			if ismember(name, {'all', 'ViolateSlices'})
				stat.ViolateSlices = nnz(num_rejects>0);
			end
		end
		if ismember(name, {'all', 'Slices'})
			type_list = [this.slice_template.Type];  % some types in the template might not present
			[~, counts] = count([this.slices.Type], type_list);
			stat.Slices = counts';
			stat.Properties.VariableDescriptions{strcmpi(stat.Properties.VariableNames, ...
				'Slices')} = strjoin(split(num2str(type_list)),'|');
		end
		if ismember(name, {'all', 'Runtime'})
			if isstruct(this.op.runtime)
				stat.Runtime = [this.op.runtime.Parallel, this.op.runtime.Serial];
			else
				stat.Runtime = [0, this.op.runtime];
			end
			stat.Properties.VariableDescriptions{strcmpi(stat.Properties.VariableNames, ...
				'Runtime')} = 'Parallel|Serial';
		end
		if ismember(name, {'all', 'Iterations'})
			stat.Iterations = this.op.iterations;
		end
		if ismember(name, {'all', 'Utilization'})
			[theta, t_link, t_node] = this.utilizationRatio(true);
			stat.Utilization = [theta.Mean, theta.Overall];
			stat.Properties.VariableDescriptions{strcmpi(stat.Properties.VariableNames, ...
				'Utilization')} = 'Mean|Overall';
			if nargout >= 3
				[~, t_node.Max, t_node.Min, t_node.Std] = this.nodeUtilization();
				output.NodeUtilization = t_node;
				[~, t_link.Max, t_link.Min, t_link.Std] = this.linkUtilization();
				output.linkUtilization = t_link;
			end
		end
	end
	
	if nargout >= 2
		type_list = [this.slice_template.Type];
		slice_stat = table;
		for t = 1:length(type_list)
			sidx = this.findSlice(type_list(t));
			if ismember(name, {'all', 'SliceProfit'})
				if isempty(sidx)
					warning off; % warning expanding table with default values.
					slice_stat{t, 'Profit'} = [0, 0, 0, 0];
					warning on;
				else
					slice_profit = zeros(length(sidx), 1);
					for s = 1:length(sidx)
						slice_profit(sidx(s)) = slices(sidx(s)).Optimizer.getProfit(options);
					end
					warning off; % warning expanding table with default values.
					slice_stat{t, 'Profit'} = [mean(slice_profit(sidx)), max(slice_profit(sidx)), ...
						min(slice_profit(sidx)), std(slice_profit(sidx))];
					warning on;
				end
			end
			if ismember(name, {'all', 'Rate'})
				if isempty(sidx)
					warning off; % warning expanding table with default values.
					slice_stat{t, 'Rate'} = [0, 0, 0, 0];
					warning on;
				else
					rate = zeros(sum([this.slices(sidx).NumberFlows]),1);
					fidx_offset = 0;
					for s = 1:length(sidx) % sidx is a row vector
						num_flow = this.slices(sidx(s)).NumberFlows;
						rate(fidx_offset + (1:num_flow)) = this.slices(sidx(s)).FlowTable.Rate;
						fidx_offset = fidx_offset + num_flow;
					end
					warning off; % warning expanding table with default values.
					slice_stat{t, 'Rate'} = [mean(rate), max(rate), min(rate), std(rate)];
					warning on;
				end
			end
		end
		if ismember(name, {'all', 'SliceProfit'})
			slice_stat.Properties.VariableDescriptions{...
				strcmpi(slice_stat.Properties.VariableNames, 'Profit')} = 'Mean|Max|Min|Std';
		end
		if ismember(name, {'all', 'Rate'})
			slice_stat.Properties.VariableDescriptions{...
				strcmpi(slice_stat.Properties.VariableNames, 'Rate')} = 'Mean|Max|Min|Std';
		end
	end
end
end