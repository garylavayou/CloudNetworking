%%
% Giving a combination of different types of slices, calculate the
% performace metrics.
%
% Network Model: Sample1/Sample2/SDRAN
% NormalSlice: slices with ordered service chains.
%
% See also <runtest06>.

%% Specification of Substrate Network
% exp_methods: character array or cell array of characters.
function results = runexp022(ioopts, exp_methods, options)
clearvars variables -except *results progress_bar
global DEBUG ITER_LIMIT;
DEBUG = false;
ITER_LIMIT = 30;
EXPNAME = sprintf('EXP022');
node_opt.Model = NetworkModel.Sample1;  % {Sample1, Sample2, SDRAN}

%% Physical Network Specification
switch node_opt.Model
  case NetworkModel.Sample1 
    % see also <runexp02>
    options.ClassName = 'NormalCloudNetwork';
    link_opt.CostModel = LinkCostOption.CapacityInverse;
		link_opt.CostUnit = 40; 
		link_opt.CapacityFactor = 30;    % scale the embedded capacity, see <loadNetworkData>.
		node_opt.CapacityModel = NodeCapacityOption.NetworkSpecified;
    node_opt.Capacity = [1000, 800, 600, 1000, 600, 800];
    node_opt.CapacityFactor = 1;     % [0.3; 0.5; 0.8; 1]
		node_opt.CostModel = NodeCostOption.CapacityInverse;
		node_opt.CostUnit = 100;        % 250 | 300 | 500
	case NetworkModel.Sample2  % data from <pre_slice_reconfigure>
		options.ClassName = 'NormalCloudNetwork';
		node_opt.CapacityModel = NodeCapacityOption.NetworkSpecified;
		node_opt.Capacity = [0, 4000, 0, 4000, 0, 4000, 0, ...
			0, 0, 6000, 0, 3000, 0, 0, 2000];  % user defined capacity;
		node_opt.CapacityFactor = 3;
		node_opt.CostModel = NodeCostOption.CapacityInverse;
		node_opt.CostUnit = 500;
		link_opt.CostModel = LinkCostOption.CapacityInverse;
		link_opt.CapacityFactor = 1000;
		link_opt.CostUnit = 50;
	case NetworkModel.SDRAN
		options.ClassName = 'NormalAccessNetwork';
		link_opt.CostUnit = 0.2;
		node_opt.CostUnit = 0.1;
end
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.RandomSeed = 20170423;

%% Specification of VNFs and Network Slices
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];
VNF_opt.Number = 10;            % number of VNF type
VNF_opt.ProcessEfficiency = [0.8, 2, 1.2, 1, 0.5, 2.5, 1.5, 1, 1, 1];  %TODO: connected with real data
VNF_opt.Ordered = true;
switch node_opt.Model
  case NetworkModel.Sample1
    type.Index = [41; 51; 61];
	case NetworkModel.Sample2
		type.Index = [42; 52; 62];
	case NetworkModel.SDRAN
		type.Index = [13; 23; 33];
end
% type.RatioSet = {[1;1;1],[4; 1; 1],[1; 4; 1], [1; 1; 4]};
type.RatioSet = {[1;0;0], [0;1;0], [0;0;1], [1;1;1]};
type.Count = 3:3:24;   % 3:3:24
type.StaticClass = {'NormalSlice'};

%% Algorithm options
options.PricingFactor = 2;      % 1 | 2 | 3 ,
options.Form = 'compact';       % {'compact'|'normal'}
options.NonzeroTolerance = 10^-3;
options.ConstraintTolerance = 10^-3;
options.PostProcessing = 'round';
options.AdmitPolicy = 'reject-flow';

%% Experiment Control
if nargin >= 1
	b_save = bitand(ioopts, hex2dec('0001'));
  b_append = bitand(ioopts, hex2dec('0002')); 
else
  b_save = false;
  b_append = true;
end
if nargin <= 1 || isempty(exp_methods)
  exp_methods = {'Dual', 'DistrScale', 'Resource', 'FixPricing', 'StaticDual'};
end
methods = table({},{}, SlicingMethod.empty(0,1), 'VariableNames', ...
	{'Name', 'Signature', 'SlicingMethod'});
% Dual: dual-pricing with dual-ADMM;
% DistrScale: scaling unit cost based on SP profit.
% Resource: adjust price based on resource usage and SP profit.
% StaticScale: online processing slices with DistrScale.
% StaticDual: online processing slices with Dual
% {StaticFix, SingleScale}
% TODO: replace with <SlicingMethod> enumeration, like <ReconfigMethod>.
if find(ismember(exp_methods, {'all', 'Dual'}),1)   % ScaleDual
  methods(end+1, :) = {'Dual', 'optimizeResourcePriceDual', SlicingMethod.DualPricing};
end
if find(ismember(exp_methods, {'all', 'Dual2'}),1)   % ScaleDual
  methods(end+1, :) = {'Dual2', 'optimizeResourcePriceDual2', SlicingMethod.DualPricing};
end
if find(ismember(exp_methods, {'all', 'DistrScale'}),1)
  methods(end+1, :) = {'DistrScale', 'optimizeResourcePriceScaling', SlicingMethod.AdjustPricing};
end
if find(ismember(exp_methods, {'all', 'DistrScale2'}),1)
  methods(end+1, :) = {'DistrScale2', 'optimizeResourcePriceScaling2', SlicingMethod.AdjustPricing};
end
if find(ismember(exp_methods, {'all', 'Resource'}),1) % Adapt
  methods(end+1, :) = {'Resource', 'optimizeResourcePrice', SlicingMethod.AdjustPricing};
  % # compares to the max-sp-profit mode.
  options.Threshold = 'off';      % 'min'|'average'|'max'|'off'
end
if find(ismember(exp_methods, {'all', 'AdaptDual'}),1)
  methods(end+1, :) = {'AdaptDual', 'pricingResourceDual', SlicingMethod.AdjustPricing};
end
if find(ismember(exp_methods, {'all', 'FixPricing'}),1)
  methods(end+1, :) = {'FixPricing', 'fixResourcePricing', SlicingMethod.AdjustPricing};
	options.PricingFactor = 5;
end
if find(ismember(exp_methods, {'all', 'StaticScale'}),1)
  methods(end+1, :) = {'StaticScale', 'staticSlicing', SlicingMethod.StaticPricing};
end
if find(ismember(exp_methods, {'all', 'StaticDual'}),1)
  methods(end+1, :) = {'StaticDual', 'staticSlicing', SlicingMethod.StaticPricing};
end
clear exp_methods;
ReplayMethods = {'Dual', 'Dual2', 'DistrScale', 'DistrScale2', 'Resource', 'SingleScale', ...
	'AdaptDual', 'FixPricing'};
NonReplayMethods = {'StaticScale', 'StaticFix', 'StaticDual', 'StaticAdaptDual'};
metrics = {'all'};
%b_single_optimal = false;

RANDSEED = 20170430;   % 20170410 20161231

%% Initialization
TOTAL_NUM = length(type.RatioSet)*length(type.Count)*height(methods);
setwaitbar;
num_config = length(type.RatioSet);
num_point = length(type.Count);
num_type = length(type.Index);
if b_save || b_append
	filename = sprintf('Results/%s_%s.mat', EXPNAME, node_opt.Model.char);
end
if b_append && exist(filename, 'file')
	load(filename, 'results', 'slice_results', 'output_results');
else
  results = struct; %#ok<*UNRCH>
	slice_results = struct;
	output_results = struct;
end

for exp_id = 1:num_config
	type.Ratio = type.RatioSet{exp_id};
	for call_id = 1:height(methods)
		method = methods{call_id, 'Name'}{1};
		invoked_method = methods{call_id, 'Signature'}{1};
		progress_bar.Name = method;
		waitbar(total_iter_num/TOTAL_NUM, progress_bar, ...
				sprintf('%s: %d/%d', invoked_method, total_iter_num, TOTAL_NUM));
		pause(0.01);
		% Network Initialization
		% Network must be recreate so as to obtain the same random number stream .
		seed_dynamic = RANDSEED;
		options.SlicingMethod = methods{call_id, 'SlicingMethod'};
		PN = instantiateclass(options.ClassName, node_opt, link_opt, VNF_opt, options);
		PN.slice_template = Slice.loadSliceTemplate(type.Index);
		PN.getOptimizer(options);
		if isequal(method, 'FixPricing')
			PN.writeLink('Price', PN.readLink('UnitCost')*options.PricingFactor);
			PN.writeDataCenter('Price', PN.readDataCenter('UnitCost')*options.PricingFactor);
		end
		number_slices = zeros(num_point,num_type);
		for point_id = 1:num_point
			fprintf('(%s) Configuration %d, Method %s, %d-th point.\n', ...
				datestr(now), exp_id, method, point_id);
			%% Specify the combination of slices
			% The combination of slices of the current point is based on that of the previous
			% point. So we should guarantee that the common part of the two point should be
			% the same, i.e. the common slices should have the same sequence and flow demand.
			total_num_slices = type.Count(point_id);
			for j = 1:(num_type-1)
				ct = round(total_num_slices * type.Ratio(j)/sum(type.Ratio));
				if ct > total_num_slices - sum(number_slices(point_id,:))
					ct = ct - 1;
				end
				number_slices(point_id,j) = ct;
			end
			number_slices(point_id,num_type) = total_num_slices - sum(number_slices(point_id,:));
			if point_id == 1
				delta_number_slices = number_slices(point_id,:);
			else
				delta_number_slices = number_slices(point_id,:) - number_slices(point_id-1,:);
				delta_number_slices(delta_number_slices<0) = 0;
			end
			rng(RANDSEED + point_id);
			delta_total = sum(delta_number_slices);
			delta_slice_sequence = unique_randi(delta_total, delta_total, 'stable');
			delta_slice_type = zeros(delta_total,1);
			t1 = 0;
			for j = 1:num_type
				t2 = t1 + delta_number_slices(j);
				delta_slice_type((t1+1):t2) = j;
				t1 = t2;
			end
			slices = Slice.empty(0,delta_total);
			if ismember(method, NonReplayMethods)
				% ~isempty(intersect(NonReplayMethods, method))
				s = 0;
				for j = 1:delta_total
					s_type = delta_slice_type(delta_slice_sequence(j));
					slice_opt = PN.slice_template(s_type);
					slice_opt.RandomSeed = seed_dynamic; seed_dynamic = seed_dynamic + 1;
					sl = PN.AddSlice(slice_opt);
					if ~isempty(sl)
						s = s + 1;
						slices(s) = sl;
					end
				end
			end
			if ismember(method, ReplayMethods)
				for j = 1:delta_total
					s_type = delta_slice_type(delta_slice_sequence(j));
					slice_opt = PN.slice_template(s_type);
					slice_opt.RandomSeed = seed_dynamic; seed_dynamic = seed_dynamic + 1;
					slices(j) = PN.AddSlice(slice_opt);
				end
			end
			fprintf('(%s) %s: adding %d slice %d.\n', ...
				datestr(now), method, length(slices), slice_opt.Type);
			%% Run optimization
			opts = struct('bCountTime', true, 'Fields', metrics);
			if contains(method, 'Dual')
				opts.TuningMethod = 'DualADMM';
				opts.PricingMethod = 'RandomizeCost';   % RandomizeCost | UniformCost
			end
			% 			if contains(method, 'Static')
			% 				opts.PricingMethod = 'NormalizeCost';   % RandomizeCost | UniformCost
			% 			end
			if ismember(method, ReplayMethods)
				slices = Slice.empty(0,1);  % optimize all slices.
			else
				opts.Invoker = method;
			end
			fprintf('(%s) %s Optimization with %s.\n', datestr(now), method, invoked_method);
			[stat, slice_stat, ...
				output_results(exp_id,point_id).(method)] = PN.(invoked_method)(slices, opts);
			% Save Results (exp_id, slice_type)
			results(exp_id).(method)(point_id, stat.Properties.VariableNames) = stat;
			for j = 1:num_type
				slice_results(exp_id,j).(method)(point_id, slice_stat.Properties.VariableNames) = ...
					slice_stat(j,:);
			end
			%% Update the progress
			total_iter_num = total_iter_num + 1;
			waitbar(total_iter_num/TOTAL_NUM, progress_bar, ...
				sprintf('%s: %d/%d', invoked_method, total_iter_num, TOTAL_NUM));
			pause(0.01);			
		end
		results(exp_id).(method).Properties.VariableDescriptions = ...
			stat.Properties.VariableDescriptions;
		for j = 1:num_type
			slice_results(exp_id,j).(method).Properties.VariableDescriptions = ...
				slice_stat.Properties.VariableDescriptions;
		end
	end
end
close(progress_bar);

%% Save experiment data
if b_save || b_append
	if b_append && exist(filename, 'file')
		new_methods = methods;
		load(filename, 'methods')
		m = setdiff(new_methods.Name, methods.Name);
		idx = strcmpi(m, new_methods.Name);
		methods = [methods; new_methods(idx,:)];
		save(filename, 'results', 'slice_results', 'output_results', 'methods', '-append');
	else
		infos = table({}, {}, 'VariableNames', {'Description', 'Details'});
		infos(end+1,:) = {'Experiment', 'Comparision of slice dimensioning methods'};
		infos(end+1,:) = {'Numbering', EXPNAME};
		infos(end+1,:) = {'Network Topology', node_opt.Model.char};
		infos(end+1,:) = ['Slice Type', type.StaticClass];
		infos(end+1,:) = {'SFC', 'Ordered'};
    s = inputdlg([], 'Please input experiment notes', [4,100]);
    if ~isempty(s)
      infos(end+1,:) = {'Notes', join(strip(string(s{1})), newline)};
    end
		save(filename, 'results', 'slice_results', 'output_results', 'infos', ...
			'node_opt', 'link_opt', 'VNF_opt', 'options', 'type', 'methods');
	end
end

%% Figure