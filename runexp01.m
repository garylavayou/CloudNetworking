%% Run a sequence of network slice
% comparison of methods, global optimization, resource pricing, static
% partitioning, and static slicing method.
%
% Network Model: Sample1/Sample2/SD_RAN

%% Specification of Substrate Network
clearvars -except progress_bar;
clear global;
global DEBUG;
DEBUG = true;
node_opt.Model = NetworkModel.Sample1; % NetworkModel.Sample2
%% Control Variables
b_single_optimal = false;
b_price_adjust = true;
b_static_slice = false;
b_dual_decomp = false;
b_price_adjust2 = false;   % N = 14 bug
b_resource_part = false;
b_part_price = false;
b_save = false;
b_plot = false;

%% Parameters
% The parameters are determined in <run_test1.m>.
node_opt.CapacityModel = NodeCapacityOption.BandwidthProportion;
node_opt.CostModel = NodeCostOption.CapacityInverse;
node_opt.CapacityFactor = 1.5;     % [0.3; 0.5; 0.8; 1]
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
switch node_opt.Model
	case NetworkModel.Sample1
		node_opt.CostUnit = 500;
		link_opt.CostUnit = 150;
		link_opt.CapacityFactor = 30;
		link_opt.RandomSeed = 20171013;
		net_opt.ClassName = 'CloudNetwork';
	case NetworkModel.Sample2
		node_opt.CostUnit = 500;
		link_opt.CostUnit = 300;
		link_opt.CapacityFactor = 1000;
		link_opt.RandomSeed = 20171017;
		net_opt.ClassName = 'CloudNetwork';
	case NetworkModel.SD_RAN
		% description = 'average profit ratio, middle cost, low price';
		node_opt.CostUnit = 0.2;
		link_opt.CostUnit = 0.1;
		link_opt.RandomSeed = 20170421;
		net_opt.ClassName = 'CloudAccessNetwork';
end

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
switch node_opt.Model
	case NetworkModel.Sample1
		VNF_opt.Number = 4;            % number of VNF type
	case NetworkModel.Sample2
		VNF_opt.Number = 6;            % number of VNF type
	case NetworkModel.SD_RAN
		VNF_opt.Number = 6;            % number of VNF type
end
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];
net_opt.AdmitPolicy = 'reject-flow';
net_opt.Form = 'compact';
net_opt.PricingFactor = 1;
net_opt.Threshold = 'min';

%% Arrival sequence configuration
% Assume that the arrival interval is fixed, *only adjust the service
% lifetime and the arrival probability* of each type of slices to control
% the number of slices in the network.
% See also <RequestEvent>.
switch node_opt.Model
	case NetworkModel.Sample1
		type.Index = [11; 21; 31];
		seed_dynamic = 20161231;
		arrival.Interval = 0.25;     % Average of arrival interval.
		arrival.Number = 100;
	case NetworkModel.Sample2
		type.Index = [12; 22; 32];
		seed_dynamic = 20170415;   % 20161231
		arrival.Interval = 1/6;     % Average of arrival interval.
		arrival.Number = 100;
	case NetworkModel.SD_RAN
		type.Index = [13; 23; 33];
		seed_dynamic = 20170410;   % 20161231
		arrival.Number = 150;
		arrival.Interval = 1/5;     % Average of arrival interval.
		arrival.StopTime = 250;
end
type.Fixed = 1;           % corresponding index in type.Index;
type.FixedCount = 3;      % Number of persistent slices: {1|2|3...}
% type.PartitionWeight = [1 2 4];  >> Move to Ex

%% Algorithm options
% options.ProfitType = {'ApproximatePrice','AccuratePrice'};        move to Ex
% options.WelfareType = {'Accurate', 'Approximate'};

%%
run_sequence;

%% Output
% plot figure
if b_plot
	switch net_opt.Model
		case NetworkModel.Sample1
			x_tick = 1:100;
			fig3_ylimits = [200, 700;
				0, 1000;
				500, 2500];
			slice_profit_limit = [0, 2500; 0, 2500; 0, 2500];
		case NetworkModel.Sample2
			x_tick = 1:160;
			fig3_ylimits = [1000, 8000;
				1000, 8000;
				1000, 8000];
			slice_profit_limit = [0, 7000; 0, 7000; 0, 7000];
		case NetworkModel.SD_RAN
			x_tick = 1:200;
			welfare_limit = [0 200000];
			cost_limit = [0 5000];
			profit_limit = [0 35000];
			slice_profit_limit = [0, 5000; 0, 1500; 0, 10000];
	end
	xlimits = [stat.times(1), stat.times(x_tick(end))];
	% fig2_limits = [1, stat.times(x_tick(end)), 600, 3000];
	data_plot;
end

% Save data
if b_save
	switch net_opt.Model
		case NetworkModel.Sample1
			description = '1-1';
			expnum = '11';
		case NetworkModel.Sample2
			description = '1-2';
			expnum = '12';
		case NetworkModel.SD_RAN
			description = '1-3';
			expnum = '13';
	end
	description = sprintf('Experiment(No. %s): handling an arrival process\n',...
		description);
	description = sprintf('%sNetwork Topology = %s\n', description, net_opt.Model.char);
	description = sprintf('%sThreshold = %s\n', desciption, net_opt.Threshold);
	description = sprintf('%sLink Cost Unit = %g\nNode Cost Unit%g\n', ...
		desciption, link_opt.CostUnit, node_opt.CostUnit);
	description = sprintf('%sPricingFactor = %g\n', net_opt.PricingFactor);
	output_file_name = strcat('Results\EXP01', exp_num, ...
		'_LC', replace(num2str(link_opt.CostUnit), '.', ''),...
		'_NC', replace(num2str(node_opt.CostUnit), '.', ''),...
		'_PF', replace(num2str(net_opt.PricingFactor, '%02d'), '.', ''),...
		'_TH', net_opt.Threshold(1:3),...
		'.mat');
	save(output_filename, 'stat', 'stat_optimal', 'stat_price', ...
		'stat_static', 'description');
end
