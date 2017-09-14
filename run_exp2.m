%% Run a sequence of network slice
% comparison of methods, global optimization, resource pricing, static partitioning, and
% static slicing method.
load('output_seq21.mat')

%% Specification of Substrate Network
% The parameters are determined in <run_test2.m>.
link_opt.delay = LinkDelayOption.Random;
link_opt.cost = LinkCostOption.CapacityInverse;
link_opt.CapacityFactor = 1000;
link_opt.CostUnit = 300;
net_opt.delta = 0.5;

node_opt.model = NetworkModel.Sample2;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;
node_opt.capacity_factor = 1.5;     % [0.3; 0.5; 0.8; 1]

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0]; 

%% Arrival sequence configuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
type.Index = [4; 5; 6]; 
type.Fixed = 1;
type.FixedCount = 3;      % Number of persistent slices: {1|2|3...}
% type.PartitionWeight = [1 2 4];
seed = 20170415;   % 20161231   -> seed_dynamic
arrival.Number = 100;
arrival.Interval = 1/6;     % Average of arrival interval.  

%% Algorithm options
options.PricingFactor = 1;
options.ProfitType = {'ApproximatePrice','AccuratePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% Control Variables
options.Display = 'off';
b_single_optimal = true;
b_price_adjust = true;
b_static_slice = true;
b_dual_decomp = false;
b_price_adjust2 = false;
b_resource_part = false;
b_part_price = false;

%%
run_sequence;

%%
x_tick = 1:150;
% fig2_limits = [1, stat.times(x_tick(end)), 600, 3000];
xlimits = [stat.times(1), stat.times(x_tick(end))];
fig3_ylimits = [1000, 8000;
    1000, 8000;
    1000, 8000];
data_plot;
%%
save('output_seq21.mat', 'profit_accurate', 'profit_approx', ...
    'profit_stat', 'rate_stat', 'runtime', 'stat', 'utilization', 'type');
