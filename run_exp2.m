%% Run a sequence of network slice
% comparison of methods, global optimization, resource pricing, static partitioning, and
% static slicing method.

%% Specification of Substrate Network
% Use |CostUnit| to set the basic cost parameters, the ratio between the node cost and
% link cost can be adjusted through |CostUnit|.
% Then we adjust the |Weight| of slices to let the resource utilization of the network be
% in a reasonable value. 
% Adjust the |capacity_factor|, so that the resource utilization of node and link is
% close.
link_opt.delay = LinkDelayOption.Random;
link_opt.cost = LinkCostOption.CapacityInverse;
link_opt.CostUnit = 80;
net_opt.delta = 0.5;

node_opt.model = NetworkModel.Sample2;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.capacity_factor = 1.5;     % [0.3; 0.5; 0.8; 1]
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 600;

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [0 20161101]; 

%% Arrival sequence confiuration
type.Index = [4; 5; 6]; 
type.Fixed = 1;
type.FixedCount = 1;      % Number of persistent slices: {1|2|3...}
% type.PartitionWeight = [1 2 4];
seed = 20161231;
arrival.Number = 10;
arrival.Interval = 0.7;     % Average of arrival interval.  

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
run_sequence_1;     % to be renamed.

data_plot;

% save('output_seq23.mat', 'profit_accurate', 'profit_approx', ...
% 'profit_stat', 'rate_stat', 'runtime', 'stat', 'utilization');