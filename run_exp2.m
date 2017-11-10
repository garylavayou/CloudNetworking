%% Run a sequence of network slice
% comparison of methods, global optimization, resource pricing, static partitioning, and
% static slicing method.
% load('output_seq21.mat')

%% Specification of Substrate Network
% The parameters are determined in <run_test2.m>.
clearvars -except progress_bar;
clear global;
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
link_opt.CostUnit = 300;
link_opt.CapacityFactor = 1000;
link_opt.RandomSeed = 20171017;
% net_opt.Delta = 0.5;

node_opt.Model = NetworkModel.Sample2;
node_opt.CapacityModel = NodeCapacityOption.BandwidthProportion;
node_opt.CostModel = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;
node_opt.CapacityFactor = 1.5;     % [0.3; 0.5; 0.8; 1]

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];
net_opt.ClassName = 'CloudNetwork';
net_opt.AdmitPolicy = 'reject-flow';
net_opt.PricingFactor = 1;
net_opt.Form = 'compact';
net_opt.Threshold = 'min';

%% Arrival sequence configuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
type.Index = [12; 22; 32];
type.Fixed = 1;
type.FixedCount = 3;      % Number of persistent slices: {1|2|3...}
% type.PartitionWeight = [1 2 4];  >> Move to Ex
seed_dynamic = 20170415;   % 20161231
arrival.Number = 100;
arrival.Interval = 1/6;     % Average of arrival interval.  

%% Control Variables
declare_info_level('Global', DisplayLevel.Off);
b_single_optimal = true;
b_price_adjust = true;
b_static_slice = true;
b_dual_decomp = false;
b_price_adjust2 = false;
b_resource_part = false;
b_part_price = false;
b_save = true;
%%
run_sequence;

%%
x_tick = 1:150;
% fig2_limits = [1, stat.times(x_tick(end)), 600, 3000];
xlimits = [stat.times(1), stat.times(x_tick(end))];
fig3_ylimits = [1000, 8000;
    1000, 8000;
    1000, 8000];
slice_profit_limit = [0, 7000; 0, 7000; 0, 7000];
data_plot;
%%
if b_save
    description = 'experiment (2-1): topology = Sample2, threshold = min.';
    save('Results\output_seq21.mat', 'stat', 'stat_optimal', 'stat_price', 'stat_static',...
        'description');
end