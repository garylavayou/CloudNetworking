%%
% test an arrival process
%% Configurable parameters in this experiment
link_opt.CostUnit = 0.2;
node_opt.CostUnit = 0.1;
options.PricingFactor = 3;      % 1 | 2 | 3 ,
options.Threshold = 'average';

%% Save the configuration
% description = 'average profit ratio, middle cost, low price';
description = sprintf('profit threshold = %s, \nlink cost unit = %g',...
    options.Threshold, link_opt.CostUnit);
description = sprintf('%s, \nnode cost unit = %g, \nprice factor = %g', ...
    description, node_opt.CostUnit, options.PricingFactor);
output_file_name = 'output_3';
output_file_name = strcat(output_file_name, ...
    '_LC', replace(num2str(link_opt.CostUnit), '.', ''),...
    '_NC', replace(num2str(node_opt.CostUnit), '.', ''),...
    '_PF', replace(num2str(options.PricingFactor, '%02d'), '.', ''),...
    '_TH', options.Threshold(1:3),...
    '.mat');

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
link_opt.RandomSeed = 20170421;
link_opt.delay = LinkDelayOption.Random;
net_opt.NetworkClass = 'CloudNetwork';
net_opt.AdmitPolicy = 'reject-flow';

VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0]; 

%% Arrival sequence confiuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
type.Index = [7; 8; 9]; 
type.Fixed = 1;
type.FixedCount = 3;      % Number of persistent slices: {1|2|3...}
seed_dynamic = 20170410;   % 20161231
arrival.Number = 150;
arrival.Interval = 1/5;     % Average of arrival interval.  
arrival.StopTime = 250;

%% Algorithm options
options.ProfitType = {'ApproximatePrice','AccuratePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% Control Variables
options.Display = 'off';
b_single_optimal = false;
b_price_adjust1 = false;
b_price_adjust2 = false;
b_static_slice = true;
b_dual_decomp = false;
b_resource_part = false;
b_part_price = false;
experiment_length = arrival.StopTime * ...
    (b_single_optimal+b_price_adjust1+b_price_adjust2+b_static_slice+b_dual_decomp...
    +b_resource_part+b_part_price);
%%
run_sequence;

%%
x_tick = 1:200;
% fig2_limits = [1, stat.times(x_tick(end)), 600, 3000];
welfare_limit = [0 200000];
cost_limit = [0 5000];
profit_limit = [0 35000];
slice_profit_limit = [0, 5000; 0, 1500; 0, 10000];
% data_plot;
%%
% save('output_seq31.mat', 'profit_accurate', 'profit_approx', ...
%     'profit_stat', 'rate_stat', 'runtime', 'stat', 'utilization', 'type');

