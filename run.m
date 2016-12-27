%% Physical Network Speification
%% Specification of Substrate Network

node_opt.model = NetworkModel.Sample1;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.capacity_factor = 1;     % [0.3; 0.5; 0.8; 1]
node_opt.cost = NodeCostOption.Uniform;
node_opt.alpha = 1; % [1; 3]
                    % the ratio of unit node cost to unit link cost.
link_opt.delay = LinkDelayOption.BandwidthPropotion;
link_opt.cost = LinkCostOption.LengthDependent;
link_opt.delay2cost = 0.1;        % 0.1
net_opt.delta = 0.5;
%% Specification of VNFs and Network Slices 

VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.StaticCostOption = NodeStaticCostOption.Random;     % 
VNF_opt.static_cost_range = [0.1 0.3];
% VNF_opt.ProcessEfficiency is not set, using random value.
VNF_opt.RandomSeed = [20161031 20161101];       % the first seed is for random static cost, 
                                                % the second is for process efficiency
number_slices = 5;              % 10
rand_seed = 20161031;
slice_opt = struct;
weight = [10000, 20000, 30000, 20000, 40000, 10000, 20000, 30000, 10000, 20000]*10;
number_vnf = [1; 1; 1];
for i = 1:number_slices
    slice_opt(i).Pattern = FlowPattern.RandomMultiFlow;
    slice_opt(i).RandomSeed = rand_seed + i;
    slice_opt(i).Weight = weight(i);
    slice_opt(i).NumberPaths = 3;
    slice_opt(i).NumberFlows = 10;  % 1 3 10
%     slice_opt(i).NumberVNFs = 3;
%     slcie_opt.DelayConstraint = [];
%     slice_opt.VNFList = [];
end
%% Construct Network
% Initialize substrate network and add network slices.

PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
for i = 1:length(slice_opt)
    PN.AddSlice(slice_opt(i));
end
% clearvars -except PN node_opt link_opt VNF_opt slice_opt net_opt
PN.plot;
link_capacity = PN.getLinkField('Capacity');
node_capacity = PN.getNodeField('Capacity');
node_staticcost = PN.getNodeField('StaticCost');
%% Optimize
% * *Global optimization*

options.Display = 'off';
% options.Method = 'single-function';
options.PricingFactor = 1;
options.PercentFactor = 0.8;
tic;
[output_optimal, ss] = PN.singleSliceOptimization(options);
toc
fprintf('\nSingle Slice Optimization\n');
fprintf('\tOptimal net social welfare is %.4e(%.4e).\n', ...
    output_optimal.welfare_accurate, output_optimal.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_optimal.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_optimal.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%% 
% * *Dual Decomposition*

options.Display = 'off';
tic
[output_dual] = PN.optimizeNetSocialWelfare1(options);
toc
fprintf('\nOptimization with Dual Decomposition\n')
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_dual.welfare_accurate, output_dual.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_dual.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_dual.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%% 
% * *Dynamic Pricing*

options.Display = 'off';
tic
[output_price] = PN.optimizeResourcePrice([], options);
toc
fprintf('\nOptimization with Adjusting Resource Price\n')
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_price.welfare_accurate, output_price.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_price.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_price.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%% 
% * *Resource Partitioning*

options.Display = 'off';
tic
[output_part] = PN.resourcePartitionOptimization([], options);
toc
fprintf('\nOptimization with Resource Partition\n');
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_part.welfare_accurate, output_part.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_part.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_part.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);