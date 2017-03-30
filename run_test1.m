%% Physical Network Speification for Sample-1
% This program used to test if the parameters are properly set.
%% Specification of Substrate Network
% Use |CostUnit| to set the basic cost parameters, the ratio between the node cost and
% link cost can be adjusted through |CostUnit|.
% Then we adjust the |Weight| of slices to let the resource utilization of the network be
% in a reasonable value. 
% Adjust the |capacity_factor|, so that the resource utilization of node and link is
% close.
link_opt.delay = LinkDelayOption.Random;
link_opt.cost = LinkCostOption.CapacityInverse;
link_opt.CostUnit = 50;
link_opt.CapacityFactor = 3;
net_opt.delta = 0.7;

node_opt.model = NetworkModel.Sample1;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.capacity_factor = 1.5;     % [0.3; 0.5; 0.8; 1]
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 300;

%% Specification of VNFs and Network Slices 
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 4;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161031 20161101];    

options.Display = 'off';
options.PricingFactor = 1;
options.ProfitType = {'AccuratePrice', 'ApproximatePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% Construct Network
% Initialize substrate network
% add network slices

PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
slice_type = 3;
PN.slice_template = Slice.loadSliceTemplate(slice_type);
link_capacity = PN.getLinkField('Capacity');
disp('link capacity:');
disp(link_capacity');
node_capacity = PN.getNodeField('Capacity');
disp('node capacity:');
disp(node_capacity');
node_staticcost = PN.getNodeField('StaticCost');
disp('node static cost:');
disp(node_staticcost');
slice_opt = PN.slice_template(1);
PN.AddSlice(slice_opt);

%% single slice optimize
[output_optimal, runtime] = PN.singleSliceOptimization(options);
fprintf('\nSingle Slice Optimization\n');
fprintf('\tOptimal net social welfare (without pricing) is %.4e(%.4e).\n', ...
    output_optimal.welfare_accurate_optimal, output_optimal.welfare_approx_optimal);
fprintf('\tOptimal net social welfare (with pricing) is %.4e(%.4e).\n', ...
    output_optimal.welfare_accurate, output_optimal.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_optimal.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_optimal.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n',PN.utilizationRatio);
fprintf('\tNode load: %.2G\n', sum(PN.getNodeField('Load')/sum(node_capacity)));
fprintf('\tLink load: %.2G\n', sum(PN.getLinkField('Load')/sum(link_capacity)));
fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
    mean(PN.getLinkField('UnitCost')), mean(PN.getNodeField('UnitCost')));
fprintf('\tRatio of unit node cost to unit link cost: %.2G.\n',...
    mean(PN.getNodeField('UnitCost'))/mean(PN.getLinkField('UnitCost')));
