%% Physical Network Specification for Sample-2
% See also <run_test2.m>. 
% This program tests with the network specified node (data center) capacity.
clear variables;
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
link_opt.RandomSeed = 20170421;

node_opt.Model = NetworkModel.Sample2;
node_opt.CostModel = NodeCostOption.CapacityInverse;
case_num = 1;
switch case_num
    case 1 % All nodes are VNF-capacble, Number of VNF-Nodes 15.
        node_opt.CapacityModel = NodeCapacityOption.BandwidthProportion;
        node_opt.CapacityFactor = 2;     
        node_opt.CostUnit = 500;
        link_opt.CapacityFactor = 1000;
        link_opt.CostUnit = 300;
    case 2 % Middle number nodes are VNF-capable, Number of VNF-Nodes 6.
        node_opt.CapacityModel = NodeCapacityOption.NetworkSpecified;
        node_opt.Capacity = [0, 4000, 0, 4000, 0, 4000, 0, ...
            0, 0, 6000, 0, 3000, 0, 0, 2000];  % user defined capacity;
        node_opt.CapacityFactor = 3;     
        node_opt.CostUnit = 500;
        link_opt.CapacityFactor = 1000;
        link_opt.CostUnit = 50;
    case 3 % Small number nodes are VNF-capable, Number of VNF-Nodes 3.
        node_opt.CapacityModel = NodeCapacityOption.NetworkSpecified;
        node_opt.Capacity = [0, 0, 5000, 0, 0, 0, 0, ...
            0, 5000, 0, 5000, 0, 0, 0, 0];  % user defined capacity;
        node_opt.CapacityFactor = 4;  
        node_opt.CostUnit = 500;
        link_opt.CapacityFactor = 1000;
        link_opt.CostUnit = 300;
end

%% Specification of VNFs and Network Slices 
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];    
net_opt.AdmitPolicy = 'reject-flow';
net_opt.PricingFactor = 1;      % used for static_slicing and single slice optimization
opopt.Form = 'compact';       % 'compact'|'normal'

%% Construct Network
% Initialize substrate network
% add network slices
% Test type: 44 54 64
slice_type = 44;

%% control variables
b_static = true;
b_optimal = true;
b_repeat = true;

%% 
if b_optimal
    net_opt.SlicingMethod = SlicingMethod.SingleNormal;
    PN = SimpleCloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN.slice_template = Slice.loadSliceTemplate(slice_type);
		PN.getOptimizer(opopt);
    link_capacity = PN.readLink('Capacity');
    node_capacity = PN.readDataCenter('Capacity');
    seed = floor(now);
    slice_opt = PN.slice_template(1);
    fprintf('\nSingle Slice Optimization:\n');
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(PN.readLink('UnitCost')), ...
        mean(PN.readDataCenter('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n\n',...
        mean(PN.readDataCenter('UnitCost'))/mean(PN.readLink('UnitCost')));
    N = 5;
    while true && N > 0
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        PN.AddSlice(slice_opt);
        %     output = PN.singleSliceOptimization();
        output = PN.singleSliceOptimization(struct('bCompact', false));
    
        fprintf('\tNumber of slices: %d.\n', PN.NumberSlices);
        fprintf('\tOptimal net social welfare (without pricing) is %.4e.\n', ...
            output.WelfareOptimal);
        fprintf('\tOptimal net social welfare (with pricing) is %.4e.\n', output.Welfare);
        fprintf('\tnet profit of each slice:\n');
        fprintf('\t\t%f\n',output.Profit(1:(end-1),:));
        fprintf('\tnet profit of substrate network:\n');
        fprintf('\t\t%f\n',output.Profit(end,:));
        fprintf('\tNetwork utilization ratio %f.\n',PN.utilizationRatio);
        fprintf('\t\t(Node utilization: %.2G)\n', sum(PN.readDataCenter('Load')/sum(node_capacity)));
        fprintf('\t\t(Link utilization: %.2G)\n\n', sum(PN.readLink('Load')/sum(link_capacity)));
        if ~b_repeat
            break;
        else
            N = N - 1;
				end
				for i = 1:PN.NumberSlices
					PN.slices{i}.initialize;
				end
    end
end
%%
if b_static
    net_opt.SlicingMethod = SlicingMethod.StaticPricing;
    PN_static = SimpleCloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN_static.slice_template = Slice.loadSliceTemplate(slice_type);
		PN_static.getOptimizer(opopt);
    link_capacity = PN_static.readLink('Capacity');
    node_capacity = PN_static.readDataCenter('Capacity');
    link_price = PN_static.getLinkCost * (1 + net_opt.PricingFactor);
    node_price = PN_static.getNodeCost * (1 + net_opt.PricingFactor);
    PN_static.writeLink('Price', link_price);
    PN_static.writeDataCenter('Price', node_price);
    slice_opt = PN_static.slice_template(1);
    if b_repeat
        slice_opt.AdmitPolicy = 'reject-slice';
        fprintf('\nStatic Slicing Repeat:\n');
    else
        fprintf('\nStatic Slicing:\n');
    end
    seed = floor(now);
    while true
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        sl = PN_static.AddSlice(slice_opt);
        if isempty(sl)
            break;
        end
        output = PN_static.staticSlicing(sl);
        fprintf('\tNumber of slices: %d.\n', PN_static.NumberSlices);
        fprintf('\tOptimal net social welfare (with pricing) is %.4e.\n', output.Welfare);
        fprintf('\tnet profit of each slice:\n');
        fprintf('\t\t%f\n',output.Profit(1:(end-1),:));
        fprintf('\tnet profit of substrate network:\n');
        fprintf('\t\t%f\n',output.Profit(end,:));
        fprintf('\tNetwork utilization ratio %f.\n',PN_static.utilizationRatio);
        fprintf('\t\t(Node utilization: %.2G)\n', ...
            sum(PN_static.readDataCenter('Load')/sum(node_capacity)));
        fprintf('\t\t(Link utilization: %.2G)\n\n', ...
            sum(PN_static.readLink('Load')/sum(link_capacity)));
        if ~b_repeat
            break;
        end
    end
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(PN_static.readLink('UnitCost')), ...
        mean(PN_static.readDataCenter('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n\n',...
        mean(PN_static.readDataCenter('UnitCost'))/...
        mean(PN_static.readLink('UnitCost')));
end
