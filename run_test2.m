%% Physical Network Specification for Sample-2
% This program used to test if the parameters are properly set.
%% Specification of Substrate Network
% # Fixed the link capacity;
% # Adjust the |CostUnit| of links and nodes, so that the average unit cost is close to a
% preset value (here, we set the value as 0.2 and 0.1);
% In addition, adjust the |CapacityFactor|, so that the resource utilization of node and
% link is close.
% # Then we adjust the |Weight| of slices to let the resource utilization of the network be
% in a reasonable value.
clear variables;
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
link_opt.CapacityFactor = 1000;
link_opt.CostUnit = 300;
link_opt.RandomSeed = 20171013;
% net_opt.delta = 0.5;

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

net_opt.PricingFactor = 1;
net_opt.AdmitPolicy = 'reject-flow';
net_opt.Form = 'compact';

%% Construct Network
% Initialize substrate network
% add network slices
% Test type: 12 22 32
slice_type = 12;

%% control variables
b_static = true;
b_optimal = true;
b_repeat = true;

%% 
if b_optimal
    net_opt.SlicingMethod = SlicingMethod.SingleNormal;
    PN = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN.slice_template = Slice.loadSliceTemplate(slice_type);
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
    end
end
%%
if b_static
    net_opt.SlicingMethod = SlicingMethod.StaticPricing;
    PN_static = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN_static.slice_template = Slice.loadSliceTemplate(slice_type);
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
