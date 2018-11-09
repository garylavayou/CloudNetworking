%% Physical Network Specification for Sample-1
% This program used to test if the parameters of substrate network (links,
% nodes, and VNFs) are properly set.
%% Specification of Substrate Network
% # Fixed the link capacity;
% # Adjust the |CostUnit| of links and nodes, so that the average unit cost is close to a
% preset value (here, we set the value as 0.5 and 0.4);
% Then we adjust the |Weight| of slices to let the resource utilization of the network be
% in a reasonable value.
% In addition, adjust the |CapacityFactor|, so that the resource utilization of node and
% link is close.
clear variables;
global DEBUG;
DEBUG = true;
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
link_opt.CostUnit = 150;     % 150
link_opt.CapacityFactor = 30;
link_opt.RandomSeed = 20171012;
% net_opt.Delta = 0.7;   
 
node_opt.Model = NetworkModel.Sample1;
node_opt.CapacityModel = NodeCapacityOption.BandwidthProportion;
node_opt.CostModel = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;    % 500
node_opt.CapacityFactor = 1.5;     % [0.3; 0.5; 0.8; 1; 1.5]

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 4;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 20161031];

netopt.PricingFactor = 2;  % {1, when |CostUnit=(150,500)|}
netopt.AdmitPolicy = 'reject-flow';
opopt.Form = 'compact';

%% Construct Network
% Initialize substrate network
% add network slices
% Test type: 11 21 31
slice_type = 11;

%% control variables
b_static = true;
b_optimal = true;
b_repeat = true;

%% 
if b_optimal
    netopt.SlicingMethod = SlicingMethod.SingleNormal;
    PN = SimpleCloudNetwork(node_opt, link_opt, VNF_opt, netopt);
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
		PN.getOptimizer(opopt);
    while true && N > 0
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        PN.AddSlice(slice_opt);
				output = PN.singleSliceOptimization();
    
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
    netopt.SlicingMethod = SlicingMethod.StaticPricing;
    PN_static = SimpleCloudNetwork(node_opt, link_opt, VNF_opt, netopt);
    PN_static.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = PN_static.readLink('Capacity');
    node_capacity = PN_static.readDataCenter('Capacity');
    link_price = PN_static.getLinkCost * (1 + netopt.PricingFactor);
    node_price = PN_static.getNodeCost * (1 + netopt.PricingFactor);
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
		PN_static.getOptimizer(opopt);
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
