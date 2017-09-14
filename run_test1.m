%% Physical Network Specification for Sample-1
% This program used to test if the parameters of substrate network (links, nodes, and
% VNFs) are properly set. 
%% Specification of Substrate Network
% # Fixed the link capacity;
% # Adjust the |CostUnit| of links and nodes, so that the average unit cost is close to a
% preset value (here, we set the value as 0.5 and 0.4);
% Then we adjust the |Weight| of slices to let the resource utilization of the network be
% in a reasonable value. 
% In addition, adjust the |capacity_factor|, so that the resource utilization of node and
% link is close.
link_opt.delay = LinkDelayOption.Random;
link_opt.cost = LinkCostOption.CapacityInverse;
link_opt.CostUnit = 150;     % 150
link_opt.CapacityFactor = 30;
net_opt.delta = 0.7;

node_opt.model = NetworkModel.Sample1;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;    % 500
node_opt.capacity_factor = 1.5;     % [0.3; 0.5; 0.8; 1; 1.5]

%% Specification of VNFs and Network Slices 
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 4;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 20161031];    

options.Display = 'off';
options.PricingFactor = 2;  % {1, when |CostUnit=(150,500)|}
options.ProfitType = {'AccuratePrice', 'ApproximatePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% Construct Network
% Initialize substrate network
% add network slices
% Test type: 1 2 3
slice_type = 1;

%% control variables
b_static_slice = false;
b_optimal = false;
b_static_slice_repeat = true;

%%
if b_static_slice
    PNs = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PNs.slice_template = Slice.loadSliceTemplate(slice_type);
    link_price = PNs.getLinkCost * (1 + options.PricingFactor);
    node_price = PNs.getNodeCost * (1 + options.PricingFactor);
    PNs.setLinkField('Price', link_price);
    PNs.setDataCenterField('Price', node_price);
    slice_opt = PNs.slice_template(1);
    slice_opt.method = 'static-slicing';
    slice_opt.admit_ploicy = 'reject-flow';
    sl = PNs.AddSlice(slice_opt);
    options.Method = 'slice-price';
    output = PNs.staticSlicing(sl, options);

    fprintf('\nStatic Slicing\n');
    fprintf('\tOptimal net social welfare (with pricing) is %.4e(%.4e).\n', ...
        output.welfare_accurate, output.welfare_approx);
    fprintf('\tnet profit of each slice:\n');
    display(output.profit(1:(end-1),:));
    fprintf('\tnet profit of substrate network:\n');
    display(output.profit(end,:));
    fprintf('\tNetwork utilization ratio %f.\n',PNs.utilizationRatio);
    fprintf('\t\t(Node utilization: %.2G)\n', sum(PNs.getDataCenterField('Load')/sum(node_capacity)));
    fprintf('\t\t(Link utilization: %.2G)\n', sum(PNs.getLinkField('Load')/sum(link_capacity)));    
end
%% single slice optimize
if b_optimal
    PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = PN.getLinkField('Capacity');
    node_capacity = PN.getDataCenterField('Capacity');
    slice_opt = PN.slice_template(1);
    PN.AddSlice(slice_opt);
    options.Method = 'normal';
    output = PN.singleSliceOptimization(options);

    fprintf('\nSingle Slice Optimization:\n');
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(PN.getLinkField('UnitCost')), mean(PN.getDataCenterField('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n',...
        mean(PN.getDataCenterField('UnitCost'))/mean(PN.getLinkField('UnitCost')));
    fprintf('\tOptimal net social welfare (without pricing) is %.4e(%.4e).\n', ...
        output.welfare_accurate_optimal, output.welfare_approx_optimal);
    fprintf('\tOptimal net social welfare (with pricing) is %.4e(%.4e).\n', ...
        output.welfare_accurate, output.welfare_approx);
    fprintf('\tnet profit of each slice:\n');
    display(output.profit(1:(end-1),:));
    fprintf('\tnet profit of substrate network:\n');
    display(output.profit(end,:));
    fprintf('\tNetwork utilization ratio %f.\n',PN.utilizationRatio);
    fprintf('\t\t(Node utilization: %.2G)\n', sum(PN.getDataCenterField('Load')/sum(node_capacity)));
    fprintf('\t\t(Link utilization: %.2G)\n', sum(PN.getLinkField('Load')/sum(link_capacity)));
end
%%
if b_static_slice_repeat
    PN_static = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN_static.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = PN_static.getLinkField('Capacity');
    node_capacity = PN_static.getDataCenterField('Capacity');
    link_price = PN_static.getLinkCost * (1 + options.PricingFactor);
    node_price = PN_static.getNodeCost * (1 + options.PricingFactor);
    PN_static.setLinkField('Price', link_price);
    PN_static.setDataCenterField('Price', node_price);
    slice_opt = PN_static.slice_template(1);
    slice_opt.method = 'static-slicing';
    slice_opt.admit_ploicy = 'reject-slice';
    seed = floor(now);
    options.Method = 'slice-price';
    fprintf('\nStatic Slicing Repeat:\n');
    while true
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        sl = PN_static.AddSlice(slice_opt);
        if isempty(sl)
            break;
        end
        output = PN_static.staticSlicing(sl, options);
        fprintf('\tNumber of slices: %d.\n', PN_static.NumberSlices);
        fprintf('\tOptimal net social welfare (with pricing) is %.4e(%.4e).\n', ...
            output.welfare_accurate, output.welfare_approx);
        fprintf('\tnet profit of each slice:\n');
        display(output.profit(1:(end-1),:));
        fprintf('\tnet profit of substrate network:\n');
        display(output.profit(end,:));
        fprintf('\tNetwork utilization ratio %f.\n',PN_static.utilizationRatio);
        fprintf('\t\t(Node utilization: %.2G)\n', sum(PN_static.getDataCenterField('Load')/sum(node_capacity)));
        fprintf('\t\t(Link utilization: %.2G)\n\n', sum(PN_static.getLinkField('Load')/sum(link_capacity)));
    end
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(PN_static.getLinkField('UnitCost')), mean(PN_static.getDataCenterField('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n\n',...
        mean(PN_static.getDataCenterField('UnitCost'))/mean(PN_static.getLinkField('UnitCost')));
end
