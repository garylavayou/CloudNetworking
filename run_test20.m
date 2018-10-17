%% Physical Network Specification for Sample-2
% Comparing the social welfare of _singleSliceOptimization_ and
% _optimizeResourcePriceNew_.
% NOTE: the 'bCompact' options has little effect on computing complexity.

%% Specification of Substrate Network
clear variables;
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
link_opt.CapacityFactor = 1000;
link_opt.CostUnit = 300;
link_opt.RandomSeed = 20171013;
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

net_opt.PricingFactor = 1;  % 0.125|0.25|0.5|1|2|4
net_opt.AdmitPolicy = 'reject-flow';
net_opt.Form = 'compact';
net_opt.Threshold = 'min';

%% Construct Network
% Initialize substrate network
% add network slices
% Test type: 12 22 32
type.Index = [12; 22; 32];
type.Fixed = [1; 2; 3];
type.FixedCount = [3; 13; 27];      % Number of persistent slices: {1|2|3...}

%% 
net_opt.SlicingMethod = SlicingMethod.SingleNormal;
PN = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
PN.slice_template = Slice.loadSliceTemplate(type.Index);
link_capacity = PN.readLink('Capacity');
node_capacity = PN.readDataCenter('Capacity');
seed_dynamic = floor(now);
num_type = length(type.Index);
num_fix_type = length(type.Fixed);
for t = 1:num_fix_type
    for s = 1:type.FixedCount(t)
        slice_opt = PN.slice_template(type.Fixed(t));
        slice_opt.RandomSeed = seed_dynamic;
        seed_dynamic = seed_dynamic + 1;
        PN.AddSlice(slice_opt);
    end
end
output_optimal = PN.singleSliceOptimization(struct('bCompact', false));
    
fprintf('\tOptimal net social welfare (without pricing) is %.4e.\n', ...
    output_optimal.WelfareOptimal);
fprintf('\tOptimal net social welfare (with pricing) is %.4e.\n', output_optimal.Welfare);
fprintf('\tnet profit of each slice:\n');
fprintf('\t\t%f\n',output_optimal.Profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
fprintf('\t\t%f\n',output_optimal.Profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n',PN.utilizationRatio);
fprintf('\t\t(Node utilization: %.2G)\n', sum(PN.readDataCenter('Load')/sum(node_capacity)));
fprintf('\t\t(Link utilization: %.2G)\n\n', sum(PN.readLink('Load')/sum(link_capacity)));

%%
% output_price = PN.optimizeResourcePriceNew();    
% fprintf('\tNet social welfare is %.4e.\n', output_price.Welfare);
% % fprintf('\tnet profit of each slice:\n');
% % fprintf('\t\t%f\n',output_optimal.Profit(1:(end-1),:));
% fprintf('\tnet profit of substrate network:\n');
% fprintf('\t\t%f\n',output_price.Profit(end,:));
% fprintf('\tNetwork utilization ratio %f.\n',PN.utilizationRatio);
% fprintf('\t\t(Node utilization: %.2G)\n', sum(PN.readDataCenter('Load')/sum(node_capacity)));
% fprintf('\t\t(Link utilization: %.2G)\n\n', sum(PN.readLink('Load')/sum(link_capacity)));

