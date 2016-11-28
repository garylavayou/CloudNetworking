%% Physical Network Speification
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
link_opt.cost = LinkCostOption.LengthDependent;
graph_data = PhysicalNetwork.LoadNetworkData(...
    NetworkModel.Sample1, link_opt, node_opt);
%% VNF Specification
VNF_data.Number = 6;            % number of VNF type
VNF_data.Model = VNFIntegrateModel.AllInOne;
VNF_data.RandomSeed = [20161031 20161101];
%% Network Slice Specification
number_slices = 3;
rand_seed = 20161031;
slice_opt = struct;
weight = [100000, 200000, 300000];
for i = 1:number_slices
    slice_opt(i).Type = SliceType.RandomSingleFlow;
    slice_opt(i).RandomSeed = rand_seed + i;
    slice_opt(i).Weight = weight(i);
%     slice_opt.NumberPaths = [];
%     slcie_opt.DelayConstraint = [];
%     slice_opt.VNFList = [];
end
%% Construct Network
PN = PhysicalNetwork(graph_data, VNF_data, slice_opt);
clearvars -except PN
% PN.plot;

%% Optimize
% PN.optimizeResourcePrice;
% PN.optimizeNetSocialWelfare;
% PN.singleSliceOptimization;
%% Test
% path_vars = [2; 1; 0.5];
% user_rate = PN.slices{1}.getFlowRate(path_vars);
% edge_load = PN.slices{1}.getLinkLoad(path_vars);
% node_vars = cell(3,1);
% for i = 1:length(node_vars)
%     node_vars{i} = sparse(4,3);
% end
% node_vars{1}(3,1) = 2;
% node_vars{1}(3,2) = 1;
% node_vars{1}(3,3) = 0.5;
% node_vars{2}(3,1) = 2;
% node_vars{2}(1,2) = 1;
% node_vars{2}(2,3) = 0.5;
% node_vars{3}(4,1) = 2;
% node_vars{3}(4,2) = 1;
% node_vars{3}(1,3) = 0.5;
% node_load = PN.slices{1}.getNodeLoad(node_vars);
% fprintf('\tFlow rate (%d): %G\n', [1:PN.slices{1}.NumberFlows; user_rate']);
% fprintf('\tEdge load (%d): %G\n', [1:PN.slices{1}.NumberVirtualLinks; edge_load']);
% fprintf('\tNode load (%d): %G\n', [1:PN.slices{1}.NumberVirtualNodes; node_load']);
