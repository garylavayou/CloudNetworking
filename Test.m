
link_uc = ss.Parent.getLinkField('UnitCost');
node_uc = ss.Parent.getNodeField('UnitCost');
link_load = ss.VirtualLinks.Load;
node_load = ss.VirtualNodes.Load;
epsilon = ss.Parent.staticNodeCost;
N = ss.Parent.NumberNodes;
theta = ss.Parent.networkUtilization(node_load, link_load);

%% single optimization
utility = sum(ss.FlowTable.Weight.*log(ss.FlowTable.Rate));
link_cost = dot(link_uc, link_load);
node_cost = dot(node_uc, node_load);
static_cost = epsilon*((N-1)*theta+1);
profit = utility - link_cost - node_cost - static_cost;
disp(profit);

%% single optimization node load
for f = PN.NumberVNFs
    for i = ss.FlowTable
end

%% dual decompoistion
link_load = PN.getLinkField('Load');
node_load = PN.getNodeField('Load');
link_cost = PN.getLinkCost;
node_cost = PN.getNodeCost;
theta = PN.networkUtilization(node_load, link_load);
static_cost = epsilon*((N-1)*theta+1);
utility = 0;
weight = [100000, 200000, 300000];
for s = 1:PN.NumberSlices
    S = PN.slices{s};
    utility = utility + weight(s)*sum(log(S.FlowTable.Rate));
end
profit = utility - link_cost - node_cost - static_cost;
disp(profit);

%%
znpf = reshape(PN.slices{1}.Variables.z, PN.slices{1}.NumberVirtualNodes, ...
    PN.slices{1}.NumberPaths, PN.slices{1}.NumberVNFs);
znpf = reshape(PN.slices{2}.Variables.z, PN.slices{2}.NumberVirtualNodes, ...
    PN.slices{2}.NumberPaths, PN.slices{2}.NumberVNFs);
znpf = reshape(PN.slices{3}.Variables.z, PN.slices{3}.NumberVirtualNodes, ...
    PN.slices{3}.NumberPaths, PN.slices{3}.NumberVNFs);
