function [ z , dual_fval ] = node_subproblem(this, lambda_pf )
NS = this.NumberSlices;
gamma = cell(NS,1);
min_gamma = inf*ones(NS,1);
node_uc = this.getNodeField('UnitCost');
node_capacity = this.getNodeField('Capacity');
epsilon = this.staticNodeCost;
delta = this.Delta;
N = this.NumberNodes;
phis_n = epsilon*delta*(N-1)/this.totalNodeCapacity;

for s = 1:NS
    NNi = this.slices{s}.NumberVirtualNodes;
    NPi = this.slices{s}.NumberPaths;
    NVi = this.slices{s}.NumberVNFs;
    
    nz = NNi*NPi;
    z_index = NPi+(1:nz);
    g_npf = zeros(NNi*NPi*NVi,1);
    for f = 1:NVi
        g_npf(z_index) = (node_uc + phis_n) .* this.slices{s}.I_node_path...
            - lambda_pf(:,f)' .* this.slices{s}.I_node_path;
        z_index = z_index + nz;
    end    
    gamma{s} = reshape(g_npf, NNi, NPi, NVi);
    min_gamma(s) = min(g_npf);
end

[m, idx] = min(min_gamma);
z = cell(NS,1);
if m>=0
    for s = 1:NS
        NNi = this.slices{s}.NumberVirtualNodes;
        NPi = this.slices{s}.NumberPaths;
        NVi = this.slices{s}.NumberVNFs;
        z{s} = zeros(NNi*NPi*NVi,1);
    end
else
    count_min = 0;
    idx = cell(NS,1);
    for n = 1
    for s = 1:NS
        idx(s) = find(gamma{s}==m);
        count_min = count_min + length(idx);
    end
    z_mean = 
end
end

