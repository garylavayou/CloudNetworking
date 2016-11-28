function [single_slice, rate, utility] = singleSliceOptimization( this )
%singleSliceOptimization Optimize the resource allocation in a single slice.
%   All service flows are placed in one slice.

% Merge slices into one single big slice.
NL = this.NumberLinks;
NN = this.NumberNodes;
NS = this.NumberSlices;
slice_data.adjacent = this.graph.Adjacent;
slice_data.link_map_S2P = (1:NL)';
slice_data.link_map_P2S = (1:NL)';
slice_data.node_map_S2P = (1:NN)';
slice_data.node_map_P2S = (1:NN)';
slice_data.flow_table = table([],[],[],[],[],[],[],'VariableNames',...
    {this.slices{1}.FlowTable.Properties.VariableNames{:,:},'Weight'});
% b_vnf = zeros(this.NumberVNFs, 1);
NF = 0;
for s = 1:NS
    NF = this.slices{s}.NumberFlows + NF;
end
flow_owner = zeros(NF, 1);
nf = 0;
for s = 1:NS
    new_table = this.slices{s}.FlowTable;
    % Map the virtual nodes to physical nodes.
    new_table.Source = this.slices{s}.VirtualNodes{new_table.Source, {'PhysicalNode'}};
    new_table.Target = this.slices{s}.VirtualNodes{new_table.Target, {'PhysicalNode'}};
    for f = 1:height(this.slices{s}.FlowTable)
        % path_list is handle object, is should be copyed to the new table.
        path_list = PathList(this.slices{s}.FlowTable{f,'Paths'});
        for p = 1:path_list.Width
            path = path_list.paths{p};
            path.node_list = this.slices{s}.VirtualNodes{path.node_list,{'PhysicalNode'}};
        end
        new_table{f,'Paths'} = path_list;
    end
    new_table.Weight = this.slices{s}.weight*ones(height(new_table),1);
    slice_data.flow_table = [slice_data.flow_table; new_table];
    flow_owner(nf+(1:this.slices{s}.NumberFlows)) = s;
    nf = nf + this.slices{s}.NumberFlows;
%     b_vnf(this.slices{i}.VNFList) = 1;
end
% slice_data.VNFList = find(b_vnf);
slice_data.VNFList = 1:this.NumberVNFs;
NV = length(slice_data.VNFList);
single_slice = Slice(slice_data);
single_slice.Parent = this;
% the flow id and path id has been allocated in each slice already, no need to reallocate.
single_slice.initializeState;

NP = single_slice.NumberPaths;
I_flow_function = zeros(NF, NV);
for f = 1:NF
    I_flow_function(f, this.slices{flow_owner(f)}.VNFList) = 1;
end
I_path_function = single_slice.I_flow_path'*I_flow_function;
num_vars = single_slice.num_vars;

%% Coefficient
% If VNF f is not used by a path p, there is no constraints on f x p.
idx = 1:NP;
for v = 1:NV
    delete_items = I_path_function(:,v)==0;
    single_slice.As_res(idx(delete_items),:) = []; 
    idx = idx + nnz(I_path_function(:,v));
end
nnz_As = nnz(single_slice.As_res) + NV*nnz(single_slice.I_node_path) ...
    + nnz(single_slice.I_edge_path);
num_lcon = single_slice.num_lcon_res + NN + NL;
As = spalloc(num_lcon, num_vars, nnz_As);
As(1:single_slice.num_lcon_res,:) = single_slice.As_res;

%% Add Node capacity constrant coefficient.
col_index = NP+(1:NN:((NP-1)*NN+1))';
col_index = repmat(col_index, 1, NV);
for v = 2:NV
    col_index(:,v) = col_index(:,v-1) + NN*NP;
end
col_index = col_index(:);
row_index = single_slice.num_lcon_res;
for n = 1:NN
    As(row_index+n, col_index) = repmat(single_slice.I_node_path(n,:),1, NV);
    col_index = col_index + 1;
end
%% Add Link capacity constrant coefficient.
row_index = (1:NL) + single_slice.num_lcon_res + NN;
As(row_index, 1:NP) = single_slice.I_edge_path; 

%%
bs = [sparse(nnz(I_path_function),1); 
    this.getNodeField('Capacity'); 
    this.getLinkField('Capacity')];
lbs = sparse(num_vars,1);

%% Feasible test of start point
% *Start Point*
% x0 = zeros(num_vars,1);
% z_min = min(this.getNodeField('Capacity'))/(NP*NV);
% x_min = min(this.getLinkField('Capacity'))/NP;
% max_alpha_f = max(this.VNFTable.ProcessEfficiency);
% if z_min >= max_alpha_f*x_min  % [x,z] is feasible
%     x0(1:NP) = x_min;
%     x0((NP+1):end) = z_min;
% else
%     x0(1:NP) = z_min/max_alpha_f;
%     x0((NP+1):end) = z_min;    
% end
x0 = zeros(num_vars,1);
x0(1:NP) = 1;
max_alpha_f = max(this.VNFTable.ProcessEfficiency);
x0((NP+1):end) = 1*NP*max_alpha_f;
if ~isempty(find(As*x0>bs+10^-6, 1))
    error('error: infeasible start point.');
end

fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.HessianFcn = @(x,lambda)PhysicalNetwork.fcnHessian(x,lambda,single_slice);
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.SpecifyConstraintGradient = false;
fmincon_opt.Display = 'notify-detailed';
% fmincon_opt.OptimalityTolerance = 10^-10;
fmincon_opt.ConstraintTolerance = 0;
[x, fval, exitflag, output] = fmincon(@(x)PhysicalNetwork.fcnNetWelfare(x,single_slice), ...
    x0, As, bs, [], [], lbs, [], [], fmincon_opt);
disp(output);
if exitflag == 1 || exitflag == 2
    if ~single_slice.checkFeasible(x)
        error('error: infeasible solution.');
    end
    %% Additional Process
    % too small components should be rounded.-
    %     x(x < 1) = 0;
    single_slice.x_path = x(1:NP);
    single_slice.z_npf = x(NP+1:end);
    %% when compute node load, z_npf corresonding to h_np = 0 has been set as zero.
    nz = single_slice.NumberVirtualNodes*this.NumberPaths;
    z_index = 1:nz;
    for v = 1:NV
        % only path use a node and path use a function, the variable z_npf is nonzero.
        mask_npf = single_slice.I_node_path.*I_path_function(:,v)';
        single_slice.z_npf(z_index) = mask_npf(:).*single_slice.z_npf(z_index);
        z_index = z_index + nz;
    end
    single_slice.Variables.x = single_slice.x_path;
    single_slice.Variables.z = single_slice.z_npf;
    single_slice.Variables.z(single_slice.Variables.z<10^-3) = 0;
    single_slice.Variables.x(single_slice.Variables.x<10^-3) = 0;
    if ~single_slice.checkFeasible
        warning('the rounding of variables with small quantity will make the solution infeasible.');
    end
    single_slice.setPathBandwidth(single_slice.Variables.x);
    single_slice.FlowTable.Rate = single_slice.getFlowRate(single_slice.Variables.x);
    single_slice.VirtualLinks.Load = single_slice.getLinkLoad(single_slice.Variables.x);
    single_slice.VirtualNodes.Load = single_slice.getNodeLoad(single_slice.Variables.z);
    rate = single_slice.FlowTable.Rate;
    utility = -fval;
%     fprintf('\tThe optimal net social welfare of the network: %G.\n', -fval);
    % delta_lambda.n = single_slice.getNodeLoad-this.getNodeField('Capacity')
    % delta_lambda.e = single_slice.getLinkLoad-this.getLinkField('Capacity')
else
    rate = [];
    utility = [];
    warning('abnormal exit with flag %d.',exitflag);
end
end