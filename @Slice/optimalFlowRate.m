%% Optimal Flow Rate
% optimize flow rate with link and node capacity constraints.
function [ net_profit ] = optimalFlowRate( this, options )
if nargin <= 1
    options.Display = 'final';
    options.Method = 'slice';
else
    if ~isfield(options, 'Display')
        options.Display = 'final';
    end
    if ~isfield(options, 'Method')
        options.Method = 'slice';
    end
end

NL = this.NumberVirtualLinks;
NN = this.NumberVirtualNodes;
NV = this.NumberVNFs;
NP = this.NumberPaths;

%% Coefficient for process-rate constriants
if strcmp(options.Method, 'normal')
    %% Coefficient
    % If VNF f is not used by a path p (the corresponding flow), there is no constraints
    % on |f x p|.    
    idx = 1:NP;
    for v = 1:NV
        vid = this.VNFList(v);
        delete_items = this.I_path_function(:,vid)==0;
        this.As_res(idx(delete_items),:) = [];
        idx = idx + nnz(this.I_path_function(:,vid));
    end
end
nnz_As = nnz(this.As_res) + NV*nnz(this.I_node_path) + nnz(this.I_edge_path);
num_lcon = this.num_lcon_res + NN + NL;
As = spalloc(num_lcon, this.num_vars, nnz_As);
As(1:this.num_lcon_res,:) = this.As_res;

%% Add link and node capacity constrant coefficient.
col_index = NP+(1:NN:((NP-1)*NN+1))';
col_index = repmat(col_index, 1, NV);
for v = 2:NV
    col_index(:,v) = col_index(:,v-1) + NN*NP;
end
col_index = col_index(:);
row_index = this.num_lcon_res;
for n = 1:NN
    As(row_index+n, col_index) = repmat(this.I_node_path(n,:),1, NV);
    col_index = col_index + 1;
end
row_index = (1:NL) + this.num_lcon_res + NN;
As(row_index, 1:NP) = this.I_edge_path; 

%% boundary
bs = [sparse(this.num_lcon_res,1); 
    this.VirtualNodes.Capacity;
    this.VirtualLinks.Capacity];
lbs = sparse(this.num_vars,1);

%% Feasible test of start point
% *Start Point*
x0 = zeros(this.num_vars,1);
z_min = min(this.VirtualNodes.Capacity)/(NP*NV);
x_min = min(this.VirtualLinks.Capacity)/NP;
switch options.Method
    case 'single-function'
        max_alpha_f = max(this.alpha_f);
    case 'normal'
        max_alpha_f = max(this.Parent.VNFTable.ProcessEfficiency);
    otherwise
        max_alpha_f = max(this.Parent.VNFTable{this.VNFList, 'ProcessEfficiency'});
end
if z_min >= max_alpha_f*x_min  % [x,z] is feasible
    x0(1:NP) = x_min;
    x0((NP+1):end) = z_min;
else
    x0(1:NP) = z_min/max_alpha_f;
    x0((NP+1):end) = z_min;    
end
% x0(1:NP) = 1;
% max_alpha_f = max(this.Parent.VNFTable.ProcessEfficiency(this.VNFList));
% x0((NP+1):end) = max_alpha_f;
if ~this.checkFeasible(x0)
    error('error: infeasible start point.');
end

fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.HessianFcn = @(x,lambda)Slice.fcnHessian(x,lambda,this);
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.SpecifyConstraintGradient = false;
if nargin >= 2 && isfield(options, 'Display')
    fmincon_opt.Display = options.Display;
else
    fmincon_opt.Display = 'notify-detailed';  %'iter';
end
[x, fval, exitflag] = fmincon(@(x)Slice.fcnNetProfit(x,this), ...
    x0, As, bs, [], [], lbs, [], [], fmincon_opt);
if exitflag < 1
    error('abnormal exit with flag %d.',exitflag);
end
if ~this.checkFeasible(x, options)
    error('error: infeasible solution.');
end

%% Additional Process
% * When compute node load, |z_npf| corresonding to |h_np = 0| has been set as zero.
% only path use a node and path use a function, the variable |z_npf| is nonzero.
% * too small components should be rounded.
this.x_path = x(1:NP);
switch options.Method
    case 'single-function'
        this.VNFList = 1:this.Parent.NumberVNFs;
        this.initializeState(this.alpha_f, 2);
        z_np = x(NP+1:end);
        z_np = reshape(this.I_node_path(:).*z_np, NN,NP);
        znpf = zeros(NN, NP, this.NumberVNFs);
        for v = 1:this.NumberVNFs
            af = this.Parent.VNFTable.ProcessEfficiency(v);
            idx_path = this.I_path_function(:,v)==1;
            p_slice = this.path_owner(idx_path);
            znpf(:, idx_path, v) = z_np(:, idx_path).*(af./this.alpha_f(p_slice))';
        end
        this.z_npf = znpf(:);
    case 'normal'
        this.VNFList = 1:this.Parent.NumberVNFs;    % recover the original number of VNFs
        znpf = x(NP+1:end);
        this.z_npf = zeros(NN*NP*this.NumberVNFs,1);
        nz = this.NumberVirtualNodes*this.NumberPaths;
        z_index = 1:nz;
        new_z_index = z_index;
        for v = 1:this.NumberVNFs
            %         if b_vnf(v)
            mask_npf = this.I_node_path.*this.I_path_function(:,v)';
            this.z_npf(new_z_index) = mask_npf(:).*znpf(z_index);
            z_index = z_index + nz;
            %         end
            new_z_index = new_z_index + nz;
        end
        this.Variables.z = this.z_npf;
    otherwise
        this.z_npf = x(NP+1:end);
end
clear znpf;
this.Variables.x = this.x_path;
this.Variables.z = this.z_npf;
this.Variables.z(this.Variables.z<10^-3) = 0;
this.Variables.x(this.Variables.x<10^-3) = 0;
options.Display = 'final';
if ~this.checkFeasible([], options)
    warning('optimalFlowRate: the rounding of variables with small quantity will make the solution infeasible.');
end
this.setPathBandwidth;
this.FlowTable.Rate = this.getFlowRate;
this.VirtualLinks.Load = this.getLinkLoad;
this.VirtualNodes.Load = this.getNodeLoad;
net_profit = -fval;
end

