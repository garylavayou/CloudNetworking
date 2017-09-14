%% Optimal Flow Rate
% optimize flow rate with link and node capacity constraints.
% Since the cpacity constraints is considered in the problem, the solution can utilize
% multipath to serve flows.
%% Methods
% * *normal*: combine all flows in one slice, and solve it as a global optimization;
% * *single-function*:
% * *slice*: solve each slice's problem independently;
% * *slice-price*: solve each slice's problem independently, with resource price known.
%% Model
% * *Accurate*
% * *Approximate*
% * *FixedCost*: used when slice's resource amount and resource prices are fixed, resource
% cost is fixed. See also <fcnSocialWelfare>, <DynamicSlice>.
function profit = optimalFlowRate( this, options )
if nargin <= 1
    options.Display = 'final';
    options.Method = 'slice';
else
    if ~isfield(options, 'Display')
        options.Display = 'final';
    end
end

NL = this.NumberVirtualLinks;
NC = this.NumberDataCenters;
NV = this.NumberVNFs;
NP = this.NumberPaths;

% Coefficient for process-rate constriants
if strcmp(options.Method, 'normal')
    %% Coefficient for global optimization
    % When all slices are combined into one slice, a VNF might not be used by all paths
    % (_i.e._ all flows). If a VNF |f| is not used by a path |p|, there is no
    % processing-rate constraints on $f \times p$(|delete_items|). To form the constraint
    % coefficient matrix, the related items in |As| should be removed. See also <Slice
    % file://E:/workspace/MATLAB/Projects/Documents/CloudNetworking/Slice.html>. 
    %
    % By the way, $z_{n,p,f}=0, \forall n$, if |p| does not use NFV |f|.
    % 
    %     delete_items = reshape(this.I_path_function, numel(this.I_path_function),1)==0;
    %     this.As_res(delete_items,:) = [];
    idx = 1:NP;
    for v = 1:NV
        vid = this.VNFList(v);
        delete_items = this.I_path_function(:,vid)==0;
        this.As_res(idx(delete_items),:) = [];
        idx = idx + nnz(this.I_path_function(:,vid));
    end
end
nnz_As = nnz(this.As_res) + NV*nnz(this.I_node_path) + nnz(this.I_edge_path);
num_lcon = this.num_lcon_res + NC + NL;
As = spalloc(num_lcon, this.num_vars, nnz_As);
As(1:this.num_lcon_res,:) = this.As_res;

%% Add node and link capacity constrant coefficient
% For the node capacity constraint, the coeffiect matrix is filled row-by-row.
% On row i, the non-zero elements located at (i-1)+(1:NC:((NP-1)*NC+1)) for the first
% |NC*NP| rows, and then the first |NC*NP| row is duplicated for |NV| times, result in
% |NC*NP*NV| rows. 
col_index = NP+(1:NC:((NP-1)*NC+1))';
col_index = repmat(col_index, 1, NV);
for v = 2:NV
    col_index(:,v) = col_index(:,v-1) + NC*NP;
end
col_index = col_index(:);
row_index = this.num_lcon_res;
for n = 1:NC
    As(row_index+n, col_index) = repmat(this.I_node_path(n,:),1, NV);
    col_index = col_index + 1;
end
row_index = (1:NL) + this.num_lcon_res + NC;
As(row_index, 1:NP) = this.I_edge_path; 

%% Boundary
% * *Resource Constraints*: processing-rate, the right side is all-zero; capacity
% constraints, the right side is the capcaity of link and node;
bs = [sparse(this.num_lcon_res,1); 
    this.VirtualDataCenters.Capacity;
    this.VirtualLinks.Capacity];
lbs = sparse(this.num_vars,1);
%%%
% * *Upper Bound*: Not neccesary, to facilitate the algorithm, we give a relaxed
% upper-bound.
ub = [max(this.Parent.getLinkField('Capacity'))*ones(NP,1);...
    max(this.Parent.getDataCenterField('Capacity'))*ones(this.num_vars-NP,1)];
%%%
% * Remove the capacity constraints with infinity capacity;
idx = find(bs==inf);
As(idx,:) = [];         
bs(idx) = [];

%% Feasible test of start point
% * *Start Point*: in case that the capacity of a virtual link/node is zero, we initialize
% $z_{min}$ and $x_{min}$ as the nonzero minimum value.
x0 = zeros(this.num_vars,1);
switch options.Method
    case 'single-function'
        max_alpha_f = max(this.alpha_f);
    case 'normal'
        max_alpha_f = max(this.Parent.VNFTable.ProcessEfficiency);
    otherwise
        max_alpha_f = max(this.Parent.VNFTable{this.VNFList, 'ProcessEfficiency'});
end
% z_min = min(this.VirtualNodes.Capacity(this.VirtualNodes.Capacity>1))/(NP*NV);
% x_min = min(this.VirtualLinks.Capacity(this.VirtualLinks.Capacity>1))/NP;
% if z_min == inf || x_min == inf
%     x0(1:NP) = 1;
%     x0((NP+1):end) = max_alpha_f;
% else
%     if z_min >= max_alpha_f*x_min  % [x,z] is feasible
%         x0(1:NP) = x_min;
%         x0((NP+1):end) = z_min;
%     else
%         x0(1:NP) = z_min/max_alpha_f;
%         x0((NP+1):end) = z_min;
%     end
% end
x0(1:NP) = 1;
x0((NP+1):end) = max_alpha_f;
if ~this.checkFeasible(x0)
    error('error: infeasible start point.');
end

fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.HessianFcn = @(x,lambda)Slice.fcnHessian(x,lambda,this);
fmincon_opt.SpecifyObjectiveGradient = true;
if nargin >= 2 && isfield(options, 'Display')
    fmincon_opt.Display = options.Display;
else
    fmincon_opt.Display = 'notify-detailed';  %'iter';
end
if strfind(options.Method, 'price')
    % output
    [x, fval, exitflag] = fmincon(@(x)Slice.fcnProfit(x,this), ...
        x0, As, bs, [], [], lbs, ub, [], fmincon_opt);
elseif strcmp(options.Model, 'FixedCost')
    [x, fval, exitflag] = fmincon(@(x)Slice.fcnSocialWelfare(x,this, options.Model), ...
        x0, As, bs, [], [], lbs, ub, [], fmincon_opt);
else
    [x, fval, exitflag] = fmincon(@(x)Slice.fcnSocialWelfare(x,this,'Approximate'), ...
        x0, As, bs, [], [], lbs, ub, [], fmincon_opt);
end
% x is a local solution to the problem when exitflag is positive.
if exitflag == 0
    warning('reaching maximum number of iterations.');
elseif exitflag < 0
    error('abnormal exit with flag %d.',exitflag);
elseif exitflag ~= 1
%     disp('residual link capacity:');
%     disp(this.getLinkLoad(this.x_path) - this.VirtualLinks.Capacity);
%     disp('residual node capacity:');
%     disp(this.getNodeLoad(this.z_npf) - this.VirtualNodes.Capacity);
    warning('(exitflag = %d) local optimal solution found.', exitflag);
end
if ~this.checkFeasible(x, options)
    error('error: infeasible solution.');
end

%% Additional Process
% Under the problem formulation, if node |n| is not used by path |p|($h_{np} = 0$) or path
% |p| does not use VNF |f|, the corresponding $z_{npf}$ should be zero. However, the
% optimization results of the related $z_{npf}$ can be arbitary value, since the related 
% $z_{npf}$ is not counted in the cost items. As a result, we shoulld mallualy set
% the value of those related $z_{npf}$ to zero.
%
% On the other hand, too small components should be rounded.
this.x_path = x(1:NP);
switch options.Method
    case 'single-function'
        this.VNFList = 1:this.Parent.NumberVNFs;
        this.initializeState(this.alpha_f, 2);
        z_np = x(NP+1:end);
        z_np = reshape(this.I_node_path(:).*z_np, NC,NP);
        znpf = zeros(NC, NP, this.NumberVNFs);
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
        this.z_npf = zeros(NC*NP*this.NumberVNFs,1);
        nz = this.NumberDataCenters*this.NumberPaths;
        z_index = 1:nz;
        new_z_index = z_index;
        for v = 1:this.NumberVNFs
            %         if b_vnf(v)
            mask_npf = this.I_node_path.*this.I_path_function(:,v)';    % compatiable arithmetic operation
            this.z_npf(new_z_index) = mask_npf(:).*znpf(z_index);
            z_index = z_index + nz;
            %         end
            new_z_index = new_z_index + nz;
        end
    otherwise
        this.z_npf = x(NP+1:end);
        nz = this.NumberDataCenters*this.NumberPaths;
        z_index = 1:nz;
        for f = 1:this.NumberVNFs
            this.z_npf(z_index) = this.I_node_path(:).*this.z_npf(z_index);
            z_index = z_index + nz;
        end
end
clear znpf;
this.flow_rate = this.getFlowRate(this.x_path);
if nargout == 1    % final results
    this.Variables.x = sparse(this.x_path);
    this.Variables.z = sparse(this.z_npf);
    %     this.Variables.z(this.Variables.z<10^-3) = 0;
    %     this.Variables.x(this.Variables.x<10^-3) = 0;
    tol_zero = 10^-4;
    this.Variables.x(this.Variables.x<tol_zero*max(this.Variables.x)) = 0;
    this.Variables.z(this.Variables.z<tol_zero*max(this.Variables.z)) = 0;
    options.Display = 'final';
    if ~this.checkFeasible([], options)
        warning('optimalFlowRate: the rounding of variables with small quantity will make the solution infeasible.');
    end
    this.setPathBandwidth;
    this.FlowTable.Rate = this.getFlowRate;
    this.VirtualLinks.Load = this.getLinkLoad;
    this.VirtualDataCenters.Load = this.getNodeLoad;
    profit = -fval;
end
end