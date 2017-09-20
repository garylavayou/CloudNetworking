%% priceOptimalFlowRate
% priceOptimalFlowRate return the optimal flow rate for each flow in the slice, and the
% net profit of the slice. 
% TODO: by setting the capacity as |inf|, this method is equivalent to _optimalFlowRate_.
%% Function Prototype
%   [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
% |options|: if price is not provied in |options|, then this method will use the price
% stored in the slice.
function [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
if nargin == 2 && ~isempty(x0)
    this.x0 = x0;
else%if isempty(this.x0)
    this.x0 = zeros(this.num_vars,1);
    this.x0(1:this.NumberPaths) = 1;
    alpha_max = max(this.Parent.VNFTable.ProcessEfficiency(this.VNFList));
    this.x0((this.NumberPaths+1):end) = alpha_max;
end
if ~this.checkFeasible(this.x0)
    error('error: infeasible start point.');
end
if nargin <= 2 || ~isfield(options, 'Display')
   options.Display = 'final'; 
else
    if strfind(options.Display, 'final')
        options.Display = 'final';
    elseif strfind(options.Display, 'notify')
        options.Display = 'notify';
    end    
end
    
bs = sparse(this.num_lcon_res,1);
lbs = sparse(this.num_vars,1);
%% Set the optimization options
% * *Algorithm* : since the problem contains linear constriants and bound
% constraints, then |trust-region-reflective| method is not applicable. Hence,
% we choose the |inteirior point| method. As a result the Hessian matix should
% be computed separately.
% is directly returned from the objective function as the second derivatives.
% * *HessianFcn* : we comupte Hession using the objective function.
% Therefore, this option is set as |'objective'|.
% * *SpecifyObjectiveGradient* : the gradient can be directly computed from
% the objective function, so this option is set to |true|.
% * *SpecifyConstraintGradient* : since this problem does not contain nonlinear
% constraint, this option is set to |false|.
% * *Display information* : use |'iter'| to display itration information for
% debug. use |'notify-detailed'| to only display exception message.
fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.HessianFcn = @(x,lambda)Slice.fcnHessian(x, lambda, this, options);
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.SpecifyConstraintGradient = false;
fmincon_opt.Display = options.Display;   %'notify-detailed'; %'iter';
[x, fval, exitflag] = fmincon(@(x)Slice.fcnProfit(x, this, options), ...
    this.x0, this.As_res, bs, [], [], lbs, [], [], fmincon_opt);
if strncmp(options.Display, 'final', 5) || strncmp(options.Display, 'iter', 4) ||...
        strncmp(options.Display, 'notify', 6)
    fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
end
if exitflag == 0
    warning('reaching maximum number of iterations.');
elseif exitflag < 0
    error('abnormal exit with flag %d.',exitflag);
elseif exitflag ~= 1
    warning('local optimal solution found.');
end

%% Output solution
if ~this.checkFeasible(x, options)
    error('error: infeasible solution.');
end
this.x_path = x(1:this.NumberPaths);
this.z_npf = x((this.NumberPaths+1):end);
%%%
% When compute node load, z_npf corresonding to h_np = 0 has been set as zero.
nz = this.NumberDataCenters*this.NumberPaths;
z_index = 1:nz;
for f = 1:this.NumberVNFs
    this.z_npf(z_index) = this.I_node_path(:).*this.z_npf(z_index);
    z_index = z_index + nz;
end
% this.x_path(this.x_path<10^-3) = 0;
% this.z_npf(this.z_npf<10^-3) = 0;
tol_zero = 10^-4;
this.x_path(this.x_path<tol_zero*max(this.x_path)) = 0;
this.z_npf(this.z_npf<tol_zero*max(this.z_npf)) = 0;
options.Display = 'final';
if ~this.checkFeasible([this.x_path; this.z_npf], options)
    warning('the rounding of variables with small quantity will make the solution infeasible.');
end
this.x0 = [this.x_path; this.z_npf];

%     this.x_path = sparse();
%     z = reshape(this.x0((this.NumberPaths+1):end), ...
%         this.NumberVirtualNodes, this.NumberPaths, this.NumberVNFs);
%     this.z_node_path_vnf = cell(this.NumberVNFs,1);
%     for f = 1:this.NumberVNFs
%         this.z_node_path_vnf{f} = sparse(z(:,:,f));
%     end
node_load = zeros(this.Parent.NumberDataCenters,1);
link_load = zeros(this.Parent.NumberLinks,1);
data_center_id = this.getDCPI;
node_load(data_center_id) = this.getNodeLoad(this.z_npf);
link_load(this.VirtualLinks.PhysicalLink) = this.getLinkLoad(this.x_path);
this.flow_rate = this.getFlowRate(this.x_path);
net_profit = -fval;
%%%
% FOR DEBUG
% this.setPathBandwidth(this.x_path);
end