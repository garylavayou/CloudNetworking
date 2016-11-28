function [rate, net_profit] = optimalFlowRate(this, x0)
% optimalFlowRate return the optimal flow rate for each flow in the slice, and the net
%       profit of the slice. This function also 
if nargin == 2
    this.x0 = x0;
else%if isempty(this.x0)
    this.x0 = zeros(this.num_vars,1);
    this.x0(1:this.NumberPaths) = 1;
    max_alpha_f = max(this.Parent.VNFTable.ProcessEfficiency(this.VNFList));
    this.x0((this.NumberPaths+1):end) = 1*this.NumberPaths*max_alpha_f;
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
fmincon_opt.HessianFcn = @(x,lambda)Slice.fcnHessian(x,lambda,this);
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.SpecifyConstraintGradient = false;
fmincon_opt.Display = 'notify-detailed';
[x, fval, exitflag] = fmincon(@(x)Slice.fcnUtility(x,this), ...
    this.x0, this.As_res, bs, [], [], lbs, [], [], fmincon_opt);
fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
if exitflag == 1 || exitflag == 2
    %% Output solution
    this.x0 = x;
    index0 = (1:numel(this.I_node_path));
    index = this.NumberPaths+index0;
    % mask those components which is multiplied by h(n,p)=0.
    for f = 1:this.NumberVNFs
        this.x0(index) = x(index)'.*this.I_node_path(index0);
        index = index + numel(this.I_node_path);
    end
    this.x0(this.x0 < 1) = 0;
%     this.x_path = sparse();
%     z = reshape(this.x0((this.NumberPaths+1):end), ...
%         this.NumberVirtualNodes, this.NumberPaths, this.NumberVNFs);
%     this.z_node_path_vnf = cell(this.NumberVNFs,1);
%     for f = 1:this.NumberVNFs
%         this.z_node_path_vnf{f} = sparse(z(:,:,f));
%     end
    this.link_load = this.getLinkLoad(this.x0(1:this.NumberPaths));
    this.node_load = this.getNodeLoad(this.x0((this.NumberPaths+1):end));
    this.flow_rate = this.getFlowRate(this.x0(1:this.NumberPaths));
    rate = this.flow_rate;
    net_profit = -fval;
    % FOR DEBUG
    this.setPathBandwidth(this.x0(1:this.NumberPaths));
else
    rate = [];
    net_profit = [];
end
if exitflag ~= 1
    warning('abnormal exit with flag %d.',exitflag);
end
end