%% subproblemNetSocialWelfare
% subproblemNetSocialWelfare Summary of this function goes here
%   Detailed explanation goes here

%%
function [fval, node_load, link_load] = subproblemNetSocialWelfare( this, lambda )
%% Parameters
% |lambda.pf| is a matrix and lambda.npf is a 3-D array.
% |dg_pf|, |dg_p|, |dg_npf| is the inrement of gradient on lambda.
global DEBUG;
%% Set the feasible start point
x0 = zeros(num_vars,1);
x0(1:this.NumberPaths) = 1;
alpha_max = max(this.Parent.VNFTable.ProcessEfficiency(this.VNFList));
x0((this.NumberPaths+1):end) = alpha_max;
if ~this.checkFeasible(x0)
    error('error: infeasible start point.');
end
num_lcon_res = size(this.As_res,1);
bs = sparse(num_lcon_res,1);
lbs = sparse(num_vars,1);
% node_capacity = this.Parent.readNode('Capacity', this.Nodes.PhysicalNode);
% ubs = [inf*ones(this.NumberPaths,1);
%     repmat(node_capacity, this.NumberPaths*this.NumberVNFs, 1)];
% Algorithm option
minopts = optimoptions(@fmincon);
minopts.Algorithm = 'interior-point';
minopts.HessianFcn = @(x,lmd)SimpleSlice.fcnHessian(x,lmd,this);
minopts.SpecifyObjectiveGradient = true;
minopts.SpecifyConstraintGradient = false;
minopts.Display = 'notify';
[x, fval, exitflag] = fmincon(@(x)SimpleSlice.subproblemObjective(x, lambda, this), ...
    x0, this.As_res, bs, [], [], lbs, [], [], minopts);
% fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
if exitflag == 1 || exitflag == 2
%     x(x<10^-5) = 0;
    if ~this.checkFeasible(x)
        error('error: infeasible solution.');
    end
    this.x_path = x(1:this.NumberPaths);
    this.temp_vars.z = x((this.NumberPaths+1):end);
    %% when compute node load, z_npf corresonding to h_np = 0 has been set as zero.
    nz = this.NumberNodes*this.NumberPaths;
    z_index = 1:nz;
    for f = 1:this.NumberVNFs
        this.temp_vars.z(z_index) = this.I_dc_path(:).*this.temp_vars.z(z_index);
        z_index = z_index + nz;
    end
    %     this.x_path(this.x_path<10^-3) = 0;
    %     this.z_npf(this.z_npf<10^-3) = 0;
    this.x_path(this.x_path<(10^-4)*max(this.x_path)) = 0;
    this.z_npf(this.temp_vars.z<(10^-4)*max(this.temp_vars.z)) = 0;
    if ~this.checkFeasible([this.temp_vars.x; this.temp_vars.z])
        if ~isempty(DEBUG) && DEBUG
            warning('subproblemNetSocialWelfare: the rounding of variables %s', ...
                'with small quantity will make the solution infeasible.');
        end
    end

    % |x_path.*alpha_f'| is a compatible arithmetic operation.
    %     dg_pf = x_path.*alpha_f';
    %         z_np = reshape(z_npf(z_index), this.NumberNodes, this.NumberPaths);
    %         %% gradient of lambda(p,f)
    %         % |z_np = this.I_dc_path.*z_npf(index)|
    %         %
    %         % $$\alpha_f x_p^{(s)} -
    %         % \sum_{n\in\mathcal{N}^{(s)}}{h_{n,p}^{(s)}z_{n,p,f}^{(s)}}$$
    %         %         dg_pf(:,f) = dg_pf(:,f) - (sum(z_np, 1))';

    node_load = zeros(1, this.Parent.NumberNodes);
    link_load = zeros(1, this.Parent.NumberLinks);
    node_load(this.Nodes.PhysicalNode) = this.getNodeLoad(false, this.z_npf);
    link_load(this.Links.PhysicalLink) = this.getLinkLoad(false, this.x_path);   
%     if nargout > 3
%         g_p = -x_path;
%     end
%     if nargout > 4
%         g_npf = -reshape(z_npf, this.NumberNodes, this.NumberPaths, this.NumberVNFs);
%     end
else
    if ~isempty(DEBUG) && DEBUG
        warning('abnormal exit with flag %d.',exitflag);
    end
end
end

% function [x,z] = optimizeSingleFlow(this)
% 
% end
