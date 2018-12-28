%% subproblemNetSocialWelfare
% subproblemNetSocialWelfare Summary of this function goes here
% FIXME: the flow is dynamically partitioned on to candidate paths, this does not work in
% terms of convergence.

%%
function [fval, node_load, link_load] = subproblemNetSocialWelfare2( this, lambda )
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
minopts.HessianFcn = @(x,la)SimpleSlice.fcnHessian(x,la,this);
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
    this.z_npf = x((this.NumberPaths+1):end);
    %% when compute node load, z_npf corresonding to h_np = 0 has been set as zero.
    nz = this.NumberNodes*this.NumberPaths;
    z_index = 1:nz;
    for f = 1:this.NumberVNFs
        this.z_npf(z_index) = this.I_dc_path(:).*this.z_npf(z_index);
        z_index = z_index + nz;
    end
    clear f;
    %     this.x_path(this.x_path<10^-3) = 0;
    %     this.z_npf(this.z_npf<10^-3) = 0;
    this.x_path(this.x_path<(10^-4)*max(this.x_path)) = 0;
    this.z_npf(this.z_npf<(10^-4)*max(this.z_npf)) = 0;
    if ~this.checkFeasible([this.x_path; this.z_npf])
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
    
    znpf = reshape(this.z_npf, this.NumberNodes, this.NumberPaths, this.NumberVNFs); % FOR DEBUG
    path_cost = this.getPathCost(lambda.e, lambda.n);
    pid_offset = 0;
    for i = 1:this.NumberFlows
        [spc, idp] = sort(path_cost{i}, 'ascend');
        p_reserved = length(idp);
        for r = 2:length(idp)
            if abs(spc(1)-spc(r))/((spc(1)+spc(r))/2) > 0.1
                p_reserved = r - 1;
                break;
            end
        end
        if p_reserved > 1
            % the flow is partitioned on the reserved paths
            pid = pid_offset + idp(1:p_reserved);
            propotion = (1./spc(1:p_reserved))/sum(1./spc(1:p_reserved));
            this.x_path(pid) = sum(this.x_path(pid)).*propotion;
            for p = pid'
                [nc, nid] = this.getPathNodeCost(p, lambda.n);
                [snc, idn] = sort(nc, 'ascend');
                n_reserved = length(idn);
                for r = 1:length(idn)
                    if abs(snc(1)-snc(r))/((snc(1)+snc(r))/2) > 0.05
                        n_reserved = r - 1;
                        break;
                    end
                end
                propotion = (1./snc(1:n_reserved))/sum(1./snc(1:n_reserved));
                for v = 1:this.NumberVNFs
                    alpha_f = this.Parent.VNFTable.ProcessEfficiency(this.VNFList(v));
                    node_load = this.x_path(p)*alpha_f;
                    z_id = nid(idn(1:n_reserved))+(p-1)*this.NumberNodes...
                        +(v-1)*this.NumberNodes*this.NumberPaths;
                    this.z_npf(z_id) = node_load.*propotion;
                end
            end
        end
        pid_offset = pid_offset + length(idp);
    end  
    fval = SimpleSlice.subproblemObjective([this.x_path; this.z_npf], lambda, this);
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
