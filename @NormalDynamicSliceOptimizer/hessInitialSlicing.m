%% Hessian matrix of for initial network slice dimensioning
% This function is used with <fcnProfitReserveSlicing>, when slices are initially added
% with explicit resource reservation.
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale.
%
% See also <NormalDynamicSliceOptimizer.fcnProfitReserveSlicing>,
% <NormalDynamicSliceOptimizer.hessSlicing>, <NormalSliceOptimizer.fcnHessianCC>.
function hs = hessInitialSlicing(vars, lambda, this, options) %#ok<INUSL>
% if options.bCompact
%     full_vars = sparse(options.num_orig_vars,1);
%     full_vars(this.I_active_variables) = vars;
%     vars = full_vars;
% end
if isempty(this.pardata.Weight)
	weight = this.pardata.FlowWeight;    % for single slice;
else
	weight = this.pardata.Weight*ones(this.pardata.NumberFlows, 1);  % for multiple slices
end
prbm = this.problem;
hs = sparse(numel(vars), numel(vars));
num_varxz = sum(prbm.num_vars(1:2));
num_vars_prim = sum(prbm.num_vars(1:3));
var_r_index = (num_varxz+1):num_vars_prim;
flow_rate = vars(var_r_index);
hs(var_r_index, var_r_index) = spdiag(weight./(1+flow_rate).^2);

if nargin >= 4 && isfield(options, 'PricingPolicy')
	switch options.PricingPolicy
		case {'quadratic-price', 'quadratic'}
			var_v_index = num_vars_prim + (1:prbm.num_vars(4));
			node_load = sum(reshape(vars(var_v_index), Nsn, Nvnf),2);
			var_c_index = num_vars - Nl + (1:Nl);
			link_load = vars(var_c_index);
			% Since we pricing the resources, we should know amount of resource occupied,
			% instead of the actual load.
			[~,~,lph] = this.fcnLinkPricing(this.prices.Link, link_load);
			[~,~,nph] = this.fcnNodePricing(this.prices.Node, node_load);
			% second derviatives of resource cost on (x,z) = 0;
			% second derviatives of resource cost on (c,w(v))
			hs(var_c_index,var_c_index) = diag(lph);
			hs(var_v_index,var_v_index) = block_diag(diag(nph), Nvnf);
		otherwise
			error('%s: invalid pricing policy', calledby);
	end
end

% if options.bCompact
%     hs = hs(this.I_active_variables, this.I_active_variables);
% end
end
