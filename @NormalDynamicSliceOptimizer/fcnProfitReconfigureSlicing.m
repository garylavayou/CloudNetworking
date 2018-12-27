%% Objective function for reconfiguration of network slice dimensioning
% Calculate the objective function value and gradient.
%
% This objective function is used when re-allocate resources for the slice with
% consideration on reconfiguration cost. Due to reconfiguration cost constraint or
% resource reservation, the resource demand might be less than the capacity to avoid
% reconfiguration.
%
% In the problem, we have the following variables:
%   a) x: path variables;           b) z: VNF allocation variables;
%   c) r: flow rate;								d) v: VNF instance capacity;
%		e) tx: auxiliary vairiables for x;
%   f) tz: auxiliary vairiables for z;
%   g) tx: auxiliary vairiables for v;
%   h) c: link capacity variables;
% |tx|,|tz|, and |tv| are auxiliary variables to transform L1 norm approximation of
% reconfiugration cost to linear cost with linear constraints.
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
%
% See also <fcnProfitReserveSlicing>, <hessSlicing>.
function [profit, grad] = fcnProfitReconfigureSlicing(vars, this, options)
% if options.bCompact
%     full_vars = sparse(options.num_orig_vars,1);
%     full_vars(this.I_active_variables) = vars;
%     vars = full_vars;
% end
if isfield(options, 'unit')
	unit = options.unit;
else
	unit = 1;
end
Nsn = this.pardata.NumberServiceNodes;
Nvnf = this.pardata.NumberVNFs;

prbm = this.problem;
idx_offset = sum(prbm.num_vars(1:2));
var_r_index = (1:prbm.num_vars(3))+idx_offset;
idx_offset = idx_offset + prbm.num_vars(3);
var_v_index = (1:this.num_vars(4)) + idx_offset;
idx_offset = idx_offset + prbm.num_vars(4);
var_tx_index = (1:prbm.num_vars(5)) + idx_offset; % the number of x and tx after compress is not the same
idx_offset = idx_offset + prbm.num_vars(5);
var_tz_index = (1:prbm.num_vars(6)) + idx_offset;
idx_offset = idx_offset + prbm.num_vars(6);
var_tv_index = (1:prbm.num_vars(7)) + idx_offset;
idx_offset = idx_offset + prbm.num_vars(7);
var_c_index = (1:prbm.num_vars(8)) + idx_offset;
flow_rate = vars(var_r_index);
profit = -this.pardata.Weight*sum(fcnUtility(unit*flow_rate));
link_load = vars(var_c_index);
node_load = sum(reshape(vars(var_v_index), Nsn, Nvnf),2);
switch options.PricingPolicy
	case {'quadratic-price', 'quadratic'}
		[link_payment,link_price_grad] = this.fcnLinkPricing(this.prices.Link, link_load, unit);
		[node_payment,node_price_grad] = this.fcnNodePricing(this.prices.Node, node_load, unit);
		profit = profit + link_payment + node_payment;
	case 'linear'
		profit = profit + dot(this.prices.Link, link_load*unit) + dot(this.prices.Node, node_load*unit);
	otherwise
		error('%s: invalid pricing policy', calledby);
end
profit = profit + dot(vars(var_tx_index), this.topts.x_reconfig_cost*unit) + ...
	dot(vars(var_tz_index), this.topts.z_reconfig_cost*unit) + ...
	dot(vars(var_tv_index), this.topts.vnf_reconfig_cost*unit);
% When the 'bFinal' option is provided, return the real profit (positive).
if options.bFinal
	profit = -profit;
end

if nargout == 2
	grad = sparse(numel(vars), 1);
	%% partial derviatives on r
	grad(var_r_index) = -this.pardata.Weight./(1/unit+flow_rate);
	%% partial derviatives on w(v)
	grad(var_v_index) = repmat(node_price_grad, 1, Nvnf);
	%% partial derivatives on c
	grad(var_c_index) = link_price_grad;
	%% partial derivatives on auxiliary variables tx,tz,tv
	grad(var_tx_index) = this.topts.x_reconfig_cost*unit;
	grad(var_tz_index) = this.topts.z_reconfig_cost*unit;
	grad(var_tv_index) = this.topts.vnf_reconfig_cost*unit;
	
% 	if options.bCompact
% 		grad = grad(this.I_active_variables);
% 	end
end
end