%% Objective function for initial network slice dimensioning
% Calculate the objective function value and gradient.
%
% This objective function is used when the slices are initially added with explicit
% resource reservation (ERR), without consider reconfiguration cost.
% Under ERR, we consider aggregate resource reservation, i.e. the total demand of
% link/node resource should be less than a proportion of the total available link/node
% capacity in the slice. The capacity of a node equals to the sum of VNF instance capcity
% in that node.
% Due to resource reservation, the resource demand (x,z) might be less than allocated
% capacity (v,c).
%
% Due to resource reservation, we have optiomization variables including:
%   a) x: path variables;           b) z: VNF allocation variables;
%   c) r: flow rate;								d) v: VNF instance capacity;
%		e) c: link capacity variables;
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
%   # num_varx,num_varz,num_varv: number of variables.
%
% See also <fcnProfitReconfigureSlicing>, <hessInitialSlicing>.
function [profit, grad] = fcnProfitReserveSlicing(vars, this, options)
% if options.bCompact
%     full_vars = sparse(options.num_orig_vars,1);
%     full_vars(this.I_active_variables) = vars;
%     vars = full_vars;
% end
if ~isfield(options, 'unit')
	unit = 1;
else
	unit = options.unit;
end
Nsn = this.pardata.NumberServiceNodes;
Nvnf = this.pardata.NumberVNFs;
Nl = this.pardata.NumberLinks;

prbm = this.problem;
num_vars = sum(prbm.num_vars);
var_r_index = sum(prbm.num_vars(1:2))+(1:prbm.num_vars(3));
var_v_index = sum(prbm.num_vars(1:3)) + (1:prbm.num_vars(4));
var_c_index = num_vars-Nl+(1:Nl);
flow_rate = vars(var_r_index);
profit = -this.pardata.Weight*sum(fcnUtility(unit*flow_rate));
node_load = sum(reshape(vars(var_v_index), Nsn, Nvnf),2);
link_load = vars(var_c_index);
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
% When the 'bFinal' option is provided, return the real profit (max -f).
if options.bFinal
	profit = -profit;
end

if nargout == 2
	%% 
	% Now that we have the VNF capacity and link capacity as variables, we count the
	% resource payment with these variables. So we derive the gradient on these variables
	% instead of the flow-link (x) and flow-VNF (z) variables.
	% 
	% By default these variables is not compressed, as determined by constraint (3)(4).
	grad = sparse(numel(vars),1);
	%% partial derviatives on r
	grad(var_r_index) = -this.pardata.Weight./(1/unit+flow_rate);
	%% partial derviatives on w(v)
	grad(var_v_index) = repmat(node_price_grad, 1, Nvnf);
	%% partial derivatives on c
	grad(var_c_index) = link_price_grad;
end
end
