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
%   c) v: VNF instance capacity;    d) tx: auxiliary vairiables for x;
%   e) tz: auxiliary vairiables for z;
%   f) tx: auxiliary vairiables for v;
%   g) c: link capacity variables;
% |tx|,|tz|, and |tv| are auxiliary variables to transform L1 norm approximation of
% reconfiugration cost to linear cost with linear constraints.
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
%   # num_varx,num_varz,num_varv: number of variables.
%
% See also <fcnProfitReserveSlicing>, <hessSlicing>.
function [profit, gd] = fcnProfitReconfigureSlicing(vars, this, options)
if isfield(options, 'bCompact') && options.bCompact
    full_vars = zeros(options.num_orig_vars,1);
    full_vars(this.I_active_variables) = vars;
    vars = full_vars;
end
slice = this.hs;
Nsn = slice.NumberServiceNodes;
Nvf = slice.NumberVNFs;
Nl = slice.NumberLinks;
Np = slice.NumberPaths;
var_x_index = 1:options.num_varx;
idx_offset = options.num_varx+options.num_varz;
var_v_index = (1:options.num_varv) + idx_offset;
idx_offset = idx_offset + options.num_varv;
var_tx_index = (1:options.num_varx) + idx_offset;
idx_offset = idx_offset + options.num_varx;
var_tz_index = (1:options.num_varz) + idx_offset;
idx_offset = idx_offset + options.num_varz;
var_tv_index = (1:options.num_varv) + idx_offset;
idx_offset = idx_offset + options.num_varv;
var_c_index = (1:Nl) + idx_offset;

node_load = sum(reshape(vars(var_v_index), Nsn, Nvf),2);
link_load = vars(var_c_index);
var_path = vars(var_x_index);

flow_rate = this.getFlowRate(var_path);

switch options.PricingPolicy
    case {'quadratic-price', 'quadratic'}
        [link_payment,link_price_grad] = slice.fcnLinkPricing(this.prices.Link, link_load);
        [node_payment,node_price_grad] = slice.fcnNodePricing(this.prices.Node, node_load);
        profit = -slice.weight*sum(fcnUtility(flow_rate)) + link_payment + node_payment;
    case 'linear'
        profit = -slice.weight*sum(fcnUtility(flow_rate)) ...
            + dot(this.prices.Link, link_load) + dot(this.prices.Node, node_load);
    otherwise
        error('%s: invalid pricing policy', calledby);
end
profit = profit + dot(vars(var_tx_index), this.topts.x_reconfig_cost) + ...
    dot(vars(var_tz_index), this.topts.z_reconfig_cost) + ...
    dot(vars(var_tv_index), this.topts.vnf_reconfig_cost);
% When the 'bFinal' option is provided, return the real profit (positive).
if isfield(options, 'bFinal') && options.bFinal
    profit = -profit;
end

if nargout == 2
    nnz_grad = Np + options.num_varv + ...
        options.num_varx + options.num_varz + options.num_varv + Nl;
    gd = spalloc(length(vars), 1, nnz_grad);
    switch options.PricingPolicy
        case {'quadratic-price', 'quadratic'}
            for p = 1:Np
                i = slice.path_owner(p);
                gd(p) = -slice.weight/(1+this.I_flow_path(i,:)*var_path); %#ok<SPRIX>
            end
        case 'linear'
            for p = 1:Np
                i = slice.path_owner(p);
                gd(p) = -slice.weight/(1+this.I_flow_path(i,:)*var_path); %#ok<SPRIX>
            end
        otherwise
            error('%s: invalid pricing policy', calledby);
    end
    %% partial derviatives on w(v)
    gd(var_v_index) = repmat(node_price_grad, 1, Nvf);
    %% partial derivatives on c
    gd(var_c_index) = link_price_grad;
    %% partial derivatives on auxiliary variables tx,tz,tv
    gd(var_tx_index) = this.topts.x_reconfig_cost;
    gd(var_tz_index) = this.topts.z_reconfig_cost;
    gd(var_tv_index) = this.topts.vnf_reconfig_cost;
    
    if isfield(options, 'bCompact') && options.bCompact
        gd = gd(this.I_active_variables);
    end
end
end