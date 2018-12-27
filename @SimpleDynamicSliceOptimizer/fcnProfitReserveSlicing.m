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
%   c) v: VNF instance capacity;    d) c: link capacity variables;
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
%   # num_varx,num_varz,num_varv: number of variables.
%
% See also <fcnProfitReconfigureSlicing>, <hessInitialSlicing>.
function [profit, gd] = fcnProfitReserveSlicing(vars, this, options)
if options.bCompact
    full_vars = sparse(options.num_orig_vars,1);
    full_vars(this.I_active_variables) = vars;
    vars = full_vars;
end
slice = this.hs;
Nsn = slice.NumberServiceNodes;
Nvnf = slice.NumberVNFs;
Nl = slice.NumberLinks;
Np = slice.NumberPaths;
var_x_index = 1:this.num_vars(1);		%% USE the original variables
idx_offset = sum(this.num_vars(1:2));
var_v_index = (1:this.num_vars(3)) + idx_offset;
idx_offset = idx_offset + this.num_vars(3);
var_c_index = (1:Nl) + idx_offset;

node_load = sum(reshape(vars(var_v_index), Nsn, Nvnf),2);
link_load = vars(var_c_index);
var_path = vars(var_x_index);

flow_rate = this.getFlowRate(var_path);

switch options.PricingPolicy
    case {'quadratic-price', 'quadratic'}
        [link_payment,link_price_grad] = this.fcnLinkPricing(this.prices.Link, link_load);
        [node_payment,node_price_grad] = this.fcnNodePricing(this.prices.Node, node_load);
        profit = -slice.Weight*sum(fcnUtility(flow_rate)) + link_payment + node_payment;
    case 'linear'
        profit = -slice.Weight*sum(fcnUtility(flow_rate)) ...
            + dot(this.prices.Link, link_load) + dot(this.prices.Node, node_load);
    otherwise
        error('%s: invalid pricing policy', calledby);
end
% When the 'bFinal' option is provided, return the real profit (max -f).
if options.bFinal
    profit = -profit;
end

if nargout == 2
    nnz_grad = Np+this.num_vars(3)+Nl;
    gd = spalloc(length(vars), 1, nnz_grad);
    switch options.PricingPolicy
        case {'quadratic-price', 'quadratic'}
            for p = 1:Np
                i = slice.path_owner(p);
                gd(p) = -slice.Weight/(1+this.I_flow_path(i,:)*var_path); %#ok<SPRIX>
            end
        case 'linear'
            for p = 1:Np
                i = slice.path_owner(p);
                gd(p) = -slice.Weight/(1+this.I_flow_path(i,:)*var_path); %#ok<SPRIX>
            end
        otherwise
            error('%s: invalid pricing policy', calledby);
    end
    %% partial derviatives on w(v)
    gd(var_v_index) = repmat(node_price_grad, 1, Nvnf);
    %% partial derivatives on c
    gd(var_c_index) = link_price_grad;

    if options.bCompact
        gd = gd(this.I_active_variables);
    end
end
end
