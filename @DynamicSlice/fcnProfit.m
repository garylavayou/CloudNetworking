%%
% This method overrides <Slice.fcnProfit>, further providing the calculation of
% reconfiguration cost and updating the calculation method of resource consumption
% cost when resource reservation is specified.
%
% NOTE: since the objective function will be evaluated many times, to accelerate the
%   computation, we rewrite the method from the superclass, instead of calling the
%   superclass method.
function [profit, grad] = fcnProfit(vars, slice, options)
num_basic_vars = (options.num_varx+options.num_varz);        % number of (x,z)
num_orig_vars = num_basic_vars + options.num_varv;
num_vars = length(vars);
NP = slice.NumberPaths;
ND = slice.NumberDataCenters;
NV = slice.NumberVNFs;
NL = slice.NumberVirtualLinks;
var_path = vars(1:NP);
if num_vars == num_basic_vars
    % no reconfiguration cost constraint, node capacity equals to load.
    var_node = vars((NP+1):num_basic_vars);
    node_load = slice.getNodeLoad(var_node);
else
    % with reconfiguration cost constraint, node capacity equals to sum of VNF capacity.
    % VNF/node capacity lower-bound might be specified or not.
    var_vnf = reshape(vars(num_basic_vars+(1:options.num_varv)), ND, NV);
    node_load = sum(var_vnf,2);
end
if isempty(slice.lower_bounds)
    % without lower-bound, link capacity equals to sum of path variables.
    % there might be reconfiguration cost constraint on path variables.
    link_load = slice.getLinkLoad(var_path);
else
    % when the link capacity lower-bound is specified, we have the link capacity variable;
    link_load = vars(num_orig_vars*2+(1:NL));
end
flow_rate = slice.getFlowRate(var_path);

switch options.PricingPolicy
    case {'quadratic-price', 'quadratic'}
        [link_payment,link_price_grad] = slice.fcnLinkPricing(slice.prices.Link, link_load);
        [node_payment,node_price_grad] = slice.fcnNodePricing(slice.prices.Node, node_load);
        profit = -slice.weight*sum(fcnUtility(flow_rate)) + link_payment + node_payment;
    case 'linear'
        profit = -slice.weight*sum(fcnUtility(flow_rate)) ...
            + dot(slice.prices.Link, link_load) + dot(slice.prices.Node, node_load);
    otherwise
        error('%s: invalid pricing policy', calledby);
end

%%
% calculate reconfiguration cost and the gradient components on VNF instance
% variables and auxiliary variables.

% var_node = vars((S.NumberPaths+1):S.num_vars);
% |var_tx| and |var_tz| are auxilliary variables to transform L1 norm to
% linear constraints and objective part.
if isfield(options, 'num_varv')  % for _fastReconfigure2_, including |v|, and |tv|
    num_basic_vars = num_basic_vars + options.num_varv;
    var_tv = vars((num_basic_vars+slice.num_vars+1):num_basic_vars*2);
end
var_tx = vars(num_basic_vars+(1:NP));
var_tz = vars(num_basic_vars+NP+(1:slice.num_varz));

profit = profit + dot(var_tx, slice.topts.x_reconfig_cost) + ...
    dot(var_tz, slice.topts.z_reconfig_cost);
if isfield(options, 'num_varv') % for _fastConfigure2_
    profit = profit + dot(var_tv, slice.topts.vnf_reconfig_cost);
end

% When the 'bFinal' option is provided, return the real profit (positive).
if isfield(options, 'bFinal') && options.bFinal
    profit = -profit;
end
if nargout == 2
    % If no lower bound, the objective is the function of (x,z,v);
    % Otherwise, the objective is the function of (x,c,w,v)(w is function of v), 
    % since sum(w)<=c. 
    if isempty(slice.lower_bounds)
        grad = spalloc(num_vars, 1, NP+nnz(slice.I_node_path)*NV+num_basic_vars);
    else
        grad = spalloc(num_vars, 1, NP+ND*NV+num_basic_vars+NL);
    end
    for p = 1:slice.NumberPaths
        i = slice.path_owner(p);
        switch options.PricingPolicy
            case {'quadratic-price', 'quadratic'}
                grad(p) = -slice.weight/(1+slice.I_flow_path(i,:)*var_path); %#ok<SPRIX>
            case 'linear'
                grad(p) = -slice.weight/(1+slice.I_flow_path(i,:)*var_path); %#ok<SPRIX>
            otherwise
                error('%s: invalid pricing policy', calledby);
        end
    end
    if isempty(slice.lower_bounds)
        for p = 1:slice.NumberPaths
            switch options.PricingPolicy
                case 'quadratic-price'
                    grad(p) = grad(p) + dot(link_price_grad,slice.I_edge_path(:,p)); %#ok<SPRIX>
                otherwise
                    grad(p) = grad(p) + dot(slice.prices.Link,slice.I_edge_path(:,p)); %#ok<SPRIX>
            end
        end
        nz = (slice.NumberDataCenters*slice.NumberPaths);
        z_index = slice.NumberPaths+(1:nz);
        for f = 1:slice.NumberVNFs
            % compatible arithmetic operation: node_price is a row vector and S.I_node_path is
            % a matrix, and these two operants have the same number of rows.
            switch options.PricingPolicy
                case {'quadratic-price', 'quadratic'}
                    % |grad(z_index)| is a vector, and the right side is a matrix, the value
                    % of the matrix will be assigned to |grad(z_index)| column by column.
                    grad(z_index) = node_price_grad.*slice.I_node_path; %#ok<SPRIX>
                case 'linear'
                    grad(z_index) = slice.prices.Node.*slice.I_node_path; %#ok<SPRIX>
                otherwise
                    error('%s: invalid pricing policy', calledby);
            end
            z_index = z_index + nz;
        end
    else
        %% partial derviatives on w(v)
        v_index = slice.num_vars+(1:slice.num_varv);
        grad(v_index) = repmat(node_price_grad, 1, NV);
        %% partial derivatives on c
        c_index = (num_orig_vars*2+1):num_vars;
        grad(c_index) = link_price_grad;
    end
    grad(num_basic_vars+(1:NP)) = slice.topts.x_reconfig_cost;
    grad(num_basic_vars+NP+(1:slice.num_varz)) = slice.topts.z_reconfig_cost;
    if isfield(options, 'num_varv')
        var_offset = num_orig_vars*2;
        grad((var_offset-slice.num_varv+1):var_offset) = slice.topts.vnf_reconfig_cost;
    end
end
end

%{
function [profit, grad] = fcnProfit(vars, slice, options)
%             if isempty(vars)
%                 % |S.temp_vars| is set by <DynamicSlice.priceOptimalFlowRate>.
%                 vars = s.get_temp_vars(true);
%             end
%             if nargin <= 2
%                 options = struct;
%             end
var_xz = vars(1:slice.num_vars);
if nargout <= 1
    profit = fcnProfit@Slice(var_xz, slice, options);
else
    [profit, sub_grad] = fcnProfit@Slice(var_xz, slice, options);
end
%%
% calculate reconfiguration cost and the gradient components on VNF instance
% variables and auxiliary variables.
num_basic_vars = slice.num_vars;
num_vars = length(vars);
NP = slice.NumberPaths;
% var_node = vars((S.NumberPaths+1):S.num_vars);
% |var_tx| and |var_tz| are auxilliary variables to transform L1 norm to
% linear constraints and objective part.
if isfield(options, 'num_varv')  % for _fastReconfigure2_, including |v|, and |tv|
    num_basic_vars = num_basic_vars + options.num_varv;
    var_tv = vars((num_vars-options.num_varv+1):end);
end
var_tx = vars(num_basic_vars+(1:NP));
var_tz = vars(num_basic_vars+NP+(1:slice.num_varz));

profit = profit + dot(var_tx, slice.topts.x_reconfig_cost) + ...
    dot(var_tz, slice.topts.z_reconfig_cost);
if isfield(options, 'num_varv') % for _fastConfigure2_
    profit = profit + dot(var_tv, slice.topts.vnf_reconfig_cost);
end

% If there is only one output argument, return the real profit (positive)
if nargout <= 1
    profit = -profit;
else
    grad = spalloc(num_vars, 1, NP+num_vars/2);
    grad(1:slice.num_vars) = sub_grad;
    grad(num_basic_vars+(1:NP)) = slice.topts.x_reconfig_cost;
    grad(num_basic_vars+NP+(1:slice.num_varz)) = slice.topts.z_reconfig_cost;
    if isfield(options, 'num_varv')
        grad((num_vars-slice.num_varv+1):end) = slice.topts.vnf_reconfig_cost;
    end
end
end
%}