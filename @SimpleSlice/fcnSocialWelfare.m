%% Net social welfare of a slice
% The net social welfare is the total utility of the slice less the resources cost.
%
% The constant part of the objective fucntion can be ignored.

%% Function Declaration
%
%   [profit, grad] = fcnSocialWelfare(var_x, S, options)
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
function [profit, gd] = fcnSocialWelfare(vars, slice, options)
if isfield(options, 'bCompact') && options.bCompact
    full_vars = zeros(options.num_orig_vars,1);
    full_vars(slice.I_active_variables) = vars;
    vars = full_vars;
end

var_x = vars(1:slice.NumberPaths);
flow_rate = slice.getFlowRate(var_x);

%% Calculate the cost of network and slices
% When calculate the cost of the network as a single slice, we can use both
% _Slice.getSliceCost_ and _CloudNetwork.getNetworkCost_. Otherwise, use
% _Slice.getSliceCost_ to calculate each slice's cost. *For simplicity, we  only use
% _getSliceCost_*. 
%
% Called by _optimalFlowRate_ as objective function;
%%%
%     profit = -sum(weight.*fcnUtility(flow_rate)) + ...
%         S.Parent.getNetworkCost(node_load, link_load);
if isempty(slice.weight)        % for network as a single slice.    TODO: remove this block
    weight = slice.FlowTable.Weight;
else                        % for network as multiple slice.
    weight = slice.weight*ones(slice.NumberFlows, 1);
end
var_node = vars((slice.NumberPaths+1):slice.num_vars);
node_load = slice.getNodeLoad(var_node);
link_load = slice.getLinkLoad(var_x);  % equal to <getLinkCapacity>
profit = -sum(weight.*fcnUtility(flow_rate)) + slice.getResourceCost(node_load, link_load);

if isfield(options, 'bFinal') && options.bFinal
    profit = -profit;
else
    link_price = slice.link_unit_cost;
    node_price = slice.node_unit_cost;
    gd = spalloc(length(vars),1, slice.NumberPaths+nnz(slice.I_dc_path)*slice.NumberVNFs);
    for p = 1:slice.NumberPaths
        i = slice.path_owner(p);
        gd(p) = -weight(i)/(1+slice.I_flow_path(i,:)*var_x) + ...
            dot(link_price, slice.I_edge_path(:,p)); %#ok<SPRIX>
    end
    
    nz = slice.NumberDataCenters*slice.NumberPaths;
    z_index = slice.NumberPaths+(1:nz);
    for f = 1:slice.NumberVNFs
        % compatiable arithmetic operation
        gd(z_index) = node_price.*slice.I_dc_path; %#ok<SPRIX>
        z_index = z_index + nz;
    end
    
    if isfield(options, 'bCompact') && options.bCompact
        gd = gd(slice.I_active_variables);
    end
end