%% Net social welfare of a slice
% The net social welfare is the total utility of the slice less the resources cost.
%
% The constant part of the objective fucntion can be ignored.

%% Function Declaration
%
%   [profit, grad] = fcnSocialWelfare(var_x, S)
%
function [profit, grad] = fcnSocialWelfare(vars, slice)
var_path = vars(1:slice.NumberPaths);
flow_rate = slice.getFlowRate(var_path);

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
var_node = vars((slice.NumberPaths+1):end);
node_load = slice.getNodeLoad(var_node);
link_load = slice.getLinkLoad(var_path);  % equal to <getLinkCapacity>
profit = -sum(weight.*fcnUtility(flow_rate)) + slice.getResourceCost(node_load, link_load);

% If there is only one output argument, return the real profit (positive)
if nargout <= 1
    profit = -profit;
else
    link_price = slice.link_unit_cost;
    node_price = slice.node_unit_cost;
    grad = spalloc(length(vars),1, slice.NumberPaths+nnz(slice.I_node_path)*slice.NumberVNFs);
    for p = 1:slice.NumberPaths
        i = slice.path_owner(p);
        grad(p) = -weight(i)/(1+slice.I_flow_path(i,:)*var_path) + ...
            dot(link_price, slice.I_edge_path(:,p)); %#ok<SPRIX>
    end
    
    nz = slice.NumberDataCenters*slice.NumberPaths;
    z_index = slice.NumberPaths+(1:nz);
    for f = 1:slice.NumberVNFs
        % compatiable arithmetic operation
        grad(z_index) = node_price.*slice.I_node_path; %#ok<SPRIX>
        z_index = z_index + nz;
    end
end