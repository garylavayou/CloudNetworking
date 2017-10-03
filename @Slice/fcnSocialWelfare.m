%% Net social welfare of a slice
% The net social welfare is the total utility of the slice less the resources cost.
%
% The constant part of the objective fucntion can be ignored.

%% Function Declaration
%
%   [profit, grad] = fcnSocialWelfare(var_x, S)
%
function [profit, grad] = fcnSocialWelfare(var_x, S)
var_path = var_x(1:S.NumberPaths);
flow_rate = S.getFlowRate(var_path);

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
if isempty(S.weight)        % for network as a single slice.    TODO: remove this block
    weight = S.FlowTable.Weight;
else                        % for network as multiple slice.
    weight = S.weight*ones(S.NumberFlows, 1);
end
var_node = var_x((S.NumberPaths+1):end);
node_load = S.getNodeLoad(var_node);
link_load = S.getLinkLoad(var_path);
profit = -sum(weight.*fcnUtility(flow_rate)) + S.getResourceCost(node_load, link_load);

% If there is only one output argument, return the real profit (positive)
if nargout <= 1
    profit = -profit;
else
    link_price = S.link_unit_cost;
    node_price = S.node_unit_cost;
    grad = spalloc(length(var_x),1, S.NumberPaths+nnz(S.I_node_path)*S.NumberVNFs);
    for p = 1:S.NumberPaths
        i = S.path_owner(p);
        grad(p) = -weight(i)/(1+S.I_flow_path(i,:)*var_path) + ...
            dot(link_price, S.I_edge_path(:,p)); %#ok<SPRIX>
    end
    
    nz = S.NumberDataCenters*S.NumberPaths;
    z_index = S.NumberPaths+(1:nz);
    for f = 1:S.NumberVNFs
        % compatiable arithmetic operation
        grad(z_index) = node_price.*S.I_node_path; %#ok<SPRIX>
        z_index = z_index + nz;
    end
end