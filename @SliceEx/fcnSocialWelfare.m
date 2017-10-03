%% Net social welfare of a slice
% The net social welfare is the total utility of the slice less the resources cost.
%
% The constant part of the objective fucntion can be ignored.

%% Function Declaration
%
%   [profit, grad] = fcnSocialWelfare(var_x, S, model)
%
% * *Model*: |Accurate|  or |Approximate|(default). 
function [profit, grad] = fcnSocialWelfare(var_x, S, model)
var_path = var_x(1:S.NumberPaths);
flow_rate = S.getFlowRate(var_path);

%% Calculate the cost of network and slices
% Called by _optimalFlowRate_ as objective function;
if isempty(S.weight)        % for network as a single slice.    TODO: remove this block
    weight = S.FlowTable.Weight;
else                        % for network as multiple slice.
    weight = S.weight*ones(S.NumberFlows, 1);
end
%%%
% get the slice cost according to the model(approximate or accurate).
if nargin <= 2
    warning('model is set as Approximate.');
    model = 'Approximate';
end
var_node = var_x((S.NumberPaths+1):end);
node_load = S.getNodeLoad(var_node);
link_load = S.getLinkLoad(var_path);
profit = -sum(weight.*fcnUtility(flow_rate)) ...
    + S.getResourceCost(node_load, link_load, model);
profit = profit - S.constant_profit;

% If there is only one output argument, return the real profit (positive)
if nargout <= 1
    profit = -profit;
else
    pn = S.Parent;
    link_price = S.link_unit_cost + pn.phis_l;
    node_price = S.node_unit_cost + pn.phis_n;
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