%% Net social welfare of a slice
% The net social welfare is the total utility of the slice less the resources cost.
%
% The constant part of the objective fucntion can be ignored.

%% Function Declaration
%
%   [profit, grad] = fcnSocialWelfare(var_x, op, options)
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
function [profit, gd] = fcnSocialWelfare(vars, this, options)
if options.bCompact
    full_vars = sparse(options.num_orig_vars,1);
    full_vars(this.I_active_variables) = vars;
    vars = full_vars;
end
slice = this.hs;
Nvnf = slice.NumberVNFs;
num_varx = this.num_vars(1);
var_x = vars(1:num_varx);
flow_rate = this.getFlowRate(var_x);

%% Calculate the cost of network and slices
% When calculate the cost of the network as a single slice, we can use both
% _Slice.getSliceCost_ and _CloudNetwork.getNetworkCost_. Otherwise, use
% _Slice.getSliceCost_ to calculate each slice's cost. *For simplicity, we  only use
% _getSliceCost_*. 
%
% Called by _optimalFlowRate_ as objective function;
%%%
%     profit = -sum(weight.*fcnUtility(flow_rate)) + ...
%         S.Parent.totalCost(load);
if isempty(slice.Weight)        % for network as a single slice.    TODO: remove this block
    weight = slice.FlowTable.Weight;
else                        % for network as multiple slice.
    weight = slice.Weight*ones(slice.NumberFlows, 1);
end
var_node = vars(num_varx+(1:this.problem.num_vars(2)));
load.Node = slice.getNodeLoad(false, var_node);
load.Link = slice.getLinkLoad(false, var_x);  % equal to <getLinkCapacity>
profit = -sum(weight.*fcnUtility(flow_rate)) + slice.getResourceCost(load);

if options.bFinal
    profit = -profit;
else
	Np = slice.NumberPaths;
    link_price = slice.LinkCost;
    node_price = slice.NodeCost;
    gd = spalloc(length(vars),1, Np+nnz(this.I_dc_path)*Nvnf);
    for p = 1:Np
        i = slice.path_owner(p);
        gd(p) = -weight(i)/(1+this.I_flow_path(i,:)*var_x) + ...
            dot(link_price, this.I_edge_path(:,p)); %#ok<SPRIX>
    end
    
    nz = slice.NumberServiceNodes*Np;
    z_index = Np+(1:nz);
    for f = 1:Nvnf
        % compatiable arithmetic operation
        gd(z_index) = node_price.*this.I_dc_path; %#ok<SPRIX>
        z_index = z_index + nz;
    end
    
    if options.bCompact
        gd = gd(this.I_active_variables);
    end
end