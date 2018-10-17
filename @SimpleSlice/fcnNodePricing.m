%% Node Pricing Function
% See <fcnLinkPricing>. 
%%
function [ payment, grad, pseudo_hess ] = fcnNodePricing(this, node_price, node_load)
%% Quadratic pricing
% $b = \frac{2}{3} \rho_s$ and $a = \frac{b \cdot |S|}{V_n}$
theta = 3;
b = node_price;
node_id = this.getDCNI;
aggr_node_usage = this.Parent.AggregateNodeUsage;
a = (theta-1)*b.*aggr_node_usage(node_id)./this.Parent.getNodeField('Capacity', node_id);

if nargout == 3
    pseudo_hess = a;
    payment = [];
    grad = [];
    return;
end

payment = sum(0.5*a.*(node_load.^2) + b.*node_load);
if nargout == 2
    grad = a.*node_load + b;
end
end

