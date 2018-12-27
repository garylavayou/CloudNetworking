%% Node Pricing Function
% See <fcnLinkPricing>. 
%%
function [ payment, grad, pseudo_hess ] = fcnNodePricing(this, node_price, node_load, unit)
%% Quadratic pricing
% $b = \frac{2}{3} \rho_s$ and $a = \frac{b \cdot |S|}{V_n}$
theta = 3;
if nargin <= 3
	unit = 1;
end
if isempty(this.hs)  
	% In parallel mode, the host slice is unavailable, the data is stored in |pardata|.
	aggr_node_usage = this.pardata.AggregateNodeUsage;
	phy_node_capacity = this.pardata.PhysicalNodeCapacity;
else
	slice = this.hs;
	aggr_node_usage = slice.Parent.AggregateNodeUsage(slice.getSNPI());
	phy_node_capacity = slice.Parent.readDataCenter('Capacity',slice.getDCPI());
end
b = node_price * unit;
a = unit * (theta-1)*b.*aggr_node_usage./phy_node_capacity;

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

