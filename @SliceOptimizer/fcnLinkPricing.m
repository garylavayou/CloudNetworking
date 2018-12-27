%% Link Pricing Function
%%
function [ payment, grad, pseudo_hess ] = fcnLinkPricing(this, link_price, link_load, unit)
%% Quadratic pricing
% $b = \frac{1}{2} \rho_s$ and $a = \frac{b \cdot |S|}{V_n}$
% 
% the cost of a candidate path is usually much larger than that of the main path.
theta = 3;
% b = (2/(1+theta))*link_price;
if nargin <= 3
	unit = 1;
end
b = link_price * unit;
if isempty(this.hs)
	% In parallel mode, the host slice is unavailable, the data is stored in |pardata|.
	aggr_link_usage = this.pardata.AggregateLinkUsage;
	phy_link_capacity = this.pardata.PhysicalLinkCapacity;
else
	slice = this.hs;
	aggr_link_usage = slice.Parent.AggregateLinkUsage(slice.Links.PhysicalLink);
	phy_link_capacity = slice.Parent.readLink('Capacity', slice.Links.PhysicalLink);
end
%%%
% Since $\mathit{S}$ slices share the resource amount |V|, we consider the price when the
% load equals to 0 and $\frac{V}{\mathit{S}}$.
%
% $$ p(0) = b, p(\frac{V}{\mathit{S}}) = a\frac{V}{\mathit{S}}+b $$
%
% NOTE: we may compute prices for only a part of slices, while we do not use the residual
% capacity and the involved slices, since this might lead to very low prices (when
% residual resource is plentful, and involved slices is few). 
a = unit * (theta-1)*b.*aggr_link_usage./phy_link_capacity;

%%%
% When computing hession matrix, the first two output arguments are not needed.
% [payment, grad] = fcnLinkPricing(ρ) is called in <fcnProfit>,
% [~, ~, hess] = fcnLinkPricing(ρ) is called in <fcnHession>.
% See also <fcnNodePricing>
if nargout == 3  
    pseudo_hess = a;
    payment = [];
    grad = [];
    return
end

payment = sum(0.5*a.*(link_load.^2) + b.*link_load);
if nargout == 2
    grad = a.*link_load + b;
end

end

%% Un-determined
% If a flow is served by multiple paths, the link cost of the shortest path is usually
% notably less than that of the alternative paths, since the alternative path has more
% hops (links). Instead, the node cost of two path is similar since the pricessing demand
% is proportional to the data rate. 
%
% To encourage using the alternative paths, the link prices should increase faster than
% node prices. Besides, the cost difference between two nodes is ususally small, so it is
% more easily to migrate the processing demand to alternative nodes.

