%% Link Pricing Function
%%
function [ payment, grad, pseudo_hess ] = fcnLinkPricing(this, link_price, link_load)
%% Quadratic pricing
% $b = \frac{1}{2} \rho_s$ and $a = \frac{b \cdot |S|}{V_n}$
% 
% the cost of a candidate path is usually much larger than that of the main path.
theta = 3;
% b = (2/(1+theta))*link_price;
b = link_price;
link_id = this.VirtualLinks.PhysicalLink;
aggr_link_usage = sum(this.Parent.link_usage,2);
%%%
% Since $\mathit{S}$ slices share the resource amount |V|, we consider the price when the
% load equals to 0 and $\frac{V}{\mathit{S}}$.
%
% $$ p(0) = b, p(\frac{V}{\mathit{S}}) = a\frac{V}{\mathit{S}}+b $$
%
% 
a = (theta-1)*b.*aggr_link_usage(link_id)./this.Parent.getLinkField('Capacity', link_id);

%%%
% When computing hession matrix, the first two output arguments are not needed.
% [payment, grad] = fcnLinkPricing(¡¤) is called in <fcnProfit>,
% [~, ~, hess] = fcnLinkPricing(¡¤) is called in <fcnHession>.
% See also <fcnNodePricing>
if nargout == 3  
    pseudo_hess = diag(a);
    payment = [];
    grad = [];
    return
end

payment = sum(0.5*a.*(link_load.^2) + b.*link_load);
if nargout == 2
    grad = a.*link_load + b;
end

end

%% Deprecated
% If a flow is served by multiple paths, the cost of the shortest path is usually much
% less than the cost of other alternative paths, since the alternative path has more hops,
% and its link cost is larger than that of the shortest path. Instead, the node cost of
% two path is similar since it is proportional to the data rate. To encourage using the
% alternative paths, the unit link cost should increase more slowly than unit node cost.
% Besides, the cost difference between two nodes is ususally small, so it is more flexible
% to migrate the processing demand to alternative nodes.

