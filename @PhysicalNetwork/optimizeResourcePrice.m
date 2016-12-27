%% Optimize resource price
% first, the resource for each slice is the same.
%% TODO Resource Cost Model
function [output] = optimizeResourcePrice(this, init_price, options)
if nargin <= 2
   options.Display = 'final'; 
end
this.clearStates;

%% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_load = zeros(NN, NS);
link_load = zeros(NL, NS);
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
epsilon = this.unitStaticNodeCost;
phis_n = epsilon*this.delta*(NN-1)/this.totalNodeCapacity;
phis_l = epsilon*(1-this.delta)*(NN-1)/this.totalLinkCapacity;

%% Initial Price
if nargin >=2 && ~isempty(init_price)
    link_price = init_price.link;
    node_price = init_price.node;
else
    init_price.link = this.getLinkField('UnitCost') + phis_l;
    link_price = init_price.link;
    init_price.node = this.getNodeField('UnitCost') + phis_n;
    node_price = init_price.node;
end
beta_link = link_price.*(1/3).*link_capacity;
delta_link_price = beta_link./link_capacity;
beta_node = node_price.*(1/3).*node_capacity;
delta_node_price = beta_node./node_capacity;

number_iter = 1;
profit = zeros(NS+1,1);
new_profit = zeros(NS+1,1);
net_welfare = -inf;
alpha0 = 1;
a_link_violate = zeros(NL,1);
a_node_violate = zeros(NN,1);
while true
    alpha = alpha0;  %/number_iter;
    % announce the resource price and optimize each network slice
    for s = 1:NS
        this.slices{s}.VirtualLinks.Price = ...
            link_price(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.VirtualNodes.Price = ...
            node_price(this.slices{s}.VirtualNodes.PhysicalNode);        
        [new_profit(s), node_load(:,s), link_load(:,s)] ...
            = this.slices{s}.priceOptimalFlowRate([], options);
    end

    %% Compute the new resource price according to the resource consumption
    % First, the resource cost function is assumed to be linear, i.e.
    % $\phi_{l}(e) = \eta_{l}\times y_e$ and $\phi_{d}(n) = \eta_{n}\times
    % v_n$.
    % Second, the static node cost function is also assumed to be linear, i.e.
    % $N_{VNF}=(N-1)\theta+1$
    %
    % Under the above assumtions, the resource price is only related to the
    % residual capacity of each link and node. Given the log barrier function
    % $\beta\log{(1-\frac{y_e}{c_e})}$, the resource price is calculated as follows.
    %
    aggr_node_load = sum(node_load, 2);
    aggr_link_load = sum(link_load, 2);
    residual_link_capacity = link_capacity - aggr_link_load;
    b_violate = residual_link_capacity<=0;
    delta_link_price(b_violate) = delta_link_price(b_violate) * 2;
    new_link_price = init_price.link;
    a_link_violate = a_link_violate | b_violate;
    new_link_price(a_link_violate) = ...
        new_link_price(a_link_violate) + delta_link_price(a_link_violate);
    
    residual_node_capacity = node_capacity - aggr_node_load;
    b_violate = residual_node_capacity<=0;
    delta_node_price(b_violate) = delta_node_price(b_violate) * 2;
    new_node_price = init_price.node;
    a_node_violate = a_node_violate | b_violate;
    new_node_price(a_node_violate) = ...
        new_node_price(a_node_violate) + delta_node_price(a_node_violate);
    
%     delta_link_price(~b_violate) = beta_link(~b_violate)./residual_link_capacity(~b_violate);
%     new_link_price = link_uc + phis_l + delta_link_price;

%     delta_node_price(~b_violate) = beta_node(~b_violate)./residual_node_capacity(~b_violate);
%     new_node_price = node_uc + phis_n + delta_node_price;    

    % compute net social welfare
    new_net_welfare = 0;
    for s = 1:NS
        sl = this.slices{s};
        new_net_welfare = new_net_welfare + sl.weight*sum(log(sl.flow_rate));
    end
    new_net_welfare = new_net_welfare ...
        - this.getNodeCost(aggr_node_load) ...
        - this.getLinkCost(aggr_link_load) ...
        - (phis_n*sum(aggr_node_load)+phis_l*sum(aggr_link_load));
    % The stop condition
    %    (norm(new_link_price-link_price)/this.NumberLinks;
    stop_cond11 = norm(new_link_price-link_price)/norm(link_price);
    stop_cond12 = norm(new_node_price-node_price)/norm(node_price);
    %     norm(new_utility-utility)/NS < 10^-3;
    stop_cond2 = norm(new_profit-profit)/norm(profit);
    %     stop_cond3 = norm(new_net_welfare-net_welfare)/NS < 10^-3;
    stop_cond3 = norm(new_net_welfare-net_welfare)/norm(new_net_welfare);
    if strncmp(options.Display, 'iter', 4)
        fprintf('\tObjective value: %d.\n', new_net_welfare);
        fprintf('\t stop condition 1: (%G,%G).\n', stop_cond11,stop_cond12);
        fprintf('\t stop condition 2: %G.\n', stop_cond2);
        fprintf('\t stop condition 3: %G.\n', stop_cond3);
    end
    if (stop_cond11<10^-3&&stop_cond12<10^-3) 
        %|| stop_cond2<10^-4 %|| stop_cond3<10^-4
        break;
    else
        link_price = alpha*new_link_price+(1-alpha)*link_price;
        node_price = alpha*new_node_price+(1-alpha)*node_price;
        profit = new_profit;
        net_welfare = new_net_welfare;
    end    
    number_iter = number_iter + 1;
end
if strncmp(options.Display, 'notify', 6)
    fprintf('\tFirst stage objective value: %d.\n', new_net_welfare);
end

%% TODO
% 1. check why the solution is optimaler than the global optimal solution
% 2. refine the price adjusting procedure.
l = 0;
alpha = 1/2;
h = 1;
while true
    link_price = init_price.link;
    link_price(a_link_violate) = ...
        link_price(a_link_violate) + alpha*delta_link_price(a_link_violate);
    node_price = init_price.node;
    node_price(a_node_violate) = ...
        node_price(a_node_violate) + alpha*delta_node_price(a_node_violate);
    for s = 1:NS
        this.slices{s}.VirtualLinks.Price = ...
            link_price(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.VirtualNodes.Price = ...
            node_price(this.slices{s}.VirtualNodes.PhysicalNode);        
        [profit(s), node_load(:,s), link_load(:,s)] ...
            = this.slices{s}.priceOptimalFlowRate([], options);
    end
    aggr_node_load = sum(node_load, 2);
    aggr_link_load = sum(link_load, 2);
    residual_link_capacity = link_capacity - aggr_link_load;
    residual_node_capacity = node_capacity - aggr_node_load;
    if isempty(find(residual_link_capacity<0,1)) && isempty(find(residual_node_capacity<0,1))
        h = alpha;
        alpha = (l+alpha)/2;
        if h-l < 10^-2
            break;
        end
    else
        l = alpha;
        alpha = (h+alpha)/2;
    end
    number_iter = number_iter + 1;
end

%% Finalize substrate network
for s = 1:NS
    sl = this.slices{s};    
    sl.Variables.x = sl.x_path;
    sl.Variables.z = sl.z_npf;
    sl.VirtualNodes.Load = sl.getNodeLoad;
    sl.VirtualLinks.Load = sl.getLinkLoad;
    sl.FlowTable.Rate = sl.getFlowRate;
    sl.setPathBandwidth;
end
this.setNodeField('Load', aggr_node_load);
this.setLinkField('Load', aggr_link_load);
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);

%% output
output.link_price = link_price;
output.node_price = node_price;
output.link_load = aggr_link_load;
output.node_load = aggr_node_load;
output.welfare_approx = 0;
output.welfare_accurate = 0;
t = zeros(this.NumberSlices+1, 1);
output.profit = table(t, t, t, t, 'VariableNames',...
    {'ApproximatePercent', 'ApproximatePrice', 'AccuratePercent', 'AccuratePrice'});
clear t;
options.Model = 'Accurate';
output.flow_rate = [];
for s = 1:NS
    sl = this.slices{s};
    output.flow_rate = [output.flow_rate; sl.FlowTable.Rate];
    
    var_x = [sl.Variables.x; sl.Variables.z];
    output.profit.AccuratePrice(s) = - Slice.fcnProfit(var_x, sl);
    
    p = -Slice.fcnNetProfit(var_x, sl);
    output.welfare_approx = output.welfare_approx + p;
    output.profit.ApproximatePercent(s) = options.PercentFactor * p;
    p = -Slice.fcnNetProfit(var_x, sl, options);
    idx = sl.VirtualNodes.Load>0;
    % NOTE it is not accurate to calculate the static cost of each slice with this method.
    nid = sl.VirtualNodes.PhysicalNode;
    p = p - dot(sl.VirtualNodes.Load(idx)./this.getNodeField('Load', nid(idx)),...
        this.getNodeField('StaticCost', nid(idx)));
    output.profit.AccuratePercent(s) = options.PercentFactor * p;
    
    output.welfare_accurate = output.welfare_accurate...
        + sl.weight*sum(log(sl.FlowTable.Rate));
end
f = (1-options.PercentFactor)/options.PercentFactor;
output.profit.ApproximatePercent(end) = ...
    sum(output.profit.ApproximatePercent(1:(end-1))*f);
output.profit.AccuratePercent(end) = ...
    sum(output.profit.AccuratePercent(1:(end-1))*f);
output.profit.ApproximatePrice = output.profit.AccuratePrice;
output.profit.ApproximatePrice(end) = output.welfare_approx - ...
    sum(output.profit.ApproximatePrice(1:(end-1)));
output.welfare_accurate = output.welfare_accurate ...
    - this.getNetworkCost([],[], options.Model);
output.profit.AccuratePrice(end) = output.welfare_accurate - ...
    sum(output.profit.AccuratePrice(1:(end-1)));  
% output the optimization results
if strncmp(options.Display, 'notify', 6) || strncmp(options.Display, 'final', 5)
    fprintf('Optimization results:\n');
    fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
    fprintf('\tOptimal objective value: %d.\n', new_net_welfare);
    fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio);
end
end