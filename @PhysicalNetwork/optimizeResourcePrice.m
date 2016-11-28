%% Optimize resource price
% first, the resource for each slice is the same.
%% TODO Resource Cost Model
function optimizeResourcePrice(this, init_link_price, init_node_price)
alpha = 0.5;
node_uc = this.getNodeField('UnitCost');
link_uc = this.getLinkField('UnitCost');
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
delta = this.Delta;
epsilon = this.staticNodeCost;
N = this.NumberNodes;
phis_n = epsilon*delta*(N-1)/this.totalNodeCapacity;
phis_l = epsilon*(1-delta)*(N-1)/this.totalLinkCapacity;
if nargin >=2
    link_price = init_link_price;
else
    %     link_price = ones(this.NumberLinks,1);
    link_price = link_uc + phis_l;
    beta_link = link_price.*(1/3).*link_capacity;
    this.setLinkField('Beta', beta_link);
end
if nargin >= 3
    node_price = init_node_price;
else
    %     node_price = ones(this.NumberNodes,1);
    node_price = node_uc + phis_n;
    beta_node = node_price.*(1/3).*node_capacity;
    this.setNodeField('Beta', beta_node);
end
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);

number_iter = 1;
utility = zeros(this.NumberSlices,1);
new_utility = zeros(this.NumberSlices,1);
rate = cell(this.NumberSlices,1);
net_welfare = -inf;
while true
%     alpha = 1/number_iter;
    % announce the resource price and optimize each network slice
    for i = 1:this.NumberSlices
        this.slices{i}.VirtualLinks.Price = ...
            link_price(this.slices{i}.VirtualLinks.PhysicalLink);
        this.slices{i}.VirtualNodes.Price = ...
            node_price(this.slices{i}.VirtualNodes.PhysicalNode);        
        [rate{i}, new_utility(i)] = this.slices{i}.optimalFlowRate;
    end

    % use the optimization results of all network slices to update the physical links and
    % nodes' load 
    link_load = zeros(this.NumberLinks,1);
    node_load = zeros(this.NumberNodes,1);
    for i = 1:this.NumberSlices
        eid = this.slices{i}.VirtualLinks.PhysicalLink;
        link_load(eid) = link_load(eid) + this.slices{i}.link_load;
        nid = this.slices{i}.VirtualNodes.PhysicalNode;
        node_load(nid) = node_load(nid) + this.slices{i}.node_load;
    end

    %% Compute the new resource price according to the resource consumption
    % First, the resource cost function is assumed to be linear, i.e.
    % $\phi_{l}(e) = \eta_{l}\times y_e$ and $\phi_{d}(n) = \eta_{n}\times
    % v_n$.
    %
    % Second, the static node cost function is also assumed to be linear, i.e.
    % $N_{VNF}=(N-1)\theta+1$
    %
    % Under the above assumtions, the resource price is only related to the
    % residual capacity of each link and node.
    residual_link_capacity = link_capacity-link_load;
    b_violate = residual_link_capacity<=0;
    residual_link_capacity(b_violate) = eps;
    new_link_price = link_uc + phis_l + this.getLinkField('Beta')./residual_link_capacity;
    new_link_price(b_violate) = link_price(b_violate)*2;
    residual_node_capacity = node_capacity-node_load;
    b_violate = residual_node_capacity<=0;
    residual_node_capacity(b_violate) = eps;
    new_node_price = node_uc + phis_n + this.getNodeField('Beta')./residual_node_capacity;    
    new_node_price(b_violate) = node_price(b_violate)*2;
    % compute net social welfare
    total_profit = 0;
    for i = 1:this.NumberSlices
        total_profit = total_profit + this.slices{i}.weight*sum(log(rate{i}));
    end
    new_net_welfare = total_profit - this.getNodeCost(node_load) ...
        - this.getLinkCost(link_load) - epsilon*...
        ((N-1)*this.networkUtilization(node_load,link_load)+1);
    % The stop condition
    %    (norm(new_link_price-link_price)/this.NumberLinks;
    stop_cond11 = norm(new_link_price-link_price)/norm(new_link_price);
    stop_cond12 = norm(new_node_price-node_price)/norm(new_node_price);
    %     norm(new_utility-utility)/this.NumberSlices < 10^-3;
    stop_cond2 = norm(new_utility-utility)/norm(utility);
    %     stop_cond3 = norm(new_net_welfare-net_welfare)/this.NumberSlices < 10^-3;
    stop_cond3 = norm(new_net_welfare-net_welfare)/norm(new_net_welfare);
    fprintf('\t stop condition 1: (%G,%G).\n', stop_cond11,stop_cond12);
    fprintf('\t stop condition 2: %G.\n', stop_cond2);
    fprintf('\t stop condition 3: %G.\n', stop_cond3);
    if (stop_cond11<10^-3&&stop_cond12<10^-3) || stop_cond2<10^-4 || stop_cond3<10^-4
        break;
    else
        link_price = alpha*new_link_price+(1-alpha)*link_price;
        node_price = alpha*new_node_price+(1-alpha)*node_price;
        utility = new_utility;
        net_welfare = new_net_welfare;
    end    
    number_iter = number_iter + 1;
end
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);
this.setNodeField('Load', node_load);
this.setLinkField('Load', link_load);

% output the optimization results
disp('Optimization results:');
fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
fprintf('\tNormalized network utilization is %G.\n', this.networkUtilization);
end