%% Optimize resource price
% * *TODO* _refine the price adjustment algorithm._
% * *TODO* Resource Cost Model: linear, convex (quatratic)
%
% *Note 1*: When the approximation static cost is considered, and the network is lightly loaded, then
% approximaite static cost cannot be covered, thus the substrate network's profit is negative.
%
% *Note 2*: Slices with lower weight may obtains little resources, which results in low service rate 
% or low profit ratio. If there is policy about minimun rate or profit, then this slice may immediately 
% depart from the network.

%%
function [output, runtime] = optimizeResourcePrice(this, init_price, options)
if nargin <= 2
    options.Display = 'final';
end
this.clearStates;
if nargout == 2
    slice_runtime = 0;
    runtime.Serial = 0;
    runtime.Parallel = 0;
end

% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_load = zeros(NN, NS);
link_load = zeros(NL, NS);
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');

% Initial Price
if nargin >=2 && ~isempty(init_price)
    link_price = init_price.link;
    node_price = init_price.node;
else
    init_price.link = this.getLinkField('UnitCost') + this.phis_l;
    link_price = init_price.link;
    init_price.node = this.getNodeField('UnitCost') + this.phis_n;
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
    alpha = alpha0;  % /number_iter;
    % announce the resource price and optimize each network slice
    for s = 1:NS
        this.slices{s}.VirtualLinks.Price = ...
            link_price(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.VirtualNodes.Price = ...
            node_price(this.slices{s}.VirtualNodes.PhysicalNode);
        if nargout == 2
            tic;
        end
        [new_profit(s), node_load(:,s), link_load(:,s)] ...
            = this.slices{s}.priceOptimalFlowRate([], options);
        if nargout == 2
            t = toc;
            slice_runtime = max(slice_runtime, t);
            runtime.Serial = runtime.Serial + t;
        end
    end
    if nargout == 2
        runtime.Parallel = runtime.Parallel + slice_runtime;
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
    new_net_welfare = new_net_welfare - ...
        this.getNetworkCost(aggr_node_load, aggr_link_load);
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

l = 0;
alpha = 1/2;
h = 1;
while true && ~(nnz(a_link_violate)==0 && nnz(a_node_violate)==0)
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
        if nargout == 2
            tic;
        end
        [profit(s), node_load(:,s), link_load(:,s)] ...
            = this.slices{s}.priceOptimalFlowRate([], options);
        if nargout == 2
            t = toc;
            slice_runtime = max(slice_runtime, t);
            runtime.Serial = runtime.Serial + t;
        end
    end
    if nargout == 2
        runtime.Parallel = runtime.Parallel + slice_runtime;
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
% # The resource allocation variables, virtual node/link load, and flow rate of each
% slice.
% # After the optimization, each network slice has record the final prices.
% # Record the substrate network's node/link load, price.
for s = 1:NS
    sl = this.slices{s};
    sl.Variables.x = sl.x_path;
    sl.Variables.z = sl.z_npf;
    sl.setPathBandwidth;
    sl.VirtualNodes.Load = sl.getNodeLoad;
    sl.VirtualLinks.Load = sl.getLinkLoad;
    sl.FlowTable.Rate = sl.getFlowRate;
end
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);
this.setNodeField('Load', aggr_node_load);
this.setLinkField('Load', aggr_link_load);

%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The approximation of net social welfare, and the accurate net social welfare;
% # The net profit of each slice and the substrate network, including four results from
% different methods, i.e. |ApproximatePercent|, |ApproximatePrice|, |AccuratePercent|,
% |AccuratePrice|.
% # Flow rate of all flows in the network.
output.node_price = node_price;
output.link_price = link_price;
output.node_load = aggr_node_load;
output.link_load = aggr_link_load;
output = this.calculateOutput(output, options);

% output the optimization results
if strncmp(options.Display, 'notify', 6) || strncmp(options.Display, 'final', 5)
    fprintf('Optimization results:\n');
    fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
    fprintf('\tOptimal objective value: %d.\n', new_net_welfare);
    fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio);
end
end

%% Modifications
% * *Embedded Profit*: the announced link/node cost is greater than the true resource cost. 
% The surplus is treated as the substrate network's profit. With this feature (PhysicalNetwork.eta ~= 0)
% enabled, when network is lightly loaded, the substrate network still makes profit. Otherwise, the 
% substrate network's profit will be zero/negative if no resource tends to be over loaded. 
%
% * *Constant Profit*: this prevents negative profit of lower weighted slices, when resource is scarce.
% See also, PhysicalNetwork.constant_profit.