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
function [output, runtime] = optimizeResourcePrice0(this, init_price)
global DEBUG;
options = getstructfields(this.options, {'PricingFactor', 'PercentFactor'});
options.Stage ='temp';

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
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');

% Initial Price
if nargin >=2 && ~isempty(init_price)
    prices.Link = init_price.Link;
    prices.Node = init_price.Node;
else
    init_price.Link = this.getLinkCost;
    prices.Link = init_price.Link;
    init_price.Node = this.getNodeCost;
    prices.Node = init_price.Node;
end
beta_link = prices.Link.*(1/3).*link_capacity;
delta_link_price = beta_link./link_capacity;
beta_node = prices.Node.*(1/3).*node_capacity;
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
        this.slices{s}.prices.Link= ...
            prices.Link(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.prices.Node = ...
            prices.Node(this.slices{s}.VirtualNodes.PhysicalNode);
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
        this.slices{s}.prices.Link = [];
        this.slices{s}.prices.Node = [];
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
    load = this.getNetworkLoad([], options);
    residual_link_capacity = link_capacity - load.Link;
    b_violate = residual_link_capacity<=0;
    delta_link_price(b_violate) = delta_link_price(b_violate) * 2;
    new_link_price = init_price.Link;
    a_link_violate = a_link_violate | b_violate;
    new_link_price(a_link_violate) = ...
        new_link_price(a_link_violate) + delta_link_price(a_link_violate);
    
    residual_node_capacity = node_capacity - load.Node;
    b_violate = residual_node_capacity<=0;
    delta_node_price(b_violate) = delta_node_price(b_violate) * 2;
    new_node_price = init_price.Node;
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
        new_net_welfare = new_net_welfare + sl.weight*sum(fcnUtility(sl.flow_rate));
    end
    new_net_welfare = new_net_welfare - this.getNetworkCost(load);
    % The stop condition
    %    (norm(new_link_price-prices.Link)/this.NumberLinks;
    stop_cond11 = norm(new_link_price-prices.Link)/norm(prices.Link);
    stop_cond12 = norm(new_node_price-prices.Node)/norm(prices.Node);
    %     norm(new_utility-utility)/NS < 10^-3;
    stop_cond2 = norm(new_profit-profit)/norm(profit);
    %     stop_cond3 = norm(new_net_welfare-net_welfare)/NS < 10^-3;
    stop_cond3 = norm(new_net_welfare-net_welfare)/norm(new_net_welfare);
    if ~isempty(DEBUG) && DEBUG
        fprintf('\tObjective value: %d.\n', new_net_welfare);
        fprintf('\t stop condition 1: (%G,%G).\n', stop_cond11,stop_cond12);
        fprintf('\t stop condition 2: %G.\n', stop_cond2);
        fprintf('\t stop condition 3: %G.\n', stop_cond3);
    end
    if (stop_cond11<10^-3&&stop_cond12<10^-3)
        %|| stop_cond2<10^-4 %|| stop_cond3<10^-4
        break;
    else
        prices.Link = alpha*new_link_price+(1-alpha)*prices.Link;
        prices.Node = alpha*new_node_price+(1-alpha)*prices.Node;
        profit = new_profit;
        net_welfare = new_net_welfare;
    end
    number_iter = number_iter + 1;
end
if ~isempty(DEBUG) && DEBUG
    fprintf('\tFirst stage objective value: %d.\n', new_net_welfare);
end

l = 0;
alpha = 1/2;
h = 1;
while true && ~(nnz(a_link_violate)==0 && nnz(a_node_violate)==0)
    prices.Link = init_price.Link;
    prices.Link(a_link_violate) = ...
        prices.Link(a_link_violate) + alpha*delta_link_price(a_link_violate);
    prices.Node = init_price.Node;
    prices.Node(a_node_violate) = ...
        prices.Node(a_node_violate) + alpha*delta_node_price(a_node_violate);
    for s = 1:NS
        this.slices{s}.prices.Link = ...
            prices.Link(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.prices.Node = ...
            prices.Node(this.slices{s}.VirtualNodes.PhysicalNode);
        if nargout == 2
            tic;
        end
        [profit(s), ~, ~] = this.slices{s}.priceOptimalFlowRate([], options);
        if nargout == 2
            t = toc;
            slice_runtime = max(slice_runtime, t);
            runtime.Serial = runtime.Serial + t;
        end
        this.slices{s}.prices.Link = [];
        this.slices{s}.prices.Node = [];
    end
    if nargout == 2
        runtime.Parallel = runtime.Parallel + slice_runtime;
    end
    load = this.getNetworkLoad([], options);
    residual_link_capacity = link_capacity - load.Link;
    residual_node_capacity = node_capacity - load.Node;
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

% Finalize substrate network
this.finalize(prices);

% Calculate the output
output = this.calculateOutput();

% output the optimization results
if ~isempty(DEBUG) && DEBUG
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