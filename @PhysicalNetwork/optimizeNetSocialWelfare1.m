%% Optimize Net Social Welfare by Dual decomposition
% 
% * *NOTE*: when a point is not actually ascent, do not go backward, and continue searching from
% current point.
%%
function [output] = optimizeNetSocialWelfare1( this, options )
if nargin <= 1
    options.Display = 'final';
end

% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_load = zeros(NN, NS);
link_load = zeros(NL, NS);
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
this.clearStates;

% iteration records
iter_num = 0;
eval_num = 0;
l0 = 10;
lambda.n = l0*ones(NN, 1)*max(node_capacity)./node_capacity;
lambda.e = l0*ones(NL, 1)*max(link_capacity)./link_capacity;
increase_num = 0;
tol_fun = 10^-4;
prev_dual_fval = inf;

% Find lambda a feasible point
while true
    iter_num = iter_num + 1;
    eval_num = eval_num + 1;
    dual_fval = 0;
    for s = 1:NS
        lambda_s = struct(...
            'n', lambda.n(this.slices{s}.VirtualNodes.PhysicalNode),...
            'e', lambda.e(this.slices{s}.VirtualLinks.PhysicalLink));
        [fval, node_load(:,s), link_load(:,s)] = ...
            this.slices{s}.subproblemNetSocialWelfare(lambda_s);
        dual_fval = dual_fval + fval;
    end
    dual_fval = dual_fval - dot(lambda.n, node_capacity) - ...
        dot(lambda.e, link_capacity);
    if strncmp(options.Display,'iter', 4)
        fprintf('\tDual problem: new value: %.3e, old value: %.3e, difference: %.3e.\n', ...
            dual_fval, prev_dual_fval, dual_fval-prev_dual_fval);
    end
    prev_dual_fval = dual_fval;
    % check the primal feasibility
    aggr_node_load = sum(node_load, 2);
    b_feasible = true;
    if ~isempty(find(aggr_node_load>node_capacity,1))
        b_feasible = false;
        lambda.n = lambda.n * 2;
    end
    aggr_link_load = sum(link_load, 2);
    if ~isempty(find(aggr_link_load>link_capacity,1))
        b_feasible = false;
        lambda.e = lambda.e * 2;
    end
    if b_feasible
        break;
    end
end

% Evaluate the initial step length
delta_lambda.n = aggr_node_load - node_capacity;
delta_lambda.e = aggr_link_load - link_capacity;
idn = delta_lambda.n<0;
ide = delta_lambda.e<0;
step_length = 1/8*min(min(-lambda.e(ide)./delta_lambda.e(ide)),...
    min(-lambda.n(idn)./delta_lambda.n(idn)));
% step_length = 0.001;

% Iteration process
while true
    iter_num = iter_num + 1;
    eval_num = eval_num + 1;
    dual_fval = 0;
    delta_lambda.n = sum(node_load, 2) - node_capacity;
    delta_lambda.e = sum(link_load, 2) - link_capacity;
    if isempty(find([delta_lambda.n>0; delta_lambda.e>0],1))
        this.saveStates(lambda);
    end
    lambda.n = lambda.n + step_length*delta_lambda.n;
    lambda.e = lambda.e + step_length*delta_lambda.e;
    lambda.n(lambda.n<0) = 0;
    lambda.e(lambda.e<0) = 0;
    % solve subproblems
    for s = 1:NS
        lambda_s = struct(...
            'n', lambda.n(this.slices{s}.VirtualNodes.PhysicalNode),...
            'e', lambda.e(this.slices{s}.VirtualLinks.PhysicalLink));
        [fval, node_load(:,s), link_load(:,s)] = ...
            this.slices{s}.subproblemNetSocialWelfare(lambda_s);
        dual_fval = dual_fval + fval;
    end
    dual_fval = dual_fval - dot(lambda.n, node_capacity) - dot(lambda.e, link_capacity);
    if strncmp(options.Display,'iter', 4)
        fprintf('\tDual problem: new value: %.3e, old value: %.3e, difference: %.3e.\n', ...
            dual_fval, prev_dual_fval, dual_fval-prev_dual_fval);
    end
    if dual_fval > prev_dual_fval
        increase_num = increase_num + 1;
        if increase_num >= 3
            step_length = step_length * (1+1/increase_num);
        end
    elseif abs(prev_dual_fval-dual_fval)/abs(dual_fval) < tol_fun
        break;
    else
        % search for better step diminish strategy.
        step_length = step_length / 2;
        increase_num = 0;
    end    
    prev_dual_fval = dual_fval;
end

%% Finalize substrate network
% # After the optimization, the resource allocation variables, flow rate, node/link load
% have been recorded (by the saveStates method). 
% # Calculate and announce the resource prices to each slice.
% # Record the substrate network's node/link load, price.
epsilon = this.unitStaticNodeCost;
phis_n = epsilon*this.delta*(NN-1)/this.totalNodeCapacity;
phis_l = epsilon*(1-this.delta)*(NN-1)/this.totalLinkCapacity;
link_price = this.getLinkField('UnitCost') + phis_l + this.lambda.e;
node_price = this.getNodeField('UnitCost') + phis_n + this.lambda.n;
for s = 1:NS
    sl = this.slices{s};    
    sl.VirtualLinks.Price = link_price(sl.VirtualLinks.PhysicalLink);
    sl.VirtualNodes.Price = node_price(sl.VirtualNodes.PhysicalNode);
    node_load(sl.VirtualNodes.PhysicalNode, s) = sl.getNodeLoad;
    link_load(sl.VirtualLinks.PhysicalLink, s) = sl.getLinkLoad;
end
aggr_node_load = sum(node_load,2);
aggr_link_load = sum(link_load,2);
this.setNodeField('Load', aggr_node_load);
this.setLinkField('Load', aggr_link_load);
this.setNodeField('Price', node_price);
this.setLinkField('Price', link_price);

%% Calculate the primal optimal value
% If strong duality holds, dual objective value equals to primal objective value.
% We calculate the primal objective value by its definition to verify if the dual can
% primal objective value are consistent.
dual_fval = -(prev_dual_fval - epsilon*(NN-1));
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
    nid = sl.VirtualNodes.PhysicalNode;
    eid = sl.VirtualLinks.PhysicalLink; 
    output.flow_rate = [output.flow_rate; sl.FlowTable.Rate];
    %%%
    % _subproblemObjective_ calculates the net profit of each slice with the given dual
    % variables. Therefore, it is not necessary to manually calculate it with price and
    % load. 
    var_x = [sl.Variables.x; sl.Variables.z];
    %%%
    % _fcnNetWelfare_ calculates each slices net welfare, the sum of which is equal to the
    % network's total net socal welfare.
    p = -Slice.fcnNetProfit(var_x, sl);
    output.welfare_approx = output.welfare_approx + p;
    output.profit.ApproximatePercent(s) = options.PercentFactor * p;

    p = -Slice.fcnNetProfit(var_x, sl, options);
    idx = sl.VirtualNodes.Load>0;
    %%% 
    % NOTE it is not accurate to calculate the static cost of each slice with this method.
    p = p - dot(sl.VirtualNodes.Load(idx)./this.getNodeField('Load', nid(idx)),...
        this.getNodeField('StaticCost', nid(idx)));
    output.profit.AccuratePercent(s) = options.PercentFactor * p;
    %%%
    % the profit of each slice is its utility less the prices pay to substrate network.
    % no matter the cost model is approximate or accurate, the profit is the same.
    % i.e. |profit.ApproximatePrice(s) = profit.ApproximatePrice(s)|
    lambda_s = struct('n', this.lambda.n(nid),'e', this.lambda.e(eid));
    %%%
    % *TODO* DEBUG the difference.
    output.profit.ApproximatePrice(s) = - Slice.subproblemObjective(var_x, lambda_s, sl);
    output.profit.AccuratePrice(s) = - Slice.fcnProfit(var_x, sl);
    
    output.welfare_accurate = output.welfare_accurate...
        + sl.weight*sum(log(sl.FlowTable.Rate));
end
f = (1-options.PercentFactor)/options.PercentFactor;
% dual_cost = dot(lambda.n, node_load) + dot(lambda.e, link_load);
output.profit.ApproximatePercent(end) = ...
    sum(output.profit.ApproximatePercent(1:(end-1))*f);
output.profit.AccuratePercent(end) = ...
    sum(output.profit.AccuratePercent(1:(end-1))*f);
output.profit.ApproximatePrice(end) = output.welfare_approx - ...
    sum(output.profit.ApproximatePrice(1:(end-1)));
output.welfare_accurate = output.welfare_accurate ...
    - this.getNetworkCost([],[], options.Model);
output.profit.AccuratePrice(end) = output.welfare_accurate - ...
    sum(output.profit.AccuratePrice(1:(end-1)));  
if ~strcmp(options.Display,'off') && ~strcmp(options.Display,'none')
    fprintf('\tOptimal solution: fx = %G, g(¦Ë) = %G.\n', output.welfare_accurate, dual_fval);
    fprintf('\tIteration number: %d, Evaluation Number: %G.\n', iter_num, eval_num);
end
end