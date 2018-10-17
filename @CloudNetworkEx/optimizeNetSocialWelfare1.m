%% Optimize Net Social Welfare by Dual decomposition
% 
% * *NOTE*: when a point is not actually ascent, do not go backward, and continue searching from
% current point.
%%
function [output] = optimizeNetSocialWelfare1( this )
global DEBUG;
options = getstructfields(this.options, {'PricingFactor', 'PercentFactor'});

% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_capacity = this.readNode('Capacity');
link_capacity = this.readLink('Capacity');
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
        [fval, ~, ~] = ...
            this.slices{s}.subproblemNetSocialWelfare(lambda_s);
        dual_fval = dual_fval + fval;
    end
    dual_fval = dual_fval - dot(lambda.n, node_capacity) - ...
        dot(lambda.e, link_capacity);
    if ~isempty(DEBUG) && DEBUG
        fprintf('\tDual problem: new value: %.3e, old value: %.3e, difference: %.3e.\n', ...
            dual_fval, prev_dual_fval, dual_fval-prev_dual_fval);
    end
    prev_dual_fval = dual_fval;
    % check the primal feasibility
    load = this.getNetworkLoad([], struct('Stage', 'temp'));
    b_feasible = true;
    if ~isempty(find(load.Node>node_capacity,1))
        b_feasible = false;
        lambda.n = lambda.n * 2;
    end
    if ~isempty(find(load.Link>link_capacity,1))
        b_feasible = false;
        lambda.e = lambda.e * 2;
    end
    if b_feasible
        break;
    end
end

% Evaluate the initial step length
delta_lambda.n = load.Node - node_capacity;
delta_lambda.e = load.Link - link_capacity;
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
    delta_lambda.n = load.Node - node_capacity;
    delta_lambda.e = load.Link - link_capacity;
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
        [fval, ~, ~] = ...
            this.slices{s}.subproblemNetSocialWelfare(lambda_s);
        dual_fval = dual_fval + fval;
    end
    dual_fval = dual_fval - dot(lambda.n, node_capacity) - dot(lambda.e, link_capacity);
    load = this.getNetworkLoad([], struct('Stage', 'temp'));
    if ~isempty(DEBUG) && DEBUG
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
link_price = this.getLinkCost + this.lambda.e;
node_price = this.getNodeCost + this.lambda.n;
for s = 1:NS
    sl = this.slices{s};    
    sl.VirtualLinks.Price = link_price(sl.VirtualLinks.PhysicalLink);
    sl.VirtualNodes.Price = node_price(sl.VirtualNodes.PhysicalNode);
end
load = this.getNetworkLoad;
this.writeDataCenter('Load', load.Node);
this.writeLink('Load', load.Link);
this.writeDataCenter('Price', node_price);
this.writeLink('Price', link_price);

%% Calculate the primal optimal value
% If strong duality holds, dual objective value equals to primal objective value.
% We calculate the primal objective value by its definition to verify if the dual can
% primal objective value are consistent.
if this.static_factor ~= 0
    dual_fval = -(prev_dual_fval - this.unitStaticNodeCost*(NN-this.static_factor));
else
    dual_fval = -prev_dual_fval;
end
%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The approximation of net social welfare, and the accurate net social welfare;
% # The net profit of each slice and the substrate network, including four results from
% different methods, i.e. |ApproximatePercent|, |ApproximatePrice|, |AccuratePercent|,
% |AccuratePrice|. 
% # Flow rate of all flows in the network.

%% TODO: replace with calculateOutput
output.node_price = node_price;
output.link_price = link_price;
output.node_load = node_load;
output.link_load = link_load;
output.welfare_approx = 0;
output.welfare_accurate = 0;
t = zeros(this.NumberSlices+1, 1);
output.profit = table(t, t, t, t, 'VariableNames',...
    {'ApproximatePercent', 'ApproximatePrice', 'AccuratePercent', 'AccuratePrice'});
clear t;
options.Model = 'Accurate';
options.bFinal = true;
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
    p = SimpleSlice.fcnNetProfit(var_x, sl);
    output.welfare_approx = output.welfare_approx + p;
    output.profit.ApproximatePercent(s) = options.PercentFactor * p;

    p = SimpleSlice.fcnNetProfit(var_x, sl, options);
    %%% 
    % NOTE it is not accurate to calculate the static cost of each slice with this method.
    if this.static_factor ~= 0
        idx = sl.VirtualNodes.Load>0;
        p = p - dot(sl.VirtualNodes.Load(idx)./this.readNode('Load', nid(idx)),...
            this.readNode('StaticCost', nid(idx)));
    end
    output.profit.AccuratePercent(s) = options.PercentFactor * p;
    %%%
    % the profit of each slice is its utility less the prices pay to substrate network.
    % no matter the cost model is approximate or accurate, the profit is the same.
    % i.e. |profit.ApproximatePrice(s) = profit.ApproximatePrice(s)|
    lambda_s = struct('n', this.lambda.n(nid),'e', this.lambda.e(eid));
    %%%
    % *TODO* DEBUG the difference.
    output.profit.ApproximatePrice(s) = SimpleSlice.subproblemObjective(var_x, lambda_s, sl);
    output.profit.AccuratePrice(s) = sl.getProfit(options);
    
    output.welfare_accurate = output.welfare_accurate...
        + sl.weight*sum(fcnUtility(sl.FlowTable.Rate));
end
if ~isempty(this.eta)
    embed_profit_approx = this.eta*this.getNetworkCost;
    embed_profit_accurate = this.getNetworkCost([], [], options.Model);
else
    embed_profit_approx = 0;
    embed_profit_accurate = 0;   
end
output.welfare_approx = output.welfare_approx + embed_profit_approx;
output.welfare_accurate = output.welfare_accurate ...
    - this.getNetworkCost([],[], options.Model) + embed_profit_accurate;

f = (1-options.PercentFactor)/options.PercentFactor;
% dual_cost = dot(lambda.n, node_load) + dot(lambda.e, link_load);
output.profit.ApproximatePercent(end) = ...
    sum(output.profit.ApproximatePercent(1:(end-1))*f) + embed_profit_approx;
output.profit.AccuratePercent(end) = ...
    sum(output.profit.AccuratePercent(1:(end-1))*f) + embed_profit_accurate;
output.profit.ApproximatePrice(end) = output.welfare_approx - ...
    sum(output.profit.ApproximatePrice(1:(end-1))) + embed_profit_approx;
output.profit.AccuratePrice(end) = output.welfare_accurate - ...
    sum(output.profit.AccuratePrice(1:(end-1))) + embed_profit_accurate;  
if ~isempty(DEBUG) && DEBUG
    fprintf('\tOptimal solution: fx = %G, g(��) = %G.\n', output.welfare_accurate, dual_fval);
    fprintf('\tIteration number: %d, Evaluation Number: %G.\n', iter_num, eval_num);
end
end