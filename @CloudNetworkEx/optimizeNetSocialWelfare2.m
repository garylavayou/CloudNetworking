%% Dual decomposition
% * *step length* when a point is not actually ascent, go backward to the last point, and
% choosing a shorter step lenght and retry. 
function [prim_fval] = optimizeNetSocialWelfare2( this, options )
if nargin <= 1
    options.Display = 'final';
end
%% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
this.clearStates;

%% iteration records
iter_num = 0;
eval_num = 0;
l0 = 10;
lambda.n = l0*ones(NN, 1)*max(node_capacity)./node_capacity;
lambda.e = l0*ones(NL, 1)*max(link_capacity)./link_capacity;
increase_num = 0;
tol_fun = 10^-5;
prev_dual_fval = inf;
%% Find lambda a feasible point
while true
    iter_num = iter_num + 1;
    eval_num = eval_num + 1;
    dual_fval = 0;
    for s = 1:NS
        lambda_s = struct(...
            'n', lambda.n(this.slices{s}.VirtualNodes.PhysicalNode),...
            'e', lambda.e(this.slices{s}.VirtualLinks.PhysicalLink));
        [fval, ~, ~] = this.slices{s}.subproblemNetSocialWelfare(lambda_s);
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
    [node_load, link_load] = this.getNetworkLoad;
    b_feasible = true;
    if ~isempty(find(node_load>node_capacity,1))
        b_feasible = false;
        lambda.n = lambda.n * 2;
    end
    if ~isempty(find(link_load>link_capacity,1))
        b_feasible = false;
        lambda.e = lambda.e * 2;
    end
    if b_feasible
        break;
    end
end
this.saveStates;

%% Evaluate the initial step length    
delta_lambda.n = node_load - node_capacity;
delta_lambda.e = link_load - link_capacity;
idn = delta_lambda.n<0;
ide = delta_lambda.e<0;
step_length = 1/32*min(min(-lambda.e(ide)./delta_lambda.e(ide)),...
    min(-lambda.n(idn)./delta_lambda.n(idn)));
% step_length = 0.001;

%% Iteration process
while true
    iter_num = iter_num + 1;
    temp_lambda.n = lambda.n + step_length*delta_lambda.n;
    temp_lambda.e = lambda.e + step_length*delta_lambda.e;
    temp_lambda.n(temp_lambda.n<0) = 0;
    temp_lambda.e(temp_lambda.e<0) = 0;
    while true
        eval_num = eval_num + 1;
        dual_fval = 0;
        %% solve subproblems
        for s = 1:NS
            lambda_s = struct(...
                'n', temp_lambda.n(this.slices{s}.VirtualNodes.PhysicalNode),...
                'e', temp_lambda.e(this.slices{s}.VirtualLinks.PhysicalLink));
            [fval, ~, ~] = this.slices{s}.subproblemNetSocialWelfare(lambda_s);
            dual_fval = dual_fval + fval;
        end
        dual_fval = dual_fval - dot(temp_lambda.n, node_capacity) - ...
            dot(temp_lambda.e, link_capacity);
        if strncmp(options.Display,'iter', 4)
            fprintf('\tDual problem: new value: %.3e, old value: %.3e, difference: %.3e.\n', ...
                dual_fval, prev_dual_fval, dual_fval-prev_dual_fval);
        end
        if dual_fval > prev_dual_fval
            lambda = temp_lambda;
            increase_num = increase_num + 1;
%             this.saveStates;
            if increase_num >= 3
                step_length = step_length * (1+1/increase_num);
            end
            break;
        elseif abs(prev_dual_fval-dual_fval)/abs(dual_fval) < tol_fun
            break;
        else
            step_length = step_length / 2;
            temp_lambda.n = lambda.n + step_length*delta_lambda.n;
            temp_lambda.e = lambda.e + step_length*delta_lambda.e;
            temp_lambda.n(temp_lambda.n<0) = 0;
            temp_lambda.e(temp_lambda.e<0) = 0;
            increase_num = 0;
        end
    end
    
    if abs(prev_dual_fval-dual_fval)/abs(dual_fval) < tol_fun
        break;
    else
        prev_dual_fval = dual_fval;
    end
    
    [node_load, link_load] = this.getNetworkLoad;
    delta_lambda.n = node_load - node_capacity;
    delta_lambda.e = link_load - link_capacity;
    if isempty(find([delta_lambda.n>0; delta_lambda.e>0],1))
        this.saveStates;
    end
end

%% calculate the primal optimal value
% If strong duality holds, dual objective value equals to primal objective value.
% We calculate the primal objective value by its definition to verify if the dual can
% primal objective value are consistent.
epsilon = this.unitStaticNodeCost;
dual_fval = -(prev_dual_fval - epsilon*(NN-this.static_factor));
prim_fval = 0;
for s = 1:NS
    prim_fval = prim_fval + this.slices{s}.weight*sum(fcnUtility(this.slices{s}.FlowTable.Rate));
end
[node_load, link_load] = this.getNetworkLoad;
this.setNodeField('Load', node_load);
this.setLinkField('Load', link_load);
prim_fval = prim_fval - this.getNetworkCost;
if strncmp(options.Display,'iter', 4) || ...
        strncmp(options.Display,'notify', 6) || strncmp(options.Display,'final', 5) 
    fprintf('\tOptimal solution: fx = %G, g(¦Ë) = %G.\n', prim_fval, dual_fval);
    fprintf('\tIteration number: %d, Evaluation Number: %G.\n', iter_num, eval_num);
end
end

