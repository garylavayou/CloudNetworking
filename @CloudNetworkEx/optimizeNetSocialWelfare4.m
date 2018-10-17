% NOTE: subgradient method
function [prim_fval] = optimizeNetSocialWelfare4( this )
%% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_capacity = this.readNode('Capacity');
link_capacity = this.readLink('Capacity');

%% iteration records
iter_num = 0;
eval_num = 0;
l0 = 10;
lambda.n = l0*ones(NN, 1)*max(node_capacity)./node_capacity;
lambda.e = l0*ones(NL, 1)*max(link_capacity)./link_capacity;
tol_fun = 10;

%% Evaluate the initial step length
dual_fval = 0;
for s = 1:NS
    lambda_s = struct(...
        'n', lambda.n(this.slices{s}.VirtualNodes.PhysicalNode),...
        'e', lambda.e(this.slices{s}.VirtualLinks.PhysicalLink));
    [fval, ~, ~] = this.slices{s}.subproblemNetSocialWelfare(lambda_s);
    dual_fval = dual_fval + fval;
end
prev_dual_fval = dual_fval - dot(lambda.n, node_capacity) - dot(lambda.e, link_capacity);
    
load = this.getNetworkLoad([], struct('Stage', 'temp'));
delta_lambda.n = load.Node-node_capacity;
delta_lambda.e = load.Link-link_capacity;
idn = delta_lambda.n<0;
ide = delta_lambda.e<0;
step_length = min(min(-lambda.e(ide)./delta_lambda.e(ide)),...
    min(-lambda.n(idn)./delta_lambda.n(idn)));
% step_length = 0.001;

%% Iteration process
while true
    iter_num = iter_num + 1;
    eval_num = eval_num + 1;
    dual_fval = 0;
    cur_step = step_length/sqrt(iter_num);
    lambda.n = lambda.n + cur_step*delta_lambda.n;
    lambda.e = lambda.e + cur_step*delta_lambda.e;
    lambda.n(lambda.n<0) = 0;
    lambda.e(lambda.e<0) = 0;
    %% solve subproblems
    for s = 1:NS
        lambda_s = struct(...
            'n', lambda.n(this.slices{s}.VirtualNodes.PhysicalNode),...
            'e', lambda.e(this.slices{s}.VirtualLinks.PhysicalLink));
        [fval, ~, ~] = this.slices{s}.subproblemNetSocialWelfare(lambda_s);
        dual_fval = dual_fval + fval;
    end
    load = this.getNetworkLoad([], struct('Stage', 'temp'));
    dual_fval = dual_fval - dot(lambda.n, node_capacity) - ...
        dot(lambda.e, link_capacity);
    fprintf('\tDual problem: new value: %.3e, old value: %.3e, difference: %.3e.\n', ...
        dual_fval, prev_dual_fval, dual_fval-prev_dual_fval);
    
    if abs(prev_dual_fval-dual_fval) < tol_fun
        break;
    else
        prev_dual_fval = dual_fval;
    end
    delta_lambda.n = load.Node-node_capacity;
    delta_lambda.e = load.Link-link_capacity;
%     if isempty(find([delta_lambda.n>0; delta_lambda.e>0],1))
%         this.saveStates;
%     end
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
load = this.getNetworkLoad;
prim_fval = prim_fval - this.getNetworkCost(load);
fprintf('Optimal solution: fx = %G, g(��) = %G.\n', prim_fval, dual_fval);
fprintf('Iteration number: %d, Evaluation Number: %G.\n', iter_num, eval_num);
end



