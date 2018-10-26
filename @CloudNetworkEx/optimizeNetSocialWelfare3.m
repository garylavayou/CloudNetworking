function [prim_fval] = optimizeNetSocialWelfare3( this )
%% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
%% iteration records
iter_num = 0;
eval_num = 0;
prev_dual_fval = -inf;
step_length = 0.001;
l0 = 1;
temp_lambda.pf = cell(NS,1);
delta_lambda.pf = cell(NS,1);
alpha_f = cell(NS,1);
for s = 1:NS
    NPi = this.slices{s}.NumberPaths;
    NFi = this.slices{s}.NumberVNFs;
    temp_lambda.pf{s} = l0*ones(NPi, NFi);
    alpha_f{s} = this.VNFTable.ProcessEfficiency(this.slices{s}.VNFList);
end
lambda = temp_lambda;
increase_num = 0;
tol_fun = 10;

%% Iteration process
while true
    iter_num = iter_num + 1;
    while true
        eval_num = eval_num + 1;
        [x, dual_val1] = this.path_subproblem(temp_lambda);
        [z, dual_val2] = this.node_subproblem(temp_lambda);
        dual_fval = dual_val1 + dual_val2;
        fprintf('\tDual problem: new value: %.3e, old value: %.3e, difference: %.3e.\n', ...
            dual_fval, prev_dual_fval, dual_fval-prev_dual_fval);
        if dual_fval > prev_dual_fval
            lambda = temp_lambda;
            increase_num = increase_num + 1;
            this.saveStates;
            if increase_num >= 3
                step_length = step_length * (1+2/increase_num);
            end
            break;
        else
            step_length = step_length / 2;
            for s = 1:NS
                temp_lambda.pf{s} = lambda.pf{s} + step_length * delta_lambda.pf{s};
                temp_lambda.pf{s}(lambda.pf{s}<0) = 0;
            end
            increase_num = 0;
        end
    end
    
    if abs(prev_dual_fval-dual_fval) < tol_fun
        break;
    else
        prev_dual_fval = dual_fval;
    end
     
    for s = 1:NS
        NPi = this.slices{s}.NumberPaths;
        NFi = this.slices{s}.NumberVNFs;
        %% TODO: define an access function.
        delta_lambda.pf{s} = reshape(this.slices{s}.As_res * [x{s}; z{s}], NPi, NFi);
        temp_lambda.pf{s} = lambda.pf{s} + step_length * delta_lambda.pf{s};
        temp_lambda.pf{s}(temp_lambda.pf{s}<0) = 0;
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
load = this.getNetworkLoad;
prim_fval = prim_fval - this.totalCost(load);
fprintf('Optimal solution: fx = %G, g(¦Ë) = %G.\n', prim_fval, dual_fval);
fprintf('Iteration number: %d, Evaluation Number: %G.\n', iter_num, eval_num);
end

