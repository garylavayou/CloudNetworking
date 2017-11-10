%% Net Social Welfare
% Use dual decomposition to find the optimal solution.
% Compared with the single slice optimization method, this method can decomposit the
% single optimization into several subproblems with smaller scales corresponding to each
% slice. thus this method can facilitate to find the optimal solution. In addition, the
% subproblems can be further decomposited into even smaller scale problems, in which only
% one user flow is involved. For simplicity, we only decompose the problem on the slice
% level.
%
% Despite the advantage of the dual decomposition method, it requires the objective
% function and constraints are separable. When the cost function is linear, this
% requirement is easily meeted. However, when the cost function is strictly convex, the
% resource cost cannot be separated into independent slices.
%% Dual decomposition
% # An initial point of the dual variable is given as ¦Ë=1 (all constraints are active)
% or 0 (the Lagrangian only contains objective function).
% ¦Ë = [¦Ë(p,f,s), ¦Ë(n), ¦Ë(e)]
% # Independently solve the subproblems for each slice.
% # Update the dual variable
% # Repeat the second and third steps, until the optimal solution is reached.
% which can be varified by the function value, or the optimal (KKT) condition.

function [primal_fval] = optimizeNetSocialWelfare( this )
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
epsilon = this.unitStaticNodeCost;
% node_unit_cost = this.getNodeCost;
% phis_n = this.staticNodeCost*this.delta*(NN-1)/this.totalNodeCapacity;
% storage for subproblem results
% the results are mapped from slices to substrate network
%% Find lambda a feasible point
l0 = 10;
lambda.n = l0*ones(NN, 1);
lambda.e = l0*ones(NL, 1);
iter_num = 0;
eval_num = 0;
while true
    iter_num = iter_num + 1;
    eval_num = eval_num + 1;
    dual_fval = 0;
    for i = 1:NS
        lambda_i = struct(...
            'n', lambda.n(this.slices{i}.VirtualNodes.PhysicalNode),...
            'e', lambda.e(this.slices{i}.VirtualLinks.PhysicalLink));
        [fval, ~, ~] = ...
            this.slices{i}.subproblemNetSocialWelfare(lambda_i);
        dual_fval = dual_fval + fval;
        %         total_utility = total_utility + ...
        %             this.slices{i}.weight*sum(fcnUtility(this.slices{i}.flow_rate));
        %         delta_lambda.pf{i} = dg_pf;
        %         step_pf{i} = ...
        %             step_length*ones(1, this.slices{i}.NumberPaths, this.slices{i}.NumberVNFs);
        %         alpha_idx = delta_lambda.pf{i}<0;
        %         step_pf{i}(alpha_idx) = min(step_pf{i}(alpha_idx), ...
        %             -lambda.pf{i}(alpha_idx)./delta_lambda.pf{i}(alpha_idx));
    end
    % check the primal feasibility
    [node_load, link_load] = this.getNetworkLoad([], 'sum');
    b_feasible = true;
    if ~isempty(find(node_load>node_capacity,1))
        b_feasible = false;
        lambda.n = lambda.n*5;
    end
    if ~isempty(find(link_load>link_capacity,1))
        b_feasible = false;
        lambda.e = lambda.e*5;
    end
    if b_feasible
        break;
    end
end
old_dual_fval = dual_fval - dot(lambda.n, node_capacity) - dot(lambda.e, link_capacity);
this.saveStates;
      
step_length = 0.01;
% compute the gradient of lambda.n/e
delta_lambda.n = node_load-node_capacity;
delta_lambda.e = link_load-link_capacity;
partial_step.e = step_length*ones(NL, 1);
partial_step.n = step_length*ones(NN, 1);
th = 100;
temp_lambda.e = lambda.e;
temp_lambda.n = lambda.n;
% states
num_increase = 0;
num_ascent_reject = 0;
b_partial_search = false;
% num_feasible_reject = 0;
while true
    iter_num = iter_num + 1;
    % dual variables update.
    temp_lambda.n = lambda.n + step_length.*delta_lambda.n;
    temp_lambda.n(temp_lambda.n<0) = 0;
    temp_lambda.e = lambda.e + step_length.*delta_lambda.e;
    temp_lambda.e(temp_lambda.e<0) = 0;
    %     step_length = default_step/sqrt(iter_num);     % step for update lambda
    while true
        dual_fval = 0;
        eval_num = eval_num + 1;
        for i = 1:NS
            lambda_i = struct(...
                'n', temp_lambda.n(this.slices{i}.VirtualNodes.PhysicalNode),...
                'e', temp_lambda.e(this.slices{i}.VirtualLinks.PhysicalLink));
            [fval, ~, ~] = ...
                this.slices{i}.subproblemNetSocialWelfare(lambda_i);
            dual_fval = dual_fval + fval;
        end
        % check the ascent direction
        b_ascent = true;
        dual_fval = dual_fval - dot(temp_lambda.n, node_capacity) - ...
            dot(temp_lambda.e, link_capacity);
        % check dual problem's objective value
        diff_value = dual_fval-old_dual_fval;
        fprintf('\tDual problem: new value: %.3e, old value: %.3e, difference: %.3e.\n', ...
            dual_fval, old_dual_fval, diff_value);
        if diff_value < 0
            b_ascent = false;
        end
        
        % check the primal feasibility
        [node_load, link_load] = this.getNetworkLoad([], 'sum');
        ide = link_load>link_capacity;
        idn = node_load>node_capacity;
        b_feasible = true;
        if ~isempty(find(ide,1)) || ~isempty(find(idn,1))
            b_feasible = false;
        end
             
        if b_ascent && b_feasible
            fprintf('\tAccept the step length %G.\n', step_length);
            lambda.e = temp_lambda.e;
            lambda.n = temp_lambda.n;
            old_dual_fval = dual_fval;
            num_increase = num_increase + 1;
            if num_increase >= 3
                step_length = step_length*2;
                num_increase = 0;
                %             else
                %                 step_length = 0.0001/iter_num;
            end
            b_partial_search = false;
            this.saveStates;
            num_ascent_reject = 0;
            break;   
        elseif ~b_ascent && b_feasible
            % While the primal problem is feasible, the dual objective cannot further
            % make any progress. Thus, when the step is smller enough, the current
            % solution can be treated as optimal solution.
            if b_partial_search
                ide = partial_step.e<0;
                idn = partial_step.n<0;
                partial_step.e(ide) = partial_step.e(ide)*2;
                partial_step.e(~ide) = partial_step.e(~ide)/2;
                partial_step.n(idn) = partial_step.n(idn)*2;
                partial_step.n(~idn) = partial_step.n(~idn)/2;
                temp_lambda.e = lambda.e + partial_step.e.*delta_lambda.e;
                temp_lambda.n = lambda.n + partial_step.n.*delta_lambda.n;
            else
                step_length = step_length/2;
                if step_length/norm([lambda.e; lambda.n]) < 10^-8   % stop condition
                    fprintf('\tCannot move forward in the ascent direction ');
                    fprintf('with the given step length tolerance.\n');
                    break;
                end
                temp_lambda.e = lambda.e + step_length.*delta_lambda.e;
                temp_lambda.n = lambda.n + step_length.*delta_lambda.n;
            end
%             b_partial_search = false;
            num_ascent_reject = num_ascent_reject + 1;
            num_increase = 0;
        elseif ~b_ascent && ~b_feasible
            if b_partial_search
                partial_step.e(ide) = partial_step.e(ide)*2;
                partial_step.e(~ide) = partial_step.e(~ide)/2;
                partial_step.n(idn) = partial_step.n(idn)*2;
                partial_step.n(~idn) = partial_step.n(~idn)/2;
                temp_lambda.e = lambda.e + partial_step.e.*delta_lambda.e;
                temp_lambda.n = lambda.n + partial_step.n.*delta_lambda.n;
            else
                step_length = step_length/2;
                if step_length/norm([lambda.e; lambda.n]) < 10^-10
                    fprintf('\tCannot move forward in the ascent direction ');
                    fprintf('with the given step length tolerance.\n');
                    break;
                end
                temp_lambda.e = lambda.e + step_length.*delta_lambda.e;
                temp_lambda.n = lambda.n + step_length.*delta_lambda.n;
            end
            num_ascent_reject = num_ascent_reject + 1;
            num_increase = 0;
        else % b_ascent & ~ b_feasible
            % it is an ascent step, but not a feasible direction.
            % two situation might apply: (1) the step is too large, thus the new point is
            % not a feasible solution. (2) when the step is small enough, the ascent
            % direction cannot be used to update the dual variable.
            %             num_feasible_reject = num_feasbile_reject + 1;
            % those links and nodes which violate the capacity constraints
            if b_partial_search
                %                     if ~isempty(find(ide>0,1))
                partial_step.e(ide) = partial_step.e(ide)*2;
                partial_step.e(~ide) = partial_step.e(~ide)/2;
                %                     end
                %                     if ~isempty(find(idn>0,1))
                partial_step.n(idn) = partial_step.n(idn)*2;
                partial_step.n(~idn) = partial_step.n(~idn)/2;
                %                     end
                %                     if norm([partial_step.e; partial_step.n])/norm([lambda.e; lambda.n])< 10^-6
                %                         break;
                %                     end
                if norm([partial_step.e(~ide); partial_step.n(~idn)])/...
                        norm([lambda.e; lambda.n]) < 10^-10
                    % cannot find feasible direction
                    break;
                end
                temp_lambda.e = lambda.e + partial_step.e.*delta_lambda.e;
                temp_lambda.n = lambda.n + partial_step.n.*delta_lambda.n;
                num_ascent_reject = 0;
            else
                if step_length/norm([lambda.e; lambda.n]) < 10^-6
                    b_partial_search = true;
                    fprintf('\tCannot move forward in the ascent direction with the given step length tolerance. ');
                    fprintf('Start the partial search to find a feasible direction.\n');
                    % move forward in a feasible direction
                    % recover the step.
                    step_length = step_length*2^num_ascent_reject;
                    partial_step.e = ones(NL, 1)*step_length;%*;
                    partial_step.n = ones(NN, 1)*step_length;%*(2^min(num_ascent_reject,3));
                    partial_step.e(ide) = partial_step.e(ide).*...
                        (-link_load(ide)./link_capacity(ide));        
                    partial_step.n(idn) = partial_step.n(idn).*...
                        (-node_load(idn)./node_capacity(idn)-1);
                    temp_lambda.e = lambda.e + partial_step.e.*delta_lambda.e;
                    temp_lambda.n = lambda.n + partial_step.n.*delta_lambda.n;
                    num_ascent_reject = 0;
                else
                    step_length = step_length/2;
                    temp_lambda.e = lambda.e + step_length.*delta_lambda.e;
                    temp_lambda.n = lambda.n + step_length.*delta_lambda.n;
                    num_ascent_reject = num_ascent_reject + 1;
                end
            end
            
            num_increase = 0;
        end
        temp_lambda.e(temp_lambda.e<0) = 0;
        temp_lambda.n(temp_lambda.n<0) = 0;
    end
    
    if abs(diff_value) < th || b_partial_search
        break;
    end
    
    % compute the gradient of lambda.n/e
    delta_lambda.n = node_load-node_capacity;
    delta_lambda.e = link_load-link_capacity;

end

% Output the primal problem's solution
dual_fval = -(old_dual_fval - epsilon*(NN-this.static_factor));
utility = 0;
for i = 1:NS
    utility = utility + this.slices{i}.weight*sum(fcnUtility(this.slices{i}.FlowTable.Rate));
end
[node_load, link_load] = this.getNetworkLoad;
primal_fval = utility - this.getNetworkCost(node_load ,link_load);
fprintf('\tOptimal solution: dual objective: %G, primal objective %G.\n', ...
    dual_fval, primal_fval);
end

%     step_n = step_length*ones(NN,1);
%     alpha_idx = delta_lambda.n<0;
%     step_n(alpha_idx) = ...
%         min(step_n(alpha_idx), -lambda.n(alpha_idx)./delta_lambda.n(alpha_idx));

%% Find the maximum feasible step length in the dual ascent direction.
% # Determine maximum feasible step length (less than the default step length)
% for ¦Ë>=0. the dual variable after the updates must still be in the feasible domain.
% # Determine maximum step length (less than the feasible step length) of ¦Ë_n
% for g(z)>=0.
% # Determine maximum step length of ¦Ë_pf for g(z)>=0.
% After finding the maximum step length for ¦Ë_n, we calculate the step length for ¦Ë_pf
% of each slice, and then update the dual variables using the new step length.
%
% *Note*: finding an conmmon feasible step length for all dual variable is not right,
% since in some direction the dual variable has reached the boundaries of ¦Ë or g(z,¦Ë).
%     for i = 1:NS
%         NPi = this.slices{i}.NumberPaths;
%         NFi = this.slices{i}.NumberVNFs;
%         NNi = this.slices{i}.NumberVirtualNodes;
%         nid = this.slices{i}.VirtualNodes.PhysicalNode;
%         lambda_pf = repmat(reshape(lambda.pf{i}, 1, NPi, NFi), NNi, 1, 1);
%         delta_lambda_pf = reshape(delta_lambda.pf{i}, 1, NPi, NFi);
%         delta_lambda_pf = repmat(delta_lambda_pf, NNi,1,1);
%         % compatible addition , gradz and delta_gradz result in 3D arraies, with dimension of N*P*F
%         gradz = node_unit_cost(nid) + phis_n + lambda.n(nid)-lambda_pf;
%         delta_gradz = delta_lambda.n(nid) - delta_lambda_pf;
%         H_npf = repmat(full(this.slices{i}.I_node_path), 1, 1, NFi);
%
%         alpha_npf = min(step_n(nid).*ones(NNi, NPi, NFi), repmat(step_pf{i}, NNi, 1, 1));
%         % if delta_gradz>=0, then any feasible step lengths still guarantees than gradz is postive.
%         % otherwise, a appropriate value of step length must be calculated.
%         alpha_idx = (delta_gradz<0) & (H_npf~=0);
%         alpha_npf(alpha_idx) = -gradz(alpha_idx)./delta_gradz(alpha_idx);
%         a1 = min(min(alpha_npf, [], 3), [], 2);		% find the maximum (max-min) step length for ¦Ë_n in slice i.
%         % new alpha_npf might be large than the feasible alpha_npf
%         % compare the maximum step length with the prior calculated value and and update it.
%         step_n(nid) = min(step_n(nid), a1);
%     end
%     lambda.n = lambda.n + step_n.*delta_lambda.n;
%
%     for i = 1:NS
%         NPi = this.slices{i}.NumberPaths;
%         NFi = this.slices{i}.NumberVNFs;
%         NNi = this.slices{i}.NumberVirtualNodes;
%         nid = this.slices{i}.VirtualNodes.PhysicalNode;
%         lambda_pf = repmat(reshape(lambda.pf{i}, 1, NPi, NFi), NNi, 1, 1);
%         delta_lambda_pf = reshape(delta_lambda.pf{i}, 1, NPi, NFi);
%         delta_lambda_pf = repmat(delta_lambda_pf, NNi,1,1);
%         % lambda.n(nid) has been updated
%         gradz = node_unit_cost(nid) + phis_n + lambda.n(nid)-lambda_pf;
%         H_npf = repmat(full(this.slices{i}.I_node_path), 1, 1, NFi);
%         % if delta_lambda_pf<= 0, any positive step length for ¦Ë_pf is feasible.
%         % otherwise, we calculate a feasible step length
%         alpha_npf = repmat(step_pf{i}, NNi, 1, 1);
%         alpha_idx = (delta_lambda_pf>0) & (H_npf~=0);
%         alpha_npf(alpha_idx) = gradz(alpha_idx)./delta_lambda_pf(alpha_idx);
%         a2 = min(step_pf{i}, min(alpha_npf, [], 1));
%         a2 = reshape(a2, NPi, NFi);
%         lambda.pf{i} = lambda.pf{i} + a2.*delta_lambda.pf{i};
%     end