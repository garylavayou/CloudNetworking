%% brute-force search method to find the optimal solution
% **binary search:the state space is too large for a network with 10 nodes, 10 users and 6
% functions. 

opt_objective = inf;
x0 = zeros(n_clients, n_nodes, n_type_functions);
y0 = zeros(n_nodes, n_type_functions);
x_mask = zeros(n_clients, n_nodes, n_type_functions);
y_mask = zeros(n_nodes, n_type_functions);
for c=1:n_clients
    for i=1:n_type_functions
        if b_client_function(c,i)==0    % client c does not use function i.
            x_mask(c,:,i)=1;
        end
    end
end
for u = 1:n_nodes
    for i=1:n_type_functions
        if nf_res(u,i) > node_res(u)    
            % function i's resource requirement on node u is larger than the capacity of
            % node u 
            y_mask(u,i) = 1;
        end
    end
end
b_functions = sum(b_client_function,1)>0;
for i = 1:n_type_functions
    if ~b_functions(i)
        y_mask(:,i) = 1;
    end
end
num_saved_cases = sum(sum(sum(x_mask)))+sum(sum(y_mask));

y = y0;
y(1,:) = b_functions-y_mask(1,:);
b_feasible_solution = false;
while true
    b_node_capacity = true;
    if ~isempty(find(sum(y,1)<b_functions,1))
        b_node_capacity = false;
    end
    if b_node_capacity
        for u = 1:n_nodes
            if dot(y(u,:), nf_res(u,:)) > node_res
                b_node_capacity = false;
                break;
            end
        end
    end
    if b_node_capacity
        x = x0;
        while true
            b_service_dependence = true;
            for c = 1:n_clients
                if ~isempty( find(reshape(x(c,:,:), n_nodes, n_type_functions)>y,1) )
                    b_service_dependence = false;
                    break;
                end
            end
            if b_service_dependence
                b_service_chain = true;
                for c = 1:n_clients
                    for i = 1:n_type_functions
                        if b_client_function(c,i)~=0 && sum(x(c,:,i)) < 1
                            b_service_chain = false;
                            break;
                        end
                    end
                    if ~b_service_chain
                        break;
                    end
                end
                if b_service_chain
                    total_cost = sum(sum(sum(x,3).*dist_cost)) + sum(sum(y.*node_cost));
                    if total_cost < opt_objective
                        opt_objective = total_cost;
                        x_opt = x;
                        y_opt = y;
                    end
                    if ~b_feasible_solution
                        b_feasible_solution = true;
                        x_feasible = x;
                        y_feasible = y;
                        feasible_objective = total_cost;
                    end
                end
            end
            if isempty(find(x_mask+x==0,1))
                break;
            end
            x = IncreaseBinaryVector(x,x_mask);
        end
    end
    if isempty(find(y_mask+y==0,1))
        break;
    end
    y = IncreaseBinaryVector(y,y_mask);
end
