%% Parameters
traffic_seed = 20160921;
traffic_density = 1;
num_candidate_path = 1;
link_unit_price = 10;
node_unit_cost = 50;
n_type_functions = 6;      % the number of types of virtual network functions
% unit for w(u,i) and w(u)
nf_res_unit = 10;
node_res_unit = 10;
%%
net = FlatNetwork('sample3', traffic_seed, traffic_density);
%%
n_nodes = net.Size;
n_clients = length(net.TrafficProfiles{1,1});
b_client_function = zeros(n_clients, n_type_functions);
for i=1:n_clients
    num_fun_i = randi([2 3],1);
    function_index = unique_randi(n_type_functions, num_fun_i);
    b_client_function(i,function_index) = 1;
end
% cost of transport traffic: d(c,u)
dist_cost = zeros(n_clients, n_nodes);    
% TODO: all-pair path is needed
for i=1:n_clients
    src = net.TrafficProfiles{1,1}(i,1);
    path_list = net.graph.ShortestPathTree(src);
    p = 1;
    for j=1:n_nodes
        if j ~= src
            dist_cost(i, j) = length(path_list{p})*link_unit_price;
            p = p + 1;
        end
    end
end
% unit cost to use a node's resource: p(u,i)
node_cost = repmat(randi([3 6],1, n_type_functions)*node_unit_cost, n_nodes, 1);
node_efficiency = (rand(n_nodes,1)+1);     % 1~2
% network function's resource requirement on each node: w(u,i)
nf_res = repmat(randi([1 5],1, n_type_functions)*nf_res_unit, n_nodes, 1);
for i = 1:net.Size
    nf_res(i,:) = round(nf_res(i,:)*node_efficiency(i));
end
% amount of resource on the nodes: w(u)
node_res = randi([3 10], n_nodes, 1) * node_res_unit;