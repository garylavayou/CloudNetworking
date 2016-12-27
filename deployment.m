%% Deployment of Network Function
% Determine the association relation between VNF(D-GW) to the NFV enabled nodes, and the
% relationship between D-GW and the end users.
% 
% The problem is formulated as a Mixed Integer Quadratic Constrained Program (see
% 2016-06-30 - Technical Report). To solve this problem, we first relax the problem to a
% version only contain continuous variables. Then this relaxed version is solved by ADMM.

% Note: with appropriate parameter $\rho=0.5$, the iteration of ADMM will converge (with
% iteration number of 80), while the distribute algorithm based on ADMM cannot converge.

%% Network
traffic_seed = [9058,1270,9134,6324,976,2785,5469,9576,9649, 8148];
% servers = [1;4];
servers = [2; 8; 13];
test_network = CloudNetwork('sample2', servers, traffic_seed,[],[]);
L = test_network.LinkNumber;
Ns = test_network.NumberServer;
%% Remove unused nodes
% All nodes are generating traffic.

%% Candidate Path
% *Note*: since the candidate path is seleceted according to haw far is the server node to
% the source node, a path that transits two server nodes does not exist.
test_network.BuildCandidate(4, 1);
P = test_network.path_set.Number;
source_nodes = test_network.source_demands{1}(:,1);
K = length(source_nodes);
d_vec = test_network.source_demands{1}(:,2);
c_vec = zeros(K*Ns+K,1);
index = 1:(Ns+1):((Ns+1)*(K-1)+1);
c_vec(index) = d_vec;

%% Link-path/Flow-path incidence matrix
As = zeros(L, P);
Fs = zeros(K, P);
for k = 1:K
    src = source_nodes(k);
    for sid = 1:test_network.NumberServer
        dest = test_network.server_id(sid);
        path_num = test_network.path_set.Count(src,dest);
        for i = 1:path_num
            path = test_network.path_set.Element{src,dest}(i);
            pid = path.id;
            for l = 1:(path.Length-1)
                e = path.Link(l);
                eid = test_network.graph.Inverse(e(1), e(2));
                As(eid, pid) = 1; 
            end
            Fs(k, pid) = 1;
        end
    end
end
As = sparse(As);
Fs = sparse(Fs);

%% Remove unused links
% idle_link = find(sum(As,2)==0);
% ei = [test_network.graph.Head(idle_link), test_network.graph.Tail(idle_link)];
% test_network.PhysicGraph = Graph(test_network.graph);
% for i = 1:length(ei)
%     test_network.graph.Adjacent(ei(i,1),ei(i,2)) = 0;
% end
% test_network.graph.Renew;

%% Server-link incidence matrix
Hs = zeros(test_network.NumberServer, L);
for sid = 1:test_network.NumberServer
    server = test_network.server_id(sid);
    tail = find(test_network.graph.Adjacent(server,:)~=0);
    for i = 1:length(tail)
        eid = test_network.graph.Inverse(server, tail(i));
        Hs(sid,eid) = 1;
    end
    head = find(test_network.graph.Adjacent(:, server)~=0);
    for i = 1:length(head)
        eid = test_network.graph.Inverse(head(i), server);
        Hs(sid,eid) = 1;
    end
end
Hs = sparse(Hs);

%% Coefficients of ADMM
A_flow = cell(K,1);
A_flow_reverse = cell(K,1);
for k = 1:K
    A_flow{k} = As(:,Fs(k,:)==1);
    A_flow_reverse{k} = [sum(eye(size(A_flow{k},2))/A_flow{k},1); Hs];
    ft = find(abs(A_flow_reverse{k})<10^-5);
    A_flow_reverse{k}(ft) = 0;
    n_path = size(A_flow{k},2);
    r_k = randi([1,10],n_path,1);
    fprintf('Test-Info: path variables ');
    fprintf('%g ', r_k);
    fprintf('\n');
    x_k = A_flow{k}*r_k;
    fprintf('Test-Info: edge variables ');
    fprintf('%g ', x_k);
    fprintf('\n');
    d_k = ones(1,size(A_flow{k},2))*(eye(size(A_flow{k},2))/A_flow{k})*x_k;
    fprintf('Test-Info: demand %g \n', d_k);
end
BigA = block_diag(A_flow_reverse);
BigB = block_diag([zeros(1, Ns); -eye(Ns)], K);

%%
[l1, l2] = size(BigA);
A_hat = spalloc(l1, l2, nnz(BigA));
index_1 = 1:L:((K-1)*L+1);
index_2 = 1:K;
Ahat_column = cell(L,1);
for j = 1:L
    Ahat_column{j} = BigA(:, index_1);
    A_hat(:,index_2) = Ahat_column{j};  %#ok
    index_1 = index_1 + 1;
    index_2 = index_2 + K;
end
[l1, l2] = size(BigB);
B_hat = spalloc(l1, l2, nnz(BigB));
index_1 = 1:Ns:((K-1)*Ns+1);
index_2 = 1:K;
Bhat_column = cell(Ns,1);
for i = 1:Ns
    Bhat_column{i} = BigB(:, index_1);
    B_hat(:,index_2) = Bhat_column{i};  %#ok
    index_1 = index_1 + 1;
    index_2 = index_2 + K;
end

%% Initialization of ADMM
rho = 0.5;
x_iter = cell(L,1);
for j = 1:L
    x_iter{j} = zeros(K,1);
end
z_iter = cell(Ns,1);
for i=1:Ns
    z_iter{i} = zeros(K,1);
end
u_iter = -rho*c_vec;

delta = 0;
for j = 1:L
    delta = delta + Ahat_column{j}*x_iter{j};
end
for i = 1:Ns
    delta = delta + Bhat_column{i}*z_iter{i};
end
delta = delta - c_vec;
stop_cond = norm(delta);

a_link = ones(L,1);
a_server = ones(Ns,1);
A = ones(1,K);
lb = zeros(K,1);

%% Iteration of ADMM
t = 80;
while t
    x_iter_old = x_iter;
    for j = 1:L
        H = rho*(Ahat_column{j}'*Ahat_column{j});
        f = ones(1,K)*a_link(j);
        for j1 = 1:L
            if j1==j
                continue;
            end
            f = f + 0.5*rho*x_iter_old{j1}'*Ahat_column{j1}'*Ahat_column{j};
%             f = f + 0.5*rho*x_iter{j1}'*Ahat_column{j1}'*Ahat_column{j};
        end
        b = 0;
        for i = 1:Ns
            b = b + Bhat_column{i}*z_iter{i};
        end
        b = b - c_vec + u_iter;
        f = f + rho*b'*Ahat_column{j};
        Cj = test_network.graph.Capacity(j);
        if nnz(H)==0
            [x_iter{j}, fval] = linprog(f, A, Cj, [],[], lb,[]);
        else
            [x_iter{j}, fval] = quadprog(H, f, A, Cj, [],[], lb,[]);
        end
    end
    
    z_iter_old = z_iter;
    for i = 1:Ns
        H = rho*(Bhat_column{i}'*Bhat_column{i});
        f = ones(1,K)*a_server(i);
        for i1=1:Ns
            if i==i1
                continue;
            end
            f = f + 0.5*rho*z_iter_old{i1}'*Bhat_column{i1}'*Bhat_column{i};
%             f = f + 0.5*rho*z_iter{i1}'*Bhat_column{i1}'*Bhat_column{i};
        end
        b = 0;
        for j = 1:L
            b = b + Ahat_column{j}*x_iter{j};
        end
        b = b - c_vec + u_iter;
        f = f + rho*b'*Bhat_column{i};
        Vi = test_network.ServerCapacity(i);
        if nnz(H)==0
            [z_iter{i}, fval] = linprog(f, A, Vi, [],[], lb,[]);
        else
            [z_iter{i}, fval] = quadprog(H, f, A, Vi, [],[], lb,[]);
        end
    end
    
    delta = 0;
    for j = 1:L
        delta = delta + Ahat_column{j}*x_iter{j};
    end
    for i = 1:Ns
        delta = delta + Bhat_column{i}*z_iter{i};
    end
    delta = delta - c_vec;
    stop_cond = norm(delta);
    if stop_cond < 10^-3
        break;
    end
    u_iter = u_iter + rho*delta;
    t = t - 1;
end