%% Profit of the slice 
% Evaluate the function value and gradient of the objective.
%
% When prices are specified (explicitly or implicitly), the profit is the total utility
% of the slice minus the total payment to the slice provider; 
% 
% When this method is called as the objective function of optimization, price is provided
% in |options|. On the other hand, when this function is used to evaluate objective value,
% price is provided by the slice |S| by default (although it could also be provided by
% |options|). 
%%
function [profit, grad] = fcnProfit(var_x, S, options)
if isempty(var_x)
    var_path = S.x_path;
    var_node = S.z_npf;
else
    var_path = var_x(1:S.NumberPaths);
    var_node = var_x((S.NumberPaths+1):end);
end
link_load = S.getLinkLoad(var_path);
node_load = S.getNodeLoad(var_node);
flow_rate = S.getFlowRate(var_path);

% determine the prices of link and node.
pricing_policy = '';
if nargin <= 2
    link_price = S.VirtualLinks.Price;
    node_price = S.VirtualDataCenters.Price;
else
    if isfield(options, 'LinkPrice')
        link_price = options.LinkPrice;
    else
        link_price = S.VirtualLinks.Price;
    end
    if isfield(options, 'NodePrice')
        node_price = options.NodePrice;
    else
        node_price = S.VirtualDataCenters.Price;
    end
    if isfield(options, 'PricingPolicy')
        pricing_policy = options.PricingPolicy;
    end
end

switch pricing_policy
    case 'quadratic-price'
        [link_payment,link_price_grad] = S.fcnLinkPricing(link_price, link_load);
        [node_payment,node_price_grad] = S.fcnNodePricing(node_price, node_load);
        profit = -S.weight*sum(fcnUtility(flow_rate)) + link_payment + node_payment;        
    otherwise
        profit = -S.weight*sum(fcnUtility(flow_rate)) ...
            + dot(link_price, link_load) + dot(node_price, node_load);
end
% If there is only one output argument, return the real profit (positive)
if nargout <= 1
    profit = -profit;
else
    %% Gradient value of the objective function
    % *The number of non-zero elements in the gradient vector*:
    %  the gradient on path variable is non-zeros, so there is |P| components;
    %  whether the gradient on node variable is zeros is depend on the node-path
    %  incidence matrix, so the number of non-zero elements is less than i.e.
    %  |nnz(I_node_path)*F|.
    %
    % Gradient of user utility on |x(p)| is given by
    %
    % $$\frac{\partial f}{\partial x(p)} = -\frac{w}{1+\sum_{p_0\in\mathcal{P}_{i_p}}{x_{p_0}}}=
    %   -\frac{w}{1+\sum_{p_0\in\mathcal{P}_{i_p}}{q_{i_p,p_0}\cdot x_{p_0}}} =
    %   -\frac{w}{1+\sum_{p_0\in\mathcal{P}}{q_{i_p,p_0}\cdot x_{p_0}}}$$
    %
    % Gradient of link resource cost on |x(p)| is given by
    % 
    % $$\frac{\partial \rho_L}{\partial x_p}=
    %   \sum_{e\in p}{\frac{\partial \rho_e}{\partial x_p}} =  
    %   \sum_{e\in p}{\left(\frac{\partial \rho_e}{\partial y_e}
    %   \cdot\frac{\partial y_e}{\partial x_p}\right)} = 
    %   \sum_{e\in p}{\rho^{'}_e\cdot g_{e,p}}$$
    %
    % Therefore, the objective function's gradient is given by
    %
    % $$ \frac{\partial f}{\partial x(p)} =
    %    -\frac{w}{1+\sum_{p_0\in\mathcal{P}}{q_{i_p,p_0}\cdot x_{p_0}}} + 
    %    \sum_{e\in p}{\rho^{'}_e\cdot g_{e,p}},~~ \forall p\in\mathcal{P}$$
    grad = spalloc(length(var_x),1, S.NumberPaths+nnz(S.I_node_path)*S.NumberVNFs);
    for p = 1:S.NumberPaths
        i = S.path_owner(p);
        switch pricing_policy
            case 'quadratic-price'
                grad(p) = -S.weight/(1+S.I_flow_path(i,:)*var_path) +  ...
                    dot(link_price_grad,S.I_edge_path(:,p)); %#ok<SPRIX>
            otherwise
               grad(p) = -S.weight/(1+S.I_flow_path(i,:)*var_path) +  ...
                   dot(link_price,S.I_edge_path(:,p)); %#ok<SPRIX>
        end     
    end
    
    %%%
    % and gradient on |z(n,p,f)| is given by
    %
    % $$ \frac{\partial f}{\partial z(n,p,f)} = \rho_n^{'}\cdot h_{n,p} $$
    %
    % To obtain the derivatives on |z|, each step we fix the VNF index |f|. i.e.
    %
    % $$\frac{\partial f}{\partial z(n,p,f)}=\rho_n^{'}\cdot h_{n,p}, 
    %   \forall n \in \mathcal{N}, p \in \mathcal{P}$$
    %
    % So, we can compute the derivatives for all |n| and |p| with fixed |f|, by the price
    % vector and the incident matrix $h_{n,p}$. The results can be directly converted into
    % a column vector. 
    nz = (S.NumberDataCenters*S.NumberPaths);
    z_index = S.NumberPaths+(1:nz);
    for f = 1:S.NumberVNFs
        % compatible arithmetic operation: node_price is a row vector and S.I_node_path is
        % a matrix, and these two operants have the same number of rows.
        switch pricing_policy
            case 'quadratic-price'
                % |grad(z_index)| is a vector, and the right side is a matrix, the value
                % of the matrix will be assigned to |grad(z_index)| column by column.
                grad(z_index) = node_price_grad.*S.I_node_path; %#ok<SPRIX>
            otherwise
                grad(z_index) = node_price.*S.I_node_path; %#ok<SPRIX>
        end
        z_index = z_index + nz;
    end
end
end