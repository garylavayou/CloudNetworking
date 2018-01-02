%% Profit of the slice
% Evaluating the function value and gradient of the objective. 
%
% When prices are specified (explicitly or implicitly), the profit is the total utility of the slice
% minus the total payment to the slice provider; 
% 
%% Input Arguments
% * |vars|: supplied by caller(see also _fmincon_, <Slice.getProfit>); 
% * |slice|: providing slice information, including:
%       *LinkPrice*:
%       *NodePrice*:
% * |options|:
%       *PricingPolicy*:
% NOTE: The second and third arguments must be provided. To accelarate the computation of _fmincon_,
%     we do not perform any arguments checking.
function [profit, grad] = fcnProfit(vars, slice, options)
% determine varriables.
var_path = vars(1:slice.NumberPaths);
var_node = vars((slice.NumberPaths+1):end);
link_load = slice.getLinkLoad(var_path);    % equal to <getLinkCapacity>
node_load = slice.getNodeLoad(var_node);
flow_rate = slice.getFlowRate(var_path);

switch options.PricingPolicy
    case 'quadratic-price'
        [link_payment,link_price_grad] = slice.fcnLinkPricing(slice.prices.Link, link_load);
        [node_payment,node_price_grad] = slice.fcnNodePricing(slice.prices.Node, node_load);
        profit = -slice.weight*sum(fcnUtility(flow_rate)) + link_payment + node_payment;
    otherwise
        profit = -slice.weight*sum(fcnUtility(flow_rate)) ...
            + dot(slice.prices.Link, link_load) + dot(slice.prices.Node, node_load);
end

% When the 'bFinal' option is provided, return the real profit (positive).
if isfield(options, 'bFinal') && options.bFinal
    profit = -profit;
end
if nargout == 2
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
    grad = spalloc(slice.num_vars, 1, ...
        slice.NumberPaths+nnz(slice.I_node_path)*slice.NumberVNFs);
    for p = 1:slice.NumberPaths
        i = slice.path_owner(p);
        switch options.PricingPolicy
            case 'quadratic-price'
                grad(p) = -slice.weight/(1+slice.I_flow_path(i,:)*var_path) +  ...
                    dot(link_price_grad,slice.I_edge_path(:,p)); %#ok<SPRIX>
            otherwise
                grad(p) = -slice.weight/(1+slice.I_flow_path(i,:)*var_path) +  ...
                    dot(slice.prices.Link,slice.I_edge_path(:,p)); %#ok<SPRIX>
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
    nz = (slice.NumberDataCenters*slice.NumberPaths);
    z_index = slice.NumberPaths+(1:nz);
    for f = 1:slice.NumberVNFs
        % compatible arithmetic operation: node_price is a row vector and S.I_node_path is
        % a matrix, and these two operants have the same number of rows.
        switch options.PricingPolicy
            case 'quadratic-price'
                % |grad(z_index)| is a vector, and the right side is a matrix, the value
                % of the matrix will be assigned to |grad(z_index)| column by column.
                grad(z_index) = node_price_grad.*slice.I_node_path; %#ok<SPRIX>
            otherwise
                grad(z_index) = slice.prices.Node.*slice.I_node_path; %#ok<SPRIX>
        end
        z_index = z_index + nz;
    end
end
end