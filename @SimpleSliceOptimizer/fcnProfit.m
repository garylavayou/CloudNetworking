%% Profit of the slice
% Evaluating the function value and gradient of the objective. 
%
% When prices are specified (explicitly or implicitly), the profit is the total utility of
% the slice minus the total payment to the slice provider; 
% 
%% Input Arguments
% * |vars|: supplied by caller(see also _fmincon_, <Slice.getProfit>); 
% * |slice|: providing slice information, including:
%       *LinkPrice*:
%       *NodePrice*:
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
%   # PricingPolicy:
% NOTE: The second and third arguments must be provided. To accelarate the computation of _fmincon_,
%     we do not perform any arguments checking.
function [profit, gd] = fcnProfit(vars, this, options)
if options.bCompact
    full_vars = sparse(options.num_orig_vars,1);
    full_vars(this.I_active_variables) = vars;
    vars = full_vars;
end
Np = this.pardata.NumberPaths;
Nvnf = this.pardata.NumberVNFs;
Nsn = this.pardata.NumberServiceNodes;
weight = this.pardata.Weight;

var_x = vars(1:Np);
var_z = vars((Np+1):(Np+this.num_vars(2)));
link_load = this.getLinkLoad(false, var_x);    % equal to <getLinkCapacity>
node_load = this.getNodeLoad(false, var_z);
flow_rate = this.getFlowRate(var_x);

switch options.PricingPolicy
    case {'quadratic-price', 'quadratic'}
        [link_payment,link_price_grad] = this.fcnLinkPricing(this.prices.Link, link_load);
        [node_payment,node_price_grad] = this.fcnNodePricing(this.prices.Node, node_load);
        profit = -weight*sum(fcnUtility(flow_rate)) + link_payment + node_payment;
    case 'linear'
        profit = -weight*sum(fcnUtility(flow_rate)) ...
            + dot(this.prices.Link, link_load) + dot(this.prices.Node, node_load);
    otherwise
        error('%s: invalid pricing policy', calledby);
end

% When the 'bFinal' option is provided, return the real profit (positive).
if options.bFinal
    profit = -profit;
end
if nargout == 2
    %% Gradient value of the objective function
    % *The number of non-zero elements in the gradient vector*:
    %  the gradient on path variable is non-zeros, so there is |P| components;
    %  whether the gradient on node variable is zeros is depend on the node-path
    %  incidence matrix, so the number of non-zero elements is less than i.e.
    %  |nnz(I_dc_path)*F|.
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
    gd = spalloc(this.num_vars(2), 1, Np+nnz(this.I_dc_path)*Nvnf);
    for p = 1:Np
        i = this.pardata.PathOwner(p);
        switch options.PricingPolicy
            case {'quadratic-price', 'quadratic'}
                gd(p) = -weight/(1+this.I_flow_path(i,:)*var_x) +  ...
                    dot(link_price_grad,this.I_edge_path(:,p)); %#ok<SPRIX>
            case 'linear'
                gd(p) = -weight/(1+this.I_flow_path(i,:)*var_x) +  ...
                    dot(this.prices.Link,this.I_edge_path(:,p)); %#ok<SPRIX>
            otherwise
                error('%s: invalid pricing policy', calledby);
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
    nz = (Nsn*Np);
    z_index = Np+(1:nz);
    for f = 1:Nvnf
        % compatible arithmetic operation: node_price is a row vector and S.I_dc_path is
        % a matrix, and these two operants have the same number of rows.
        switch options.PricingPolicy
            case 'quadratic-price'
                % |grad(z_index)| is a vector, and the right side is a matrix, the value
                % of the matrix will be assigned to |grad(z_index)| column by column.
                gd(z_index) = node_price_grad.*this.I_dc_path; %#ok<SPRIX>
            case 'linear'
                gd(z_index) = this.prices.Node.*this.I_dc_path; %#ok<SPRIX>
            otherwise
                error('%s: invalid pricing policy', calledby);
        end
        z_index = z_index + nz;
    end
    
    if options.bCompact
        gd = gd(this.I_active_variables);
    end
end
end