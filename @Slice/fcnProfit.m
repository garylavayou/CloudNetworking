%% Evalute the objective function and gradient
%
%%
function [profit, grad]= fcnProfit(var_x, S)
var_path = var_x(1:S.NumberPaths);
var_node = var_x((S.NumberPaths+1):end);
link_price = S.VirtualLinks.Price;
node_price = S.VirtualNodes.Price;

profit = -S.weight*sum(log(S.getFlowRate(var_path))) ...
    + dot(link_price, S.getLinkLoad(var_path)) ...
    + dot(node_price, S.getNodeLoad(var_node));


%% Gradient value of the objective function
% The upper bound number of non-zero elements in the gradient vector:
%  the gradient on path variable is nonzeros, so there is |P| components;
%  whether the gradient on node variable is zeros is depend on the node-path
%  incidence matrix, so the number of non-zero elements is less than i.e.
%  |nnz(I_node_path)*F|.  
%
% Gradient on |x(p)| is given by
%
% $$\frac{\partial f}{\partial x(p)} = -\frac{1}{\sum_{p_0\in\mathcal{P}_{i_p}}{x_{p_0}}}=
%   -\frac{1}{\sum_{p_0\in\mathcal{P}_{i_p}}{q_{i_p,p_0}\cdot x_{p_0}}} =
%   -\sum_{i}{\frac{q_{i,p}}{\sum_{p_0\in\mathcal{P}}{q_{i,p}\cdot x_{p_0}}}}$$
%
% and gradient on |z(n,p,f)| is given by
%
% $$ \frac{\partial f}{\partial z(n,p,f)} = \rho_n\cdot h_{n,p} $$
%
% To obtain the derivatives on |z|, each step we fix the function index |f|. i.e.
%
% $$\frac{\partial f}{\partial z(n,p,f)}=h_{n,p}, \forall n \in \mathcal{N}$$
%
% So, we can compute the derivatives for all |n| and |p| with fixed |f|, which consist a
% matrix $h_{n,p}$, and can be converted into a column vector.
if nargout > 1
    grad = spalloc(length(var_x),1, S.NumberPaths+nnz(S.I_node_path)*S.NumberVNFs);
    for p = 1:S.NumberPaths
        i = S.path_owner(p);
        grad(p) = -S.weight/(S.I_flow_path(i,:)*var_path) + ...
            dot(link_price,S.I_edge_path(:,p)); %#ok<SPRIX>
    end
    
    nz = (S.NumberVirtualNodes*S.NumberPaths);
    z_index = S.NumberPaths+(1:nz);
    for f = 1:S.NumberVNFs
        % compatible arithmetic operation: node_price is a row vector and S.I_node_path is
        % a matrix, and these two operants have the same number of rows.
        grad(z_index) = node_price.*S.I_node_path; %#ok<SPRIX>
        z_index = z_index + nz;
    end
end
end