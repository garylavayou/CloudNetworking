%% Hessian matrix of the Largrangian
% Called by _singleSliceOptimization_;
%%
function hess = fcnHessian(var_x, ~, S, options)
if isempty(S.weight)
    weight = S.FlowTable.Weight;    % for single slice;
else
    weight = S.weight*ones(S.NumberFlows, 1);  % for multiple slices
end
NP = S.NumberPaths;
hess = spalloc(length(var_x),length(var_x), NP^2);
var_path = var_x(1:NP);
for p = 1:NP
    i = S.path_owner(p);
    f = S.I_flow_path(:,p)~=0;
    hess(p,1:NP) = weight(f)*...
        S.I_flow_path(i,:)/(1+(S.I_flow_path(i,:)*var_path))^2; %#ok<SPRIX>
end
%%%
% Since the problem only contains linear constraint, the hessian matrix of the
% Largrangian is equal to the seconderivatives of the objective function, and the
% Largrangian multipliers $\lambda$ takes no effect.
%
% the Hessian matrix contains only $P^2$ nonzeros elements on the diagonal,
% which is the second derviatives on path variables.

if nargin==4
    if ~isfield(options, 'PricingPolicy')
        return;
    end
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
    switch options.PricingPolicy
        case 'quadratic-price'
            %% quadratic-price
            [~,~,lph] = S.fcnLinkPricing(link_price, S.getLinkLoad(var_path));
            h1 = (S.I_edge_path') * lph * S.I_edge_path;
            hess(1:NP, 1:NP) = hess(1:NP, 1:NP) + h1;
            var_node = var_x((NP+1):end);
            [~,~,nph] = S.fcnNodePricing(node_price, S.getNodeLoad(var_node));
            NC = S.NumberDataCenters;
            h2 = spalloc(NP*NC, NP*NC, NP*NP*NC);
            z_index1 = 1:NC;
            for p1 = 1:NP
                z_index2 = 1:NC;
                for p2 = 1:p1
                    h2(z_index1,z_index2) = ...
                        diag(S.I_node_path(:,p1) .* nph .* S.I_node_path(:,p2));
                    z_index2 = z_index2 + NC;
                end
                z_index1 = z_index1 + NC;
            end
            h2 = h2 + (tril(h2,-1))';
            h2 = repmat(h2, S.NumberVNFs, S.NumberVNFs);
            hess((NP+1):end, (NP+1):end) = hess((NP+1):end, (NP+1):end) + h2;
        otherwise
    end
end
end
