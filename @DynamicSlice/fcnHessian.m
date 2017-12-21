%% Hessian matrix of the Lagrangian
% Override <Slice.fcnHessian>, add support for reconfiguration cost in slice dimensioning.
% 
% When resources are reserved, the resource cost is the function of node/link
% capacity, instead of the link/node load. (when no resource reservation,
% load and capacity are the same. In this case, we can add a group of
% equations to make the resource cost also be a function of resource
% capacity, instead of directly function to resource load).    
%
% NOTE: since the hessian matrix will be evaluated many times, to accelerate the
%   computation, we rewrite the method from the superclass, instead of calling the
%   superclass method.
function hess = fcnHessian(vars, lambda, slice, options) %#ok<INUSL>
if isempty(slice.weight)
    weight = slice.FlowTable.Weight;    % for single slice;
else
    weight = slice.weight*ones(slice.NumberFlows, 1);  % for multiple slices
end
NP = slice.NumberPaths;
hess = spalloc(length(vars),length(vars), NP^2);
var_path = vars(1:NP);
for p = 1:NP
    i = slice.path_owner(p);
    hess(p,1:NP) = weight(i)*...
        slice.I_flow_path(i,:)/(1+(slice.I_flow_path(i,:)*var_path))^2; %#ok<SPRIX>
end

if nargin >= 4
    if ~isfield(options, 'PricingPolicy')
        return;
    end
    switch options.PricingPolicy
        case 'quadratic-price'
            var_node = vars((NP+1):slice.num_vars);
            [~,~,lph] = slice.fcnLinkPricing(slice.prices.Link, slice.getLinkLoad(var_path));
            [~,~,nph] = slice.fcnNodePricing(slice.prices.Node, slice.getNodeLoad(var_node));
            NC = slice.NumberDataCenters;
            if isempty(slice.lower_bounds)
                h1 = (slice.I_edge_path') * lph * slice.I_edge_path;
                hess(1:NP, 1:NP) = hess(1:NP, 1:NP) + h1;
                h2 = spalloc(NP*NC, NP*NC, NP*NP*NC);
                z_index1 = 1:NC;
                for p1 = 1:NP
                    z_index2 = 1:NC;
                    for p2 = 1:p1
                        % we only calculate the lower triangle, and use the symetric property
                        % to fill the upper triangle.
                        h2(z_index1,z_index2) = ...
                            diag(slice.I_node_path(:,p1) .* nph .* slice.I_node_path(:,p2)); %#ok<SPRIX>
                        z_index2 = z_index2 + NC;
                    end
                    z_index1 = z_index1 + NC;
                end
                h2 = h2 + (tril(h2,-1))';   % fill the upper triangle since the
                h2 = repmat(h2, slice.NumberVNFs, slice.NumberVNFs);
                hess((NP+1):slice.num_vars, (NP+1):slice.num_vars) = h2;
            else
                % second derviatives of resource cost on (x,z) = 0;
                % second derviatives of resource cost on (c,w(v))
                var_offset = options.num_orig_vars*2;
                c_index = var_offset+(1:slice.NumberVirtualLinks);
                hess(c_index,c_index) = lph;
                v_index = slice.num_vars+(1:slice.num_varv);
                hess(v_index,v_index) = block_diag(diag(nph), slice.NumberVNFs);
            end
        otherwise
    end
end
end
