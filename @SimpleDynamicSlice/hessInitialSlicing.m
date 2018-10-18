%% Hessian matrix of for initial network slice dimensioning
% This function is used with <fcnProfitReserveSlicing>, when slices are initially added
% with explicit resource reservation.
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale.
%
% See also <DynamicSlice.fcnProfitReserveSlicing>, <DynamicSlice.hessSlicing>.
function hs = hessInitialSlicing(vars, lambda, this, options) %#ok<INUSL>
if isfield(options, 'bCompact') && options.bCompact
    full_vars = zeros(options.num_orig_vars,1);
    full_vars(this.I_active_variables) = vars;
    vars = full_vars;
end

if isempty(this.weight)
    weight = this.FlowTable.Weight;    % for single slice;
else
    weight = this.weight*ones(this.NumberFlows, 1);  % for multiple slices
end
NP = this.NumberPaths;
hs = spalloc(length(vars),length(vars), NP^2);
var_path = vars(1:NP);
for p = 1:NP
    i = this.path_owner(p);
    hs(p,1:NP) = weight(i)*...
        this.I_flow_path(i,:)/(1+(this.I_flow_path(i,:)*var_path))^2; %#ok<SPRIX>
end

if nargin >= 4 && isfield(options, 'PricingPolicy')
    switch options.PricingPolicy
        case {'quadratic-price', 'quadratic'}
            ND = this.NumberServiceNodes;
            NV = this.NumberVNFs;
            NL = this.NumberLinks;
            var_v_index = (options.num_varx+options.num_varz) + (1:options.num_varv);
            node_load = sum(reshape(vars(var_v_index), ND, NV),2);
            var_c_index = options.num_varx+options.num_varz+options.num_varv + (1:NL);
            link_load = vars(var_c_index);
            % Since we pricing the resources, we should know amount of resource occupied,
            % instead of the actual load.
            [~,~,lph] = this.fcnLinkPricing(this.prices.Link, link_load);
            [~,~,nph] = this.fcnNodePricing(this.prices.Node, node_load);
            % second derviatives of resource cost on (x,z) = 0;
            % second derviatives of resource cost on (c,w(v))
            hs(var_c_index,var_c_index) = diag(lph);
            hs(var_v_index,var_v_index) = block_diag(diag(nph), NV);
        otherwise
            error('%s: invalid pricing policy', calledby);
    end
end

if isfield(options, 'bCompact') && options.bCompact
    hs = hs(this.I_active_variables, this.I_active_variables);
end
end
