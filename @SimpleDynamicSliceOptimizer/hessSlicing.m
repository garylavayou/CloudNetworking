%% Hessian matrix for network slice dimensioning with reconfiguration cost
% This function is used with <fcnProfitReconfigureSlicing>.
% Different from <hessInitialSlicing>, there are auxiliary variables (tx,tz,tv).
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale.
%
% See also <DynamicSlice.fcnProfitReconfigureSlicing>, <DynamicSlice.hessInitialSlicing>.
function hs = hessSlicing(vars, lambda, this, options) %#ok<INUSL>
if isfield(options, 'bCompact') && options.bCompact
    full_vars = zeros(options.num_orig_vars,1);
    full_vars(this.I_active_variables) = vars;
    vars = full_vars;
end

slice = this.hs;
if isempty(this.weight)
    weight = slice.FlowTable.Weight;    % for single slice;
else
    weight = slice.weight*ones(slice.NumberFlows, 1);  % for multiple slices
end
Np = slice.NumberPaths;
Nsn = slice.NumberServiceNodes;
Nvf = slice.NumberVNFs;
Nl = slice.NumberLinks;
hs = spalloc(length(vars),length(vars), Np^2 + Nl + Nvf*Nsn);
var_path = vars(1:Np);
for p = 1:Np
    i = slice.path_owner(p);
    hs(p,1:Np) = weight(i)*...
        this.I_flow_path(i,:)/(1+(this.I_flow_path(i,:)*var_path))^2; %#ok<SPRIX>
end

if nargin >= 4 && isfield(options, 'PricingPolicy')
    switch options.PricingPolicy
        case {'quadratic-price', 'quadratic'}
            var_v_index = (options.num_varx+options.num_varz) + (1:options.num_varv);
            node_load = sum(reshape(vars(var_v_index), Nsn, Nvf),2);
            var_c_index = (options.num_varx+options.num_varz+options.num_varv)*2 + (1:Nl);
            link_load = vars(var_c_index);
            % Since we pricing the resources, we should know amount of resource occupied,
            % instead of the actual load.
            [~,~,lph] = slice.fcnLinkPricing(this.prices.Link, link_load);
            [~,~,nph] = slice.fcnNodePricing(this.prices.Node, node_load);
            % second derviatives of resource cost on (x,z) = 0;
            % second derviatives of resource cost on (c,w(v))
            hs(var_c_index,var_c_index) = diag(lph);
            hs(var_v_index,var_v_index) = block_diag(diag(nph), Nvf);
        otherwise
            error('%s: invalid pricing policy', calledby);
    end
end

if isfield(options, 'bCompact') && options.bCompact
    hs = hs(this.I_active_variables, this.I_active_variables);
end
end
