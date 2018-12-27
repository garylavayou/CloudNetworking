%% Hessian matrix for network slice dimensioning with reconfiguration cost
% This function is used with <fcnProfitReconfigureSlicing>.
% Different from <hessInitialSlicing>, there are auxiliary variables (tx,tz,tv).
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale.
%
% See also <DynamicSlice.fcnProfitReconfigureSlicing>, <DynamicSlice.hessInitialSlicing>.
function hs = hessSlicing(vars, lambda, this, options) %#ok<INUSL>
if options.bCompact
    full_vars = sparse(options.num_orig_vars,1);
    full_vars(this.I_active_variables) = vars;
    vars = full_vars;
end

if isempty(this.pardata.Weight)
    weight = this.pardata.FlowWeight;    % for single slice;
else
    weight = this.pardata.Weight*ones(this.pardata.NumberFlows, 1);  % for multiple slices
end
Np = this.pardata.NumberPaths;
Nsn = this.pardata.NumberServiceNodes;
Nvnf = this.pardata.NumberVNFs;
Nl = this.pardata.NumberLinks;
hs = spalloc(length(vars),length(vars), Np^2 + Nl + Nvnf*Nsn);
var_path = vars(1:Np);
for p = 1:Np
    i = this.pardata.PathOwner(p);
    hs(p,1:Np) = weight(i)*...
        this.I_flow_path(i,:)/(1+(this.I_flow_path(i,:)*var_path))^2; %#ok<SPRIX>
end

if nargin >= 4 && isfield(options, 'PricingPolicy')
    switch options.PricingPolicy
        case {'quadratic-price', 'quadratic'}
            var_v_index = sum(this.num_vars(1:2)) + (1:this.num_vars(3));
            node_load = sum(reshape(vars(var_v_index), Nsn, Nvnf),2);
            var_c_index = sum(this.num_vars(1:3))*2 + (1:Nl);
            link_load = vars(var_c_index);
            % Since we pricing the resources, we should know amount of resource occupied,
            % instead of the actual load.
            [~,~,lph] = this.fcnLinkPricing(this.prices.Link, link_load);
            [~,~,nph] = this.fcnNodePricing(this.prices.Node, node_load);
            % second derviatives of resource cost on (x,z) = 0;
            % second derviatives of resource cost on (c,w(v))
            hs(var_c_index,var_c_index) = diag(lph);
            hs(var_v_index,var_v_index) = block_diag(diag(nph), Nvnf);
        otherwise
            error('%s: invalid pricing policy', calledby);
    end
end

if options.bCompact
    hs = hs(this.I_active_variables, this.I_active_variables);
end
end
