%% Static Network Slicing
% In the static slicing method, once the resource is allocated to a slice,
% the allocation scheme is not changed during its lifetime.
%
% *NOTE*: when link and node resources are exhausted, some slice request
% might be rejected.
%%
% |options|: 
%		*SlicingMethod*: _StaticPricing_
%
% *TODO*: we can adjust the unit price according to the residual capacity.
function output = staticSlicing(this, slice)
options = getstructfields(this.options, ...
	{'ConstraintTolerance','SlicingMethod','Form'}, 'ignore');
if options.SlicingMethod.IsPricing  % options for _optimalFlowRate_.
    options.PricingPolicy = 'linear';
end

if nargin>=2 && ~isempty(slice)
    %% Allocate Resource to the new arrival slice
    % The residual capacity of the substrate network is available to the slice.
    slice.VirtualLinks.Price = this.getLinkField('Price',slice.VirtualLinks.PhysicalLink);
    slice.VirtualDataCenters.Price = this.getDataCenterField('Price',slice.getDCPI);
    % ss = slice.copy;
    slice.VirtualDataCenters.Capacity = this.getDataCenterField('ResidualCapacity', slice.getDCPI);
    slice.VirtualLinks.Capacity = ...
        this.getLinkField('ResidualCapacity', slice.VirtualLinks.PhysicalLink);
    slice.prices.Link = slice.VirtualLinks.Price;
    slice.prices.Node = slice.VirtualDataCenters.Price;
    [~] = slice.optimalFlowRate(options);
    %% Finalize the new slice and the substrate network
    % # After the optimization, the resource allocation variables, flow rate, virtual
    % node/link load of the last slice have been recorded.
    % # Calculate and announce the resource prices to the new slice. The price is fixed in
    % the static slicing method, so the price has been calculated in advance.
    % # Record the substrate network's node/link load, price. When a slice arrive or
    % depart, the network load changes.
    slice.VirtualLinks.Capacity = slice.VirtualLinks.Load;
    slice.VirtualDataCenters.Capacity = slice.VirtualDataCenters.Load;
    slice.prices.Link = [];
    slice.prices.Node = [];
end

[node_load, link_load] = this.getNetworkLoad;
this.setDataCenterField('Load', node_load);
this.setLinkField('Load', link_load);

% Calculate the output
options.Slices = this.slices;
output = this.calculateOutput([], options);
end
