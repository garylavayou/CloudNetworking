%% Static Network Slicing
% In the static slicing method, once the resource is allocated to a slice, the allocation
% scheme is not changed during its lifetime.
%
% *NOTE*: when link and node resources are exhausted, some slice request might be
% rejected.
%%
% |options|: if _Method_ is 'slice', then we dimension the slice without pricing,
% otherwise, the fixed pricing policy is adopted (slice-price). [Deprecated]
%
% *TODO*: we can adjust the unit price according to the residual capacity.
function output = staticSlicing(this, slice)

if nargin>=2 && ~isempty(slice)
    %% Allocate Resource to the new arrival slice
    % The residual capacity of the substrate network is available to the slice.
    slice.VirtualLinks.Price = this.getLinkField('Price',slice.VirtualLinks.PhysicalLink);
    slice.VirtualDataCenters.Price = this.getDataCenterField('Price',slice.getDCPI);
    ss = slice.copy;
    ss.VirtualDataCenters.Capacity = this.getDataCenterField('ResidualCapacity', ss.getDCPI);
    ss.VirtualLinks.Capacity = ...
        this.getLinkField('ResidualCapacity', ss.VirtualLinks.PhysicalLink);
    [~] = ss.optimalFlowRate(getstructfields(this.options, ...
        {'Method', 'ConstraintTolerance'}));
    %% Finalize the new slice and the substrate network
    % # After the optimization, the resource allocation variables, flow rate, virtual
    % node/link load of the last slice have been recorded.
    % # Calculate and announce the resource prices to the new slice. The price is fixed in
    % the static slicing method, so the price has been calculated in advance.
    % # Record the substrate network's node/link load, price. When a slice arrive or
    % depart, the network load changes.
    slice.x_path = ss.x_path;
    slice.z_npf = ss.z_npf;
    slice.Variables = ss.Variables;
    slice.setPathBandwidth;
    slice.FlowTable.Rate = ss.FlowTable.Rate;
    slice.VirtualLinks.Load = ss.VirtualLinks.Load;
    slice.VirtualDataCenters.Load = ss.VirtualDataCenters.Load;
    slice.VirtualLinks.Capacity = slice.VirtualLinks.Load;
    slice.VirtualDataCenters.Capacity = slice.VirtualDataCenters.Load;
end

[node_load, link_load] = getNetworkLoad(this);
this.setDataCenterField('Load', node_load);
this.setLinkField('Load', link_load);

% Calculate the output
output = this.calculateOutput([]);
end