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
%		results = staticSlicing(this, slices, options)
function varargout = staticSlicing(this, slices, options)
defaultopts = Dictionary('bFixedPrice', false, 'bCountTime', false, 'Fields', 'all',...
	'Invoker', 'StaticScale');
if nargin <= 2
	options = defaultopts;
else
	options = structmerge(defaultopts, options, 'silent');	% merge and override.
end

if nargin >= 2 && ~isempty(slices)
	if options.bFixedPrice
		if options.bCountTime
			t_start = tic;
		end
		assert(~isempty(find([this.readLink('Price');this.readDataCenter('Price')]>0,1)),...
			'error: price is not initialized.');
		Ns = length(slices);
		if Ns == 1
			sl = slices;
			%% Allocate Resource to the new arrival slice
			% The residual capacity of the substrate network is available to the slice.
			sl.Links.Price = this.readLink('Price',sl.Links.PhysicalLink);
			sl.ServiceNodes.Price = this.readDataCenter('Price',sl.getDCPI);
			% ss = slice.copy;
			sl.ServiceNodes.Capacity = this.readDataCenter('ResidualCapacity', sl.getDCPI());
			sl.Links.Capacity = ...
				this.readLink('ResidualCapacity', sl.Links.PhysicalLink);
			sl.Optimizer.setProblem('LinkPrice', sl.Links.Price,...
				'NodePrice', sl.ServiceNodes.Price);  % no need to reset the Problem, no change to network information.
			options = Dictionary(...
				'SlicingMethod', SlicingMethod.AdjustPricing, ...
				'isFinalize', false,...
				'bInitialize', true);
			sl.Optimizer.optimalFlowRate(options);
			%% Finalize the new slice and the substrate network
			% # After the optimization, the resource allocation variables, flow rate, virtual
			% node/link load of the last slice have been recorded.
			% # Calculate and announce the resource prices to the new slice. The price is fixed in
			% the static slicing method, so the price has been calculated in advance.
			% # Record the substrate network's node/link load, price. When a slice arrive or
			% depart, the network load changes.
			sl.finalize();
			sl.Optimizer.setProblem('Price', []);
			load = this.getNetworkLoad();
			this.writeDataCenter('Load', load.Node);
			this.writeLink('Load', load.Link);
			if options.bCountTime
				t_stop = toc(t_start); 
				this.op.runtime = strcut('Parallel', t_stop, 'Serial', t_stop);
			end
			[varargout{1:nargout}] = this.calculateOutput([], options.fields, options); % [1] = this.slices
		else
			[varargout{1:nargout}] = this.singleSliceOptimization(slices, options);
		end
	else
		switch options.Invoker
			case 'StaticScale'
				invoked_method = 'optimizeResourcePriceScaling';
			case 'StaticDual'
				invoked_method = 'optimizeResourcePriceDual';
			case 'StaticResource'
				invoked_method = 'optimizeResourcePrice';
		end
		options = rmfield(options, 'Invoker');
		this.(invoked_method)(slices, options);  % not calculate the output
		[varargout{1:nargout}] = this.calculateOutput(this.slices, options.Fields, ...
			struct('PricingPolicy', 'linear'));
	end
else
	if options.bCountTime
		this.op.runtime = struct('Parallel', 0, 'Serial', 0);
		this.op.iterations = 0;
	end
	[varargout{1:nargout}] = this.calculateOutput(this.slices, options.Fields, ...
		struct('PricingPolicy', 'linear'));
end
end
