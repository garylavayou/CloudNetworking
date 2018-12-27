%% DynamicNetwork
% <DynamicNetwork> enables handling of network events.
%
% Subclass should impelement _onAddingSlice_ and _onRemovingSlice_ to handle the network
% dynamics: deciding whether to admit a slice or not, and how to allocate/release
% reosurces.
%% TODO
% (2) resource reallocation: when flows were uncovered by the serving slice (due to new
%     origining or handover),  performing slice dimensioning considering reconfiguration
%     cost.
% (3) decide when to perform resource reallocation, if not unexpected flows or handovers.
%     Options including:
%     (a) event-based: after every N events, performs a reconfiguration;
%     (b) period-based: after every T time, performs a reconfiguration;
%     (c) profit threshold based: after the profit of fast reconfiguration is lower than a
%         threshold, performs a reconfiguration;
%     (d) profit threshold based with prediction: additionaly predict the threshold*.

%% Implementation Issues
% Implementation of this interface should also inherit from <PhysicalNetwork>.
%
classdef (Abstract) DynamicNetwork < PhysicalNetwork & EventSender & EventReceiver
	properties (Constant, Access = protected)
			MIN_NUM_CONFIG = 10;
	end
	properties
		%%
		% for *DeferDimensioning*;
		pending_slices;
	end
	
	events
		FlowArrive;
		FlowDepart;
		AddSliceSucceed;
		AddSliceFailed;
		RemoveSliceSucceed;
		RemoveSliceFailed;          % NOT used.
		AddFlowSucceed;
		AddFlowFailed;
		RemoveFlowSucceed;
		RemoveFlowFailed;          % NOT used.
		ReallocateSlice;
	end
	
	%% Constructor
	methods
		function this = DynamicNetwork(node_opt, link_opt, VNF_opt, net_opt)
			this@PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
			defaultopts = struct(...
				'VNFReconfigCoefficient', 5, ...
				'UnitReconfigureCost', 1 ...
				);
			defaultopts = structupdate(defaultopts, net_opt);
			this.options = setdefault(this.options, defaultopts);
			go = IDynamicSliceOptimizer.GLOBAL_OPTIONS;
			go.eta = this.options.UnitReconfigureCost;
			% assert(isempty(h.eta) || isempty(this.options.UnitReconfigureCost))
			this.pending_slices = ListArray('Slice');
		end
	end
	
	%% Event handler
	methods
		function h = eventhandler(this, source, eventData)
			global DEBUG; %#ok<NUSED>
			% target = eventData.targets;
			% where target should be the <DynamicNetwork> object.
			%% TODO: adding dynamic slice dimensioning for flow arrival/departure scale.
			h = [];         % creatempty('Slice');
			ev = eventData.event;
			switch eventData.EventName
				case 'SliceArrive'
					et = ev.Entity;
					sl = this.AddSlice(et.Options);
					% 2. Allocate flow id: should it be allocated by the Network or
					% allocated by the Event Dispatcher.
					% Inform the <EventDispatcher> to add existing flows and the flow
					% entity builder.
					if isempty(sl)
						notify(this, 'AddSliceFailed');
					else
						sl.Identifier = et.SliceIdentifier;
						if isa(sl, 'DynamicSlice') && sl.isDynamicFlow()
							this.AddListener(sl, {'FlowArrive', 'FlowDepart'}, @sl.eventhandler);
							sl.AddListener(this, {'AddFlowSucceed', 'AddFlowFailed', ...
								'RemoveFlowSucceed', 'RemoveFlowFailed', ...
								'RequestDimensioning', 'DeferDimensioning'}, @this.eventhandler);
						end
						data = DispatchEventData(ev, sl);
						% The |entity| and |slice| information are needed to create flow
						% entitities.
						notify(this, 'AddSliceSucceed', data);
					end
					if nargout >= 1
						h = sl;
					end
				case 'SliceDepart'
					% Remove the slice if 'mandatorydepart' is enable; otherwise, if
					% 'naturaldepart' is enabled, we set a flag. and until the last one
					% flow of this slice departs, we remove the slice.
					sl = ev.userdata;
					if strcmpi(source.flow_depart_option, 'naturaldepart')
						sl.bOnDepart = true;
					else
						this.RemoveSlice(sl);
						data = FlowEventData(ev, sl, []);
						notify(this, 'RemoveSliceSucceed', data);
					end
					if nargout >= 1
						h = sl;
					end
				case 'FlowArrive'
					%%%
					% Notify slice that a flow arrives.
					% |SliceIdentifier| of <SliceEntity> is equal to the |Identifier| of
					% <Slice>.
					% NOTE: handle may also be used in place of |Identifier|.
					% TODO: pass the slice identifier directly through event data.
					slice_id = ev.Entity.Parent.SliceIdentifier;
					sl = this.slices(this.findSlice(slice_id, 'Identifier'));
					[ft, phy_adj] = this.createflow(sl);
					data = FlowEventData(ev, sl, ft, phy_adj);
					%                     data.targets = this.FindSlice(et.Parent.SliceIdentifier);
					notify(this, 'FlowArrive', data);
				case 'FlowDepart'
					% notify slice
					% TODO: remove flow entity child form slice entity
					slice_id = ev.Entity.Parent.SliceIdentifier;
					sl = this.slices(this.findSlice(slice_id, 'Identifier'));
					flow_id = ev.Entity.GlobalIdentifier;
					data = FlowEventData(ev, sl, flow_id);
					notify(this, 'FlowDepart', data);
					if strcmpi(source.flow_depart_option, 'naturaldepart')
						if sl.bOnDepart && sl.NumberFlows == 0
							% Finally remove the slice from the network
							this.RemoveSlice(sl);
							% then remove the slice entity form the entity list.
							data = FlowEventData(ev, sl, []);
							notify(this, 'RemoveSliceSucceed', data);
						end
					end
				case 'AddFlowSucceed'
					% allocate flow id
					% notify <EventDispatcher> to assign identifier to flow entry
					fidx = eventData.flow;
					sl = eventData.slice;
					identifier = this.flow_identifier_generator.next(length(fidx));
					sl.FlowTable{fidx, 'Identifier' } = identifier;
					eventData = DispatchEventData(eventData.event, identifier);
					notify(this, 'AddFlowSucceed', eventData);
				case 'AddFlowFailed'
					% notify <EventDispatcher> to remove the invalid flow entry.
				case 'RemoveFlowSucceed'
					data = EventData(eventData.entity);
					notify(this, 'RemoveFlowSucceed', data);
				case 'RemoveFlowFailed'
				case 'RequestDimensioning'
					sl = source;
					this.pending_slices.Add(sl);
					for i = 1:this.pending_slices.Length
						sl = this.pending_slices{i};
						sl.Optimizer.update_reconfig_costinfo([], true);
					end
					output = this.optimizeResourcePrice(this.pending_slices{:});
					for i = 1:this.pending_slices.Length
						source.Results.Profit = output.Profit(i);
						source.Results.Value = 0;   % TODO: if there are other return values.
					end
					this.pending_slices.Clear();
				case 'DeferDimensioning'
					sl = source;
					this.pending_slices.Add(sl);
					%% TODO
					% decide when to perform dimensioning
					if this.pending_slices.Length >= 3
						output = this.optimizeResourcePrice(this.pending_slices{:});
						for i = 1:this.pending_slices.Length
							sl = this.pending_slices{i};
							sl.Results.Value = 0;
							sl.Results.Profit = output.Profit(i);
						end
						this.pending_slices.Clear;
					end
				otherwise
					error('error: cannot handle event %s.', eventData.EventName);
			end
		end
	end
	
	%% Public methods
	methods
		function sl = AddSlice(this, slice_opt, varargin)
			% We select the method of <DynamicNetwork> to perform adding slice.
			%             AddSlice@CloudNetwork(this, slice_opt, varargin{:});
			sl = AddSlice@PhysicalNetwork(this, slice_opt, varargin{:});
			
			this.onAddingSlice(sl);
		end
		
		function sl = RemoveSlice(this, arg1)
			sl = RemoveSlice@PhysicalNetwork(this, arg1);
			this.optimizeResourcePrice();
			warning('warning: optimization after removing slice.');
		end
		
		function [output, runtime] = optimizeResourcePrice(this, slices, options)
			if nargin <= 1 || isempty(slices)
				slices = this.slices;       % all slices are involved in slice dimensioning
			end
			if nargin <= 2
				options = Dictionary();
			else
				options = Dictionary(options);
			end
			if nargout >= 2
				runtime = 0;
			end
			output = Dictionary();
			
			b_idle_slices = false(length(slices),1);
			for i = 1:length(slices)
				if slices(i).NumberFlows == 0
					b_idle_slices(i) = true;
				end
			end
			idle_slices = slices(b_idle_slices);
			normal_slices = slices(~b_idle_slices);
			num_slices = length(normal_slices);
			if num_slices > 0
				if ~isfield(options, 'InitPrice') || isempty(options.InitPrice)
					% set the initial prices for the slices that need to be re-dimensioned.
					link_prices = zeros(this.NumberLinks, num_slices);
					node_prices = zeros(this.NumberDataCenters, num_slices);
					for i = 1:num_slices
						link_id = normal_slices(i).Links.PhysicalLink;
						link_prices(link_id, i) = normal_slices(i).Links.Price;
						node_id = normal_slices(i).getDCPI;
						node_prices(node_id, i) = normal_slices(i).ServiceNodes.Price;
					end
					options.InitPrice.Link = zeros(this.NumberLinks,1);
					options.InitPrice.Node = zeros(this.NumberDataCenters,1);
					for i = 1:this.NumberLinks
						lp = link_prices(i, link_prices(i,:)~=0);
						if isempty(lp)
							options.InitPrice.Link(i) = 0;
						else
							options.InitPrice.Link(i) = min(lp);
						end
					end
					options.InitPrice.Link = (1/2)*options.InitPrice.Link;
					for i = 1:this.NumberDataCenters
						np = node_prices(i, node_prices(i,:)~=0);
						if isempty(np)
							options.InitPrice.Node(i) = 0;
						else
							options.InitPrice.Node(i) = min(np);
						end
					end
					options.InitPrice.Node = (1/2)*options.InitPrice.Node;
				end
				
				if nargout >= 2
					[output, runtime] = optimizeResourcePrice@PhysicalNetwork...
						(this, normal_slices, options);
				elseif nargout == 1
					output = optimizeResourcePrice@PhysicalNetwork...
						(this, normal_slices, options);
				else
					optimizeResourcePrice@PhysicalNetwork(this, normal_slices, options);
				end
				%% Reconfiguration Cost Model (optional)
				% profit of Slice Customer: utility - resource consumption payment - reconfiguration cost;
				% profit of Slice Provider: resource consumption payment - resource consumption cost =
				%       (resource consumption payment + reconfiguration cost - resource consumption cost
				%       - reconfiguration cost);
				% net social welfare: utility - resource consumption cost - reconfiguration cost.
				%
				% NOTE: the model should be refined to dexcribe the reconfiguration cost.
				%
				% In <CloudNetwork.optimizeResourcePrice> we did not calculate the
				% reconfiguration cost for slices and the net social welfare
				% (see <CloudNetwork.calculateOutput> and <Slice.getProfit>). So we need
				% to append this part of cost.
				% Reconfiguration cost for slices is additionally calculate in
				% <executeMethod>.
				%                 for i = 1:num_slices
				%                     if isa(normal_slices(i), 'DynamicSlice')
				%                         reconfig_cost = normal_slices(i).get_reconfig_cost();
				%                         output.Welfare = output.Welfare - reconfig_cost;
				%                     end
				%                 end
			end
			if ~isempty(idle_slices)
				% recycle all resources
				prices.Link = this.readLink('Price');
				prices.Node = this.readDataCenter('Price');
				for i = 1:length(idle_slices)
					idle_slices(i).finalize(prices);
				end
				% since all resources are released, the profit (reconfiguration cost not
				% included) and cost is zero.
				profit_table = zeros(length(slices), 1);
				if ~isempty(output)
					profit_table(~b_idle_slices) = output.Profit(1:(end-1));
					output.Profit = [profit_table; output.Profit(end)];
				else
					% both slice and network have no profit.
					output = struct('Profit', [profit_table; 0]);
				end
			end
		end
		
		function rc = getReconfigurationCost(this, slices)
			if nargin <= 1
				slices = this.slices;
			end
			
			rc = 0;
			for i = 1:length(slices)
				% Only <DynamicSlice> can perform redimensioning
				% When |b_dim=true|, the slice is going through redimensioning. After
				% redimensioning, |b_dim| is set to 'false'.
				if isa(slices(i), 'DynamicSlice') && slices(i).op.invoke_method >= 2
					if slices(i).isFinal()
						if slices(i).Optimizer.b_dim
							% slice has been finalized, but the |b_dim| flag has not been
							% cleared, we count the true reconfiguration cost.
							rc = rc + slices(i).Optimizer.get_reconfig_cost('const');
						end
						% If the |b_dim| flag has been cleared, we do not count
						% reconfiguration cost (When adding slice).
					else
						% still in optimization stage, we count the approximated
						% reconfiguration cost.
						rc = rc + slices(i).Optimizer.get_reconfig_cost('linear', false);
					end
				end
			end
		end
		
		%% Get network cost
		% Reconfiguration cost is considered. Note: only those slices involved in
		% dimensioning will count the reconfiguration cost.
		% Call this function from superclass will not calculate the reconfiguration cost
		% (unless the 'slices' argument is set improperly, the superclass method has no
		% such an argument.)
		%
		% See also <CloudNetwork.getNetworkCost>.
		function c = totalCost(this, load, reconfig_slices)
			switch nargin
				case 1
					c = totalCost@PhysicalNetwork(this);
				case {2,3}
					c = totalCost@PhysicalNetwork(this, load);
				otherwise
					error('error: unexpected number of input arguments.');
			end
			if nargin >= 3
				c = c + this.getReconfigurationCost(reconfig_slices);
			end
		end
		
		%%
		% see also <DynamicSlice.finalize>.
		% specify the |link_id| and |node_id|, if we need to inquire the removed link and node's
		% reconfiguration cost.
		function updateRedimensionCost(this, slice)
			global DEBUG; %#ok<NUSED>
			link_id = slice.Links.PhysicalLink;
			dc_id = slice.getDCPI();
			link_load = this.readLink('Load', link_id);
			dc_load = this.readDataCenter('Load', dc_id);
			b_zero_load_link = link_load==0;
			b_zero_load_dc = dc_load==0;
			link_load(b_zero_load_link) = this.readLink('Capacity', ...
				link_id(b_zero_load_link)) * (1/20);
			dc_load(b_zero_load_dc) = this.readDataCenter('Capacity', ...
				dc_id(b_zero_load_dc)) * (1/20);
			if ~isempty(slice.Optimizer.prices) && slice.Optimizer.prices.isValid({'Link','Node'})
				link_price = min(slice.Optimizer.prices.Link, this.readLink('Price', link_id));
				node_price = min(slice.Optimizer.prices.Node, this.readDataCenter('Price', dc_id));
			else
				link_price = this.readLink('Price', link_id);
				node_price = this.readDataCenter('Price', dc_id);
			end
			%% ISSUE: HOW TO DETERMINE RECONFIG COST
			[~, slice.Links.ReconfigCost] = ...
				slice.Optimizer.fcnLinkPricing(link_price, link_load);
			num_config = slice.time.DimensionInterval/slice.time.ConfigureInterval;
			num_config = max(DynamicNetwork.MIN_NUM_CONFIG, num_config/4);
			slice.time.DimensionIntervalModified = slice.time.ConfigureInterval*num_config;
			eta = IDynamicSliceOptimizer.GLOBAL_OPTIONS.eta;
			slice.Links.ReconfigCost = (eta/num_config) * slice.Links.ReconfigCost;
			[~, slice.ServiceNodes.ReconfigCost] = ...
				slice.Optimizer.fcnNodePricing(node_price, dc_load);
			slice.ServiceNodes.ReconfigCost = (eta/num_config) * slice.ServiceNodes.ReconfigCost;
		end
		
		
	end
	
	methods (Access = protected)
		%% Perform admitting control and resource allocation
		% Compute the resource allocation for the slice, and decide if this slice can
		% be admitted.
		% NOTE: currently we admit all slice request, so this function also does not
		% do anything.
		function sl = onAddingSlice(this, sl)
			this.pending_slices.Add(sl);
			if ~isa(sl, 'DynamicSlice')
				% We may add static slices to the network. In that case, we will allocate
				% resource resource mannually (e.g, calling _optimizeResourcePriceNew_) or wait
				% until a dynamic slice is added and resource allocation is triggered.
				return;
			end
			
			this.optimizeResourcePrice(this.pending_slices{:});
			this.pending_slices.Clear();
			% At the beginning, the slice is added, without consideration of
			% reconfiguration cost.
		end
		
		function tf = onRemovingSlice(this) %#ok<MANU>
			tf = true;
		end

		%%
		% |finalize| should only be called when dimensiong network slices.
		function finalize(this, prices, sub_slices)
			if nargin <= 2
				sub_slices = this.slices;
			end
			finalize@PhysicalNetwork(this, prices, sub_slices);
		end
		
		function newobj = copyElement(this)
			if this.isShallowCopyable
				newobj = copyElement@PhysicalNetwork(this);
				newobj.isShallowCopyable = false;
				newobj = copyElement@EventSender(newobj);
				newobj.isShallowCopyable = true;
			else
				newobj = this;
			end
			%% Deep Copy Issue.
			% Make a deep copy of the DeepCp object
			% *pending_slices* is just a soft link to data (slices). Therefore, we do not directly
			% call <ListArray.copy> to avoid copying the content in the List. Instead, we will
			% update the corresponding elements with new links to the data.
			newobj.pending_slices = ListArray('DynamicSlice');
			for i = 1:this.pending_slices.Length
				sid = this.FindSlice(this.pending_slices{i});
				newobj.pending_slices.Add(newobj.slices(sid));
			end
			%% Reset the listener of the new instance
			% We should reconfigure the listeners by using AddListeners outside.
			% see <DynamicNetwork>, <EventSender>, <RepeatSliceReconfiguration>.
			newobj.ClearListener();
		end
		
		%%%
		% *Create new flows*
		% Creating new flows in the slice could guarantee no extra node or link would be
		% needed. If we enable new flows from new locations, we should create the flow
		% in the network.
		% |ft|: return flow table entries.
		% overide the default action.
		% (1) simulate flows originating from un-covered stations, which triggers resource
		%     reallocation of the slice.
		%     (a) It should be controlled that only a portion of flows will originated
		%         from uncovered areas, other the reoource reallocation will be too
		%         frequent. The portion can be specified as a paramteter of the slice. We
		%         can monitor the actual portion of the unexpected flows, if the portion
		%         supercedes the specified the flow, and regenerate an expected flow.
		%% TODO
		% Add slice.Options.CandidateNodes;
		function [ft, phy_adjacent] = createflow(this, slice, numflow)
			global DEBUG;
			if nargin <= 2
				numflow = 1;
			end
			assert(isempty(fieldnames(slice.net_changes)), ...
				'error: <slice.net_changes> not reset.');
			%{
				for fi = 1:numflow
					if slice.options.Adhoc == false || ~slice.isNextFlowAdhoc
					else
					end
					%% [TODO] Update the next ad-hoc flow statistics of the slice.
				end
				if ~isempty(find(ft.Type==FlowType.Adhoc,1)) 
					%% [TODO] Update the slice topology
				end
			%}
			if slice.options.Adhoc == false || ~slice.isNextFlowAdhoc
				%% TODO: generate flow from the virtual topology
				A = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
				C = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
				for i = 1:slice.NumberLinks
					h = slice.graph.Head(i);
					t = slice.graph.Tail(i);
					ph = slice.Nodes{h, 'PhysicalNode'};
					pt = slice.Nodes{t, 'PhysicalNode'};
					A(ph, pt) = slice.graph.Adjacent(h,t); %#ok<SPRIX>
					% slice.Topology.Capacity is not available.
					% Some links may have no capacity available, while the descripter is
					% reserved. We assign a minimum value to these zero-capacity links.
					% These links may be allocated real resources if in later stage, we
					% perform dimensioning.
					C(ph, pt) = max(slice.Links{slice.graph.IndexEdge(h,t),'Capacity'},1); %#ok<SPRIX>
				end
				graph = DirectedGraph(A, C);
				slice_opt = getstructfields(slice.options, ...
					{'FlowPattern','DelayConstraint','NumberPaths'});
				slice_opt = this.updateDynamicSliceOptions(slice, slice_opt);
				if nargin <= 2
					slice_opt.NumberFlows = 1;
				else
					slice_opt.NumberFlows = numflow;
				end
				slice_opt.SlicingMethod = SlicingMethod.AdjustPricing;
				slice_opt.DuplicateFlow = true; % if we need to check duplicated flow, do it here.
				b_vailid_flow = false;
				while ~b_vailid_flow
					try
						b_vailid_flow = true;
						ft = this.generateFlowTable(graph, slice_opt);
					catch ME
						if ~isempty(DEBUG) && DEBUG
							disp(ME)
						end
						if strcmp(ME.identifier, 'PhysicalNetwork:Disconnected')
							b_vailid_flow = false;
						else
							rethrow(ME);
						end
					end
				end
				%%%
				% Update slice information as creating slice.
				ft.Properties.VariableNames = ...
					{'Source', 'Target', 'Rate', 'Delay', 'Paths'};
				ft{:,'Type'} = FlowType.Normal;
				if nargout >= 2
					phy_adjacent = [];
				end
			else
				graph = this.graph;
				slice_opt = getstructfields(slice.options, ...
					{'FlowPattern', 'DelayConstraint', 'NumberPaths', 'SlicingMethod'});
				slice_opt = this.updateDynamicSliceOptions(slice, slice_opt);
				slice_opt.NumberFlows = numflow;
				[ft, phy_adjacent] = this.generateFlowTable(graph, slice_opt);
				ft.Properties.VariableNames = ...
					{'Source', 'Target', 'Rate', 'Delay', 'Paths'};
				%%
				% Add tags to the adhoc flows
				ft{:,'Type'} = FlowType.Normal;
				for k = 1:height(ft)
					path_list = ft{k, 'Paths'};
					for p = 1:path_list.Width
						path = path_list{p};
						h = path.node_list(1:(end-1));
						t = path.node_list(2:end);
						phy_eid = this.LinkId(h,t);
						if ~isempty(find(slice.PhysicalLinkMap{phy_eid,'VirtualLink'} == 0,1))
							ft{k,'Type'} = FlowType.Adhoc;
							break;
						end
					end
				end
				if ~isempty(find(ft.Type==FlowType.Adhoc,1))
					error('error: not implemented!');
					%%
					% Back-up: if following flow processing failed, recover slice information.
					% TODO: MOVE to <DynamicSlice>: flowtable + phy_adjacent.
					slice.old_net_state.Nodes = slice.Nodes;
					slice.old_net_state.Links = slice.Links;
					slice.old_net_state.ServiceNodes = slice.ServiceNodes;
					slice.old_net_state.PhysicalNodeMap = slice.PhysicalNodeMap;
					slice.old_net_state.PhysicalLinkMap = slice.PhysicalLinkMap;
					slice.old_net_state.graph = slice.graph.copy;
					%% IMPORTANT INFORMATION ABOUT NEW NODE/LINK ORDER
					% Resource Mapping: In <Slice>(<VirtualNetwork>), we have 'Nodes',
					% 'PhyscialNodeMap', 'Links', 'PhysicalLinkMap', which include the map of
					% virtual resources to physical resource.
					%
					% Here, we add the new nodes/links to the end of the list. As a result,
					% the original virtual node/link indices keep unchanged. i.e.,
					%    [n1,n2,...nN, na, nb, ...]
					%    [e1,e2,...eL, ea, eb, ...]
					% Since new nodes are append to the end of the list, it may be arranged out of
					% physical order (order in the physical network), i.e., node with small physical
					% ID is append to end. The node mapping information can be retrived from the
					% Node table of <VirtualNetwork>.
					%
					% To keep the original link index ([idx, head, tail]) unchanged, we should
					% index links by [head, tail] pairs. As a result, the appended links have larger
					% indices. See also <DirectedGraph>.<Update>.
					%
					% Another solution is to insert the node and links in physical order.
					%  * New node index is determined by physical node index;
					%  * But, we need to record the changes of orignal node/link index, e.g.
					%           old nodes:      1   2   3   4   ... 18        (virtual index)
					%           new nodes:      1   2   +   3   ... +   18    (virtual index)
					%    '+' represent the newly added nodes in physical order. Thus the new
					%    indices of the origin/new nodes becomes
					%           new index:  1   2   4   5   ... 20   | 3	19
					%    and the new node set's old indices is
					%           old index:  1   2   0   3   ... 0   18  (0 means no old index)
					%   the old-new node mapping is similar to the physical-virtual node mapping,
					%   which is used to remap the solution of the last stage (x,z,v). So that we
					%   can compare the old and new solution correctly.
					% On the other hand, the new links are append to the end of the list.
					%
					% The state matrices and the variables (x,z,v) should be augmented, and the
					% orginal ones are embedded in. Then append the incremental information.
					pre_num_nodes = slice.NumberNodes;
					pre_num_edges = slice.NumberLinks;
					pre_num_dcs = slice.NumberDataCenters;
					pre_phy_node_id = slice.Nodes.PhysicalNode;
					pre_phy_head = slice.Nodes{slice.graph.Head, 'PhysicalNode'};
					pre_phy_tail = slice.Nodes{slice.graph.Tail, 'PhysicalNode'};
					b_phy_node = transpose(sum(phy_adjacent,1)~=0) | sum(phy_adjacent,2)~=0;
					b_phy_node(pre_phy_node_id) = 0;
					new_phy_node_id = find(b_phy_node);
					new_num_nodes = numel(new_phy_node_id);
					new_node_index = pre_num_nodes + (1:new_num_nodes)';
					slice.Nodes{new_node_index, :} = 0;
					slice.Nodes{new_node_index, 'PhysicalNode'} = new_phy_node_id;
					% TO BE REMOVED: slice.PhysicalNodeMap{new_phy_node_id, 'VirtualNode'} = new_node_index;
					new_dc_node_index = pre_num_nodes + ...
						find(this.readNode('Capacity', new_phy_node_id) > 0);
					num_new_dcs = length(new_dc_node_index);
					new_dc_index = pre_num_dcs+(1:num_new_dcs)';
					slice.ServiceNodes{new_dc_index, :} = 0;
					slice.ServiceNodes{new_dc_index, 'VirtualNode'} = new_dc_node_index;
					slice.Nodes{new_dc_node_index, 'DataCenter'} = new_dc_index;
					
					%  mask the existing edges in the incident matrix, get new links.
					for i = 1:length(pre_phy_head)
						phy_adjacent(pre_phy_head(i), pre_phy_tail(i)) = 0;
					end
					[new_phy_head, new_phy_tail] = find(phy_adjacent);
					link_map_s2p = this.graph.IndexEdge(new_phy_head,new_phy_tail);
					new_num_edges = length(new_phy_head);
					new_edge_index = pre_num_edges + (1:new_num_edges)';
					slice.Links{new_edge_index, :} = 0;
					slice.Links{new_edge_index, 'PhysicalLink'} = link_map_s2p;
					% TO BE REMOVED: slice.PhysicalLinkMap{link_map_s2p, 'VirtualLink'} = new_edge_index;
					
					% construct new adjacent matrix for the update graph
					new_vhead = slice.PhysicalNodeMap(new_phy_head);
					new_vtail = slice.PhysicalNodeMap(new_phy_tail);
					props.Weight = this.readLink('Weight', ...
						slice.Links{new_edge_index ,'PhysicalLink'});
					%% DISCUSS
					% Adopting the fisrt topology update scheme, the variable changes can be handled
					% more easily.
					slice.graph.Update(new_vhead, new_vtail, props);
					
					%%
					% recorde changes
					% when removing components, the fields with true value correspond to components
					% being removed.
					% Slice will select method, accoding to whether |net_changes| is empty.
					slice.net_changes.NodeIndex = new_node_index;
					slice.net_changes.EdgeIndex = new_edge_index;
					slice.net_changes.DCIndex = new_dc_index;
					
					%                 this.updateRedimensionCost(slice);
					%% Post Processing
					% after the flow is created, the information will be processed by
					% <DynamicSlice>.<OnAddlingFlow>.
					% 
				end
			end
			%% update paths
			% [TODO]: MOVE to <DynamicSlice>: flowtable + phy_adjacent.
			for k = 1:height(ft)
				path_list = ft{k, 'Paths'};
				for p = 1:path_list.Width
					path_list{p}.node_list = slice.PhysicalNodeMap(path_list{p}.node_list);
					path_list{p}.id = this.path_identifier_generator.next;
				end
			end
		end

		% 		function slice_opt = preAddingSlice(this, slice_opt)
		% 			slice_opt = preAddingSlice@PhysicalNetwork(this, slice_opt);
		% 			slice_opt = setdefault(slice_opt);
		% 		end
		
		%%
		% [Experimental]Reconfiguration cost is counted into the slice provider's revenue.
		% Profit of Slice Provider: resource consumption payment - resource consumption cost =
		%       (resource consumption payment + reconfiguration cost - resource consumption cost
		%       - reconfiguration cost);
		function [profit, revenue] = getSliceProviderProfit(this, slices, prices, options)
			if nargin <= 3
				options = Dictionary();
			end
			if nargin <= 1 || isempty(slices)
				slices = this.slices;
			end
			
			[profit, revenue] = getSliceProviderProfit@PhysicalNetwork(this, slices, prices, options);
			
			% In the superclass (PhysicalNetwork) method, the reconfiguration cost is not
			% counted.
			% Therefore, the revenue should be added with reconfiguration cost, while the
			% profit is not changed.
			reconfig_cost = this.getReconfigurationCost(slices);
			revenue = revenue + reconfig_cost;
		end
		
		function argout = calculateOutput(this, argin, options)
			if nargin <= 1 
				argin = Dictionary();
			end
			if nargin <= 2
				options = Dictionary();
			end
			argout = calculateOutput@PhysicalNetwork(this, argin, options);
			
			%%
			% the base method does not count the reconfiguration cost;
			options = getstructfields(options, 'Slices', 'default', {this.slices});
			rc = this.getReconfigurationCost(options.Slices);
			argout.Welfare = argout.Welfare - rc;
			argout.Profit(end) = argout.Profit(end) - rc;
		end
		
		%%%
		% subclass can override this method, Called by _createflow_ method.
		function slice_opt = updateDynamicSliceOptions(this, slice, slice_opt)
			switch slice.options.FlowPattern
				case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
					if slice.options.Adhoc == false || ~slice.isNextFlowAdhoc
						slice_opt.NodeSet = slice.Nodes.PhysicalNode;
					else
						if isfield(slice.options, 'NodeSet')
							% the slice has limited coverage area
							slice_opt.NodeSet = slice.options.NodeSet;  
						else
							% the slice demand can emerge from any node of the network
							slice_opt.NodeSet = 1:slice.Parent.NumberNodes; 
						end
					end
				otherwise
					error('error: cannot handle the flow pattern <%s>.', ...
						slice.options.FlowPattern.char);
			end
			structmerge(slice_opt, getstructfields(slice.options, 'MiddleNodes', ...
						'default-ignore', this.DataCenters.Node));
		end
	end
end

