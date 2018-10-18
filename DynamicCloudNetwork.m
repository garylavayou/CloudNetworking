classdef DynamicCloudNetwork < PhysicalNetwork & DynamicNetwork
    
    properties
        %%
        % for *DeferDimensioning*;
        pending_slices;
    end
    
    events
        ReallocateSlice;
    end
    
    methods
        function h = eventhandler(this, source, eventData)
            switch eventData.EventName
                case 'RequestDimensioning'
                    sl = source;
                    this.pending_slices.Add(sl);
                    output = ...
                        this.optimizeResourcePriceNew([], this.pending_slices{:});
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
                        output = this.optimizeResourcePriceNew([], this.pending_slices{:});
                        for i = 1:this.pending_slices.Length
                            sl = this.pending_slices(i);
                            sl.Results.Value = 0;
                            sl.Results.Profit = output.Profit(i);
                        end
                        this.pending_slices.Clear;
                    end
                otherwise
                    h = eventhandler@DynamicNetwork(this, source, eventData);
            end
            switch eventData.EventName
                case 'SliceArrive'
                    if ~isempty(h)      % h is the added slice
                        sl= h;
                        sl.AddListener(this, {'RequestDimensioning', ...
                            'DeferDimensioning'}, @this.eventhandler);
                    end
                case 'RequestDimensioning'
                    
                otherwise
            end
        end
        
    end
    
    methods
        function this = DynamicCloudNetwork(node_opt, link_opt, VNF_opt, net_opt)
            this@CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
            this@DynamicNetwork([], [], [], net_opt);
            % <DynamicNetwork> has the same field with <PhysicalNetwork> is treated as an
            % interface (without data member). Thus, through <CloudAccessNetwork> the
            % initilization has finished.
            this.pending_slices = ListArray('Slice');
        end
        function sl = AddSlice(this, slice_opt, varargin)
            slice_opt = this.preAddingSlice(slice_opt);
            % We select the method of <DynamicNetwork> to perform adding slice.
            %             AddSlice@CloudNetwork(this, slice_opt, varargin{:});
            sl = AddSlice@DynamicNetwork(this, slice_opt, varargin{:});
            
            this.onAddingSlice(sl);
        end
        
        function sl = RemoveSlice(this, arg1)
            sl = RemoveSlice@CloudNetwork(this, arg1);
            this.optimizeResourcePriceNew();
        end
        
        function [output, runtime] = optimizeResourcePriceNew(this, init_price, slices)
            if nargin <= 2
                slices = this.slices;       % all slices are involved in slice dimensioning
            end
            if nargout == 2
                runtime.Serial = 0;
                runtime.Parallel = 0;
            end
            output = struct([]);
            
            b_idle_slices = false(length(slices),1);
            for i = 1:length(slices)
                if slices{i}.NumberFlows == 0
                    b_idle_slices(i) = true;
                end
            end
            idle_slices = slices(b_idle_slices);
            normal_slices = slices(~b_idle_slices);
            num_slices = length(normal_slices);
            if num_slices > 0
                if isempty(init_price)
                    % set the initial prices for the slices that need to be re-dimensioned.
                    link_prices = zeros(this.NumberLinks, num_slices);
                    node_prices = zeros(this.NumberDataCenters, num_slices);
                    for i = 1:num_slices
                        link_id = normal_slices{i}.VirtualLinks.PhysicalLink;
                        link_prices(link_id, i) = normal_slices{i}.VirtualLinks.Price;
                        node_id = normal_slices{i}.getDCPI;
                        node_prices(node_id, i) = normal_slices{i}.VirtualDataCenters.Price;
                    end
                    init_price.Link = zeros(this.NumberLinks,1);
                    init_price.Node = zeros(this.NumberDataCenters,1);
                    for i = 1:this.NumberLinks
                        lp = link_prices(i, link_prices(i,:)~=0);
                        if isempty(lp)
                            init_price.Link(i) = 0;
                        else
                            init_price.Link(i) = min(lp);
                        end
                    end
                    init_price.Link = (1/2)*init_price.Link;
                    for i = 1:this.NumberDataCenters
                        np = node_prices(i, node_prices(i,:)~=0);
                        if isempty(np)
                            init_price.Node(i) = 0;
                        else
                            init_price.Node(i) = min(np);
                        end
                    end
                    init_price.Node = (1/2)*init_price.Node;
                end
                if nargout >= 2
                    [output, runtime] = optimizeResourcePriceNew@CloudNetwork...
                        (this, init_price, normal_slices);
                elseif nargout == 1
                    output = optimizeResourcePriceNew@CloudNetwork...
                        (this, init_price, normal_slices);
                else
                    optimizeResourcePriceNew@CloudNetwork(this, init_price, normal_slices);
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
                % In <CloudNetwork.optimizeResourcePriceNew> we did not calculate the
                % reconfiguration cost for slices and the net social welfare
                % (see <CloudNetwork.calculateOutput> and <Slice.getProfit>). So we need
                % to append this part of cost.
                % Reconfiguration cost for slices is additionally calculate in
                % <executeMethod>.
                %                 for i = 1:num_slices
                %                     if isa(normal_slices{i}, 'DynamicSlice')
                %                         reconfig_cost = normal_slices{i}.get_reconfig_cost();
                %                         output.Welfare = output.Welfare - reconfig_cost;
                %                     end
                %                 end
            end
            if ~isempty(idle_slices)
                % recycle all resources
                prices.Link = this.readLink('Price');
                prices.Node = this.readDataCenter('Price');
                for i = 1:length(idle_slices)
                    idle_slices{i}.finalize(prices);
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
                if isa(slices{i}, 'DynamicSlice') && slices{i}.invoke_method >= 2
                    if slices{i}.isFinal()
                        if slices{i}.b_dim
                            % slice has been finalized, but the |b_dim| flag has not been
                            % cleared, we count the true reconfiguration cost.
                            rc = rc + slices{i}.get_reconfig_cost('const');
                        end
                        % If the |b_dim| flag has been cleared, we do not count
                        % reconfiguration cost (When adding slice).
                    else
                        % still in optimization stage, we count the approximated
                        % reconfiguration cost.
                        rc = rc + slices{i}.get_reconfig_cost('linear', false);
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
        function c = getNetworkCost(this, node_load, link_load, reconfig_slices)
            switch nargin
                case 1
                    c = getNetworkCost@CloudNetwork(this);
                case 2
                    c = getNetworkCost@CloudNetwork(this, node_load);
                case {3,4,5}
                    c = getNetworkCost@CloudNetwork(this, node_load, link_load);
                otherwise
                    error('%s: unexpected number of input arguments.', calledby);
            end
            if nargin >= 4
                c = c + this.getReconfigurationCost(reconfig_slices);
            end   
        end
        
    end
    
    methods(Access=protected)
        %% Deep Copy
        function newobj = copyElement(this)
            % Make a shallow copy of all properties and perform deep copy for
            % superclasses
            newobj = copyElement@CloudNetwork(this);
						newobj = copyElement@DynamicNetwork(newobj);
						%% Reset the listener of the new instance
						% We should reconfigure the listeners by using AddListeners outside.
						% see <DynamicNetwork>, <EventSender>, <RepeatSliceReconfiguration>.
            %% Deep Copy Issue.
            % Make a deep copy of the DeepCp object
            % *pending_slices* is just a soft link to data (slices). Therefore, we do not directly
            % call <ListArray.copy> to avoid copying the content in the List. Instead, we will
            % update the corresponding elements with new links to the data.
            newobj.pending_slices = ListArray('DynamicSlice');
            for i = 1:this.pending_slices.Length
                sid = this.FindSlice(this.pending_slices{i});
                newobj.pending_slices.Add(newobj.slices{sid});
            end
        end
        
        function tf = onRemovingSlice(this) %#ok<MANU>
            tf = true;
        end
        function sl = createslice(this, slice_opt, varargin)
            if isfield(slice_opt, 'ClassName')
                sl = instantiateclass(slice_opt.ClassName, ...
                    rmfield(slice_opt, 'ClassName'), varargin{:});
            else
                sl = DynamicSlice(slice_opt);
            end
            this.slices{end+1} = sl;
            sl.FlowTable{:, 'Type'} = FlowType.Normal;  % Specify flow type
        end
        
        % Override the inherited function from <CloudNetwork> and <DynamicNetwork>.
        function slice_opt = preAddingSlice(this, slice_opt)
            slice_opt = structmerge(preAddingSlice@CloudNetwork(this, slice_opt),...
                preAddingSlice@DynamicNetwork(this, slice_opt));
        end
        
        function sl = onAddingSlice(this, sl)
            this.pending_slices.Add(sl);
            if ~isa(sl, 'DynamicSlice')
                % We may add static slices to the network. In that case, we will allocate
                % resource resource mannually (e.g, calling _optimizeResourcePriceNew_) or wait
                % until a dynamic slice is added and resource allocation is triggered.
                return;
            end

            this.optimizeResourcePriceNew([], this.pending_slices{:});
            this.pending_slices.Clear();
            % At the beginning, the slice is added, without consideration of
            % reconfiguration cost.
        end
        
        %%
        % |finalize| should only be called when dimensiong network slices.
        function finalize(this, prices, sub_slices)
            if nargin <= 3
                sub_slices = this.slices;
            end
            finalize@CloudNetwork(this, prices, sub_slices);
            %             finalize@DynamicNetwork(this, sub_slices);
        end
        
        % [Experimental]Reconfiguration cost is counted into the slice provider's revenue.
        % Profit of Slice Provider: resource consumption payment - resource consumption cost =
        %       (resource consumption payment + reconfiguration cost - resource consumption cost
        %       - reconfiguration cost);
        function [profit, revenue] = getSliceProviderProfit(this, prices, new_opts)
            if nargin <= 3
                new_opts = struct();
            end
            
            [profit, revenue] = getSliceProviderProfit@CloudNetwork(this, prices, new_opts);
            
            % In the superclass (ClouNetwork) method, the reconfiguration cost is not
            % counted.
            % Therefore, the revenue should be added with reconfiguration cost, while the
            % profit is not changed.
            new_opts = getstructfields(new_opts, {'Slices'}, 'default-ignore', this.slices);
            reconfig_cost = this.getReconfigurationCost(new_opts.Slices);
            revenue = revenue + reconfig_cost;
        end
        
        function argout = calculateOutput(this, argin, options)
            if nargin < 2
                options = struct;
            end
            if nargin < 3
                options = struct;
            end
            argout = calculateOutput@CloudNetwork(this, argin, options);
            
            %%
            % the base method does not count the reconfiguration cost;
            options = getstructfields(options, 'Slices', 'default', this.slices);
            rc = this.getReconfigurationCost(options.Slices);
            argout.Welfare = argout.Welfare - rc;
            argout.Profit(end) = argout.Profit(end) - rc;
        end
    end
    
    methods(Access=protected)
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
            if nargin <= 2
                numflow = 1;
            end
            assert(isempty(fieldnames(slice.net_changes)), ...
                'error: <slice.net_changes> not reset.');
            if slice.options.Adhoc ==false || ~slice.isAdhocFlow
                ft = createflow@DynamicNetwork(this, slice, numflow);
                ft{:,'Type'} = FlowType.Normal;
                if nargout >= 2
                    phy_adjacent = [];
                end
                return;
            end
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
                path_list = ft{k,{'Paths'}};
                for p = 1:path_list.Width
                    path = path_list.paths{p};
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
                %%
                % Back-up: if following flow processing failed, recover slice information.
                % TODO: MOVE to <DynamicSlice>: flowtable + phy_adjacent.
                slice.old_net_state.VirtualNodes = slice.VirtualNodes;
                slice.old_net_state.VirtualLinks = slice.VirtualLinks;
                slice.old_net_state.VirtualDataCenters = slice.VirtualDataCenters;
                slice.old_net_state.PhysicalNodeMap = slice.PhysicalNodeMap;
                slice.old_net_state.PhysicalLinkMap = slice.PhysicalLinkMap;
                slice.old_net_state.Topology = slice.Topology.copy;
                %%
                % new nodes and edges might be added to the slice;
                % In <Slice>(<VirtualNetwork>), we have 'VirtualNodes', 'PhyscialNodeMap',
                % 'VirtualLinks', 'PhysicalLinkMap', which include the map of virtual
                % resources to physical resource.
                %
                % Here, we add the new nodes/links to the end of the list. As a result,
                % the original virtual node/link indices keep unchanged. i.e.,
                %    [n1,n2,...nN, na, nb, ...]
                %    [e1,e2,...eL, ea, eb, ...]
                % Since new nodes/links are append to the end of the list, it may be
                % arranged out of physical order (order in the physical network, e.g. node
                % with small physical ID is append to end).
                % To keep the original link index ([idx, head, tail]) unchanged, we should
                % index the added links separately, instead of column indexing the whole
                % new adjacent matrix. Then append the results to the original link index.
                % *Note* that the resulted link index is no more mapped to the colum
                % indexing of the adjacent matrix.
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
                %    We still perform column-indexing for the link re-mapping procedure,
                %    and identify those new links (not appear in the orignal link set).
                %    Use the identification to generate the new index.
                %  * The new index changes are used to remap the solution of the last
                %    stage (x,z,v). So that we can compare the old and new solution
                %    correctly.
                %
                % The edge-path, node-path matrix might be augmented. Besides, with the
                % second method, we should first perform row exchanges for the matrices,
                % according to the new index. Then append the incremental part.
                %
                % Finally, the components of old solution (x,z,v) should be exchanged
                % according to the new index. Then apply the incremental part.
                pre_num_nodes = slice.NumberVirtualNodes;
                pre_num_edges = slice.NumberVirtualLinks;
                pre_num_dcs = slice.NumberDataCenters;
                pre_phy_node_id = slice.VirtualNodes.PhysicalNode;
                pre_phy_head = slice.VirtualNodes{slice.Topology.Head, 'PhysicalNode'};
                pre_phy_tail = slice.VirtualNodes{slice.Topology.Tail, 'PhysicalNode'};
                b_phy_node = transpose(sum(phy_adjacent,1)~=0) | sum(phy_adjacent,2)~=0;
                b_phy_node(pre_phy_node_id) = 0;
                new_phy_node_id = find(b_phy_node);
                new_num_nodes = numel(new_phy_node_id);
                new_node_index = pre_num_nodes + (1:new_num_nodes)';
                slice.VirtualNodes{new_node_index, :} = 0;
                slice.VirtualNodes{new_node_index, 'PhysicalNode'} = new_phy_node_id;
                slice.PhysicalNodeMap{new_phy_node_id, 'VirtualNode'} = new_node_index;
                new_dc_node_index = pre_num_nodes + ...
                    find(this.readNode('Capacity', new_phy_node_id) > 0);
                num_new_dcs = length(new_dc_node_index);
                new_dc_index = pre_num_dcs+(1:num_new_dcs)';
                slice.VirtualDataCenters{new_dc_index, :} = 0;
                slice.VirtualDataCenters{new_dc_index, 'VirtualNode'} = new_dc_node_index;
                slice.VirtualNodes{new_dc_node_index, 'DataCenter'} = new_dc_index;
                
                %  mask the existing edges in the incident matrix, get new links.
                for i = 1:length(pre_phy_head)
                    phy_adjacent(pre_phy_head(i), pre_phy_tail(i)) = 0;
                end
                [new_phy_head, new_phy_tail] = find(phy_adjacent);
                link_map_s2p = this.graph.IndexEdge(new_phy_head,new_phy_tail);
                new_num_edges = length(new_phy_head);
                new_edge_index = pre_num_edges + (1:new_num_edges)';
                slice.VirtualLinks{new_edge_index, :} = 0;
                slice.VirtualLinks{new_edge_index, 'PhysicalLink'} = link_map_s2p;
                slice.PhysicalLinkMap{link_map_s2p, 'VirtualLink'} = new_edge_index;
                
                % construct new adjacent matrix for the update graph
                new_vhead = slice.PhysicalNodeMap.VirtualNode(new_phy_head);
                new_vtail = slice.PhysicalNodeMap.VirtualNode(new_phy_tail);
                props.Weight = this.readLink('Weight', ...
                    slice.VirtualLinks{new_edge_index ,'PhysicalLink'});
                slice.Topology.Update(new_vhead, new_vtail, props);
                
                %%
                % recorde changes
                % when removing coponents, the fields with true value correspond to components
                % being removed.
                % Slice will select method, accoding to whether |net_changes| is empty.
                slice.net_changes.NodeIndex = false(slice.NumberVirtualNodes,1);
                slice.net_changes.NodeIndex(new_node_index) = true;
                slice.net_changes.EdgeIndex = false(slice.NumberVirtualLinks,1);
                slice.net_changes.EdgeIndex(new_edge_index) = true;
                slice.net_changes.DCIndex = false(slice.NumberDataCenters,1);
                slice.net_changes.DCIndex(new_dc_index) = true;
                
%                 this.updateRedimensionCost(slice);
            end
            %%
            % update paths
            for k = 1:height(ft)
                path_list = ft{k,{'Paths'}};
                for p = 1:path_list.Width
                    path = path_list.paths{p};
                    % Convert to virtual nodes
                    path.node_list = slice.PhysicalNodeMap{path.node_list,'VirtualNode'};
                    path.id = this.path_identifier_generator.next;
                end
            end
        end
        
    end
    
    methods (Access = {?DynamicCloudNetwork,?DynamicSlice})
        %%
        % see also <DynamicSlice.finalize>.
        % specify the |link_id| and |node_id|, if we need to inquire the removed link and node's
        % reconfiguration cost.
        function updateRedimensionCost(this, slice)
            global DEBUG; %#ok<NUSED>
            link_id = slice.VirtualLinks.PhysicalLink;
            node_id = slice.getSNPI;
            link_load = this.readLink('Load', link_id);
            node_load = this.readNode('Load', node_id);
            zero_load_link_id = link_load==0;
            zero_load_node_id = node_load==0;
            link_load(zero_load_link_id) = this.readLink('Capacity', ...
                link_id(zero_load_link_id)) * (1/20);
            node_load(zero_load_node_id) = this.readNode('Capacity', ...
                node_id(zero_load_node_id)) * (1/20);
            if ~isempty(slice.prices)
                link_price = min(slice.prices.Link,this.readLink('Price', link_id));
                node_price = min(slice.prices.Node, this.readNode('Price', node_id));
            else
                link_price = this.readLink('Price', link_id);
                node_price = this.readNode('Price', node_id);
            end
            %% ISSUE: HOW TO DETERMINE RECONFIG COST
            [~, slice.VirtualLinks.ReconfigCost] = ...
                slice.fcnLinkPricing(link_price, link_load);
            num_config = slice.time.DimensionInterval/slice.time.ConfigureInterval;
            MIN_NUM_CONFIG = 10;
            num_config = max(MIN_NUM_CONFIG, num_config/4);
            slice.time.DimensionIntervalModified = slice.time.ConfigureInterval*num_config;
            eta = DynamicSlice.GLOBAL_OPTIONS.get('eta');
            slice.VirtualLinks.ReconfigCost = ...
                (eta/num_config) * slice.VirtualLinks.ReconfigCost;
            [~, slice.VirtualDataCenters.ReconfigCost] = ...
                slice.fcnNodePricing(node_price, node_load);
            slice.VirtualDataCenters.ReconfigCost = ...
                (eta/num_config) * slice.VirtualDataCenters.ReconfigCost;
        end
        %         function [link_reconfig_cost, node_reconfig_cost] = ...
        %                 updateRedimensionCost(this, slice, link_id, node_id)
        %             global DEBUG; %#ok<NUSED>
        %             if nargin <= 2
        %                 link_id = slice.VirtualLinks.PhysicalLink;
        %             end
        %             if nargin <= 3
        %                 node_id = slice.getSNPI;
        %             end
        %             link_load = this.readLink('Load', link_id);
        %             node_load = this.readNode('Load', node_id);
        %             zero_load_link_id = link_load==0;
        %             zero_load_node_id = node_load==0;
        %             link_load(zero_load_link_id) = this.readLink('Capacity', ...
        %                 link_id(zero_load_link_id)) * (1/20);
        %             node_load(zero_load_node_id) = this.readNode('Capacity', ...
        %                 node_id(zero_load_node_id)) * (1/20);
        %             link_price = this.readLink('Price', link_id);
        %             node_price = this.readNode('Price', node_id);
        %             [~, link_reconfig_cost] = slice.fcnLinkPricing(link_price, link_load);
        %             [~, node_reconfig_cost] = slice.fcnNodePricing(node_price, node_load);
        %         end
    end
end

