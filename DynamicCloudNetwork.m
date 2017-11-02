classdef DynamicCloudNetwork < CloudNetwork & DynamicNetwork
    
    properties
        %%
        % for *DeferDimensioning*;
        pending_slices = ListArray('DynamicSlice');    
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
                    source.Results.Profit = output.Profit(1);
                    source.Results.Value = 0;   % if there are other return values.
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
        end

        function sl = AddSlice(this, slice_opt, varargin)
            slice_opt = this.preAddingSlice(slice_opt);
            % We select the method of <DynamicNetwork> to perform adding slice.
            %             AddSlice@CloudNetwork(this, slice_opt, varargin{:});
            sl = AddSlice@DynamicNetwork(this, slice_opt, varargin{:});
        end   
        
        function sl = RemoveSlice(this, arg1)
            sl = RemoveSlice@CloudNetwork(this, arg1);
            this.optimizeResourcePriceNew();
        end
        
    end
    
    methods(Access=protected)
        function tf = onRemovingSlice(this) %#ok<MANU>
            tf = true;
        end
        function sl = createslice(this, slice_opt)
            this.slices{end+1} = DynamicSlice(slice_opt);
            sl = this.slices{end};
            sl.FlowTable{:, 'Type'} = FlowType.Normal;  % Specify flow type
        end
        
        % Override the inherited function.
        %         function slice_opt = preAddingSlice(this, slice_opt)
        %         end
        
        function sl = onAddingSlice(this, sl)          
            if isempty(sl)
                this.RemoveSlice(sl);
                sl = sl.empty;
            else
                this.pending_slices.Add(sl);
                output = this.optimizeResourcePriceNew([], this.pending_slices{:});
                this.pending_slices.Clear();
                global g_results; %#ok<TLEV>
                % At the beginning, the slice is added, without consideration of
                % reconfiguration cost.
                stat = sl.get_reconfig_stat();
                stat.Profit = output.Profit(1);
                g_results = stat;       % The first event.
            end
        end
        %%
        % |finalize| should only be called when dimensiong network slices.
        function finalize(this, node_price, link_price, sub_slices)
            if nargin <= 3
                sub_slices = this.slices;
            end
            finalize@CloudNetwork(this, node_price, link_price, sub_slices);
            %             finalize@DynamicNetwork(this, sub_slices);
        end
    end
    
    methods(Static, Access = protected)
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
            if slice.getOption('Adhoc')==false || ~slice.isAdhocFlow
                ft = createflow@DynamicNetwork(this, slice, numflow);   
                ft{:,'Type'} = FlowType.Normal;
                if nargout >= 2
                    phy_adjacent = [];
                end
                return;
            end
            graph = this.graph;
            slice_opt.FlowPattern = slice.options.FlowPattern;
            slice_opt.DelayConstraint = slice.options.DelayConstraint;
            slice_opt = this.updateDynamicSliceOptions(slice, slice_opt);
            slice_opt.NumberFlows = numflow;
            slice_opt.NumberPaths = slice.options.NumberPaths;
            slice_opt.Method = slice.options.Method;
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
                    find(this.getNodeField('Capacity', new_phy_node_id) > 0);
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
                props.Weight = this.getLinkField('Weight', ...
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
                
                this.updateRedimensionCost(slice);
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
            node_id = slice.getDCNI;
            link_load = this.getLinkField('Load', link_id);
            node_load = this.getNodeField('Load', node_id);
            zero_load_link_id = link_load==0;
            zero_load_node_id = node_load==0;
            link_load(zero_load_link_id) = this.getLinkField('Capacity', ...
                link_id(zero_load_link_id)) * (1/20);
            node_load(zero_load_node_id) = this.getNodeField('Capacity', ...
                node_id(zero_load_node_id)) * (1/20);
            link_price = this.getLinkField('Price', link_id);
            node_price = this.getNodeField('Price', node_id);
            [~, slice.VirtualLinks.ReconfigCost] = ...
                slice.fcnLinkPricing(link_price, link_load);
            slice.VirtualLinks.ReconfigCost = ...
                DynamicSlice.THETA * slice.VirtualLinks.ReconfigCost;
            [~, slice.VirtualDataCenters.ReconfigCost] = ...
                slice.fcnNodePricing(node_price, node_load);
            slice.VirtualDataCenters.ReconfigCost = ...
                DynamicSlice.THETA * slice.VirtualDataCenters.ReconfigCost;
        end
        %         function [link_reconfig_cost, node_reconfig_cost] = ...
        %                 updateRedimensionCost(this, slice, link_id, node_id)
        %             global DEBUG; %#ok<NUSED>
        %             if nargin <= 2
        %                 link_id = slice.VirtualLinks.PhysicalLink;
        %             end
        %             if nargin <= 3
        %                 node_id = slice.getDCNI;
        %             end
        %             link_load = this.getLinkField('Load', link_id);
        %             node_load = this.getNodeField('Load', node_id);
        %             zero_load_link_id = link_load==0;
        %             zero_load_node_id = node_load==0;
        %             link_load(zero_load_link_id) = this.getLinkField('Capacity', ...
        %                 link_id(zero_load_link_id)) * (1/20);
        %             node_load(zero_load_node_id) = this.getNodeField('Capacity', ...
        %                 node_id(zero_load_node_id)) * (1/20);
        %             link_price = this.getLinkField('Price', link_id);
        %             node_price = this.getNodeField('Price', node_id);
        %             [~, link_reconfig_cost] = slice.fcnLinkPricing(link_price, link_load);
        %             [~, node_reconfig_cost] = slice.fcnNodePricing(node_price, node_load);
        %         end
    end
end

