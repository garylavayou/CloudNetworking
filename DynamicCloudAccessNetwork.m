classdef DynamicCloudAccessNetwork < CloudAccessNetwork & DynamicNetwork
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function this = DynamicCloudAccessNetwork(node_opt, link_opt, VNF_opt, net_opt)
            this@CloudAccessNetwork(node_opt, link_opt, VNF_opt, net_opt);
            %    this@DynamicNetwork(node_opt, link_opt, VNF_opt, net_opt);
            % <DynamicNetwork> has the same field with <PhysicalNetwork> is treated as an
            % interface. Thus, through <CloudAccessNetwork> the initilization has
            % finished.
        end

        function sl = AddSlice(this, slice_opt, varargin)
            this.bypass = {'CloudNetwork'};  % or 'DynamicNetwork' is OK.
            if ~isfield(slice_opt, 'method') || isempty(slice_opt.method)
                slice_opt.method = 'dynamic-slicing';
            end
            AddSlice@CloudNetwork(this, slice_opt, varargin{:});
            sl = AddSlice@DynamicNetwork(this, slice_opt, varargin{:});
        end   
        
        function sl = RemoveSlice(this, arg1)
            sl = RemoveSlice@CloudNetwork(this, arg1);
            this.optimizeResourcePriceNew([], options);
        end
        
    end
    
    methods(Access=protected)
        function tf = onRemovingSlice(this) %#ok<MANU>
            tf = true;
        end
        function sl = createslice(this, slice_opt)
            this.slices{end+1} = DynamicAccessSlice(slice_opt);
            sl = this.slices{end};
        end
        
        function sl = onAddingSlice(this, sl)          
            %             output2 = this.singleSliceOptimization(this.options);
            %             vnf2 = sl.Vnf;
            output = this.optimizeResourcePriceNew([], this.options);
            %             vnf1 = sl.Vnf;
            if isempty(sl)
                this.RemoveSlice(sl);
                sl = sl.empty;
            else
                sl = this.slices{end};
                global g_results event_num; %#ok<TLEV>
                % At the beginning, the slice is added, without consideration of
                % reconfiguration cost.
                g_results.profit(event_num,1) = output.profit.AccuratePrice(1);
                g_results.solution(event_num,1) = sl.Variables;
                [g_results.cost(event_num,1),...
                    g_results.num_reconfig(event_num,1),...
                    g_results.rat_reconfig(event_num,1)]...
                    = sl.get_reconfig_cost;
                g_results.num_flows(event_num,1) = sl.NumberFlows;
            end
        end
        function finalize(this, node_price, link_price)
            finalize@CloudAccessNetwork(this, node_price, link_price);
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.setVnfCapacity;
                %                 if strcmp(this.options.PricingPolicy, 'quadratic-price')
                [~, sl.edge_reconfig_cost] = sl.fcnLinkPricing(...
                    sl.VirtualLinks.Price, sl.VirtualLinks.Load);
                sl.edge_reconfig_cost = DynamicSlice.THETA*sl.edge_reconfig_cost;
                % here the |node_price| is the price of all data centers.
                [~, sl.node_reconfig_cost] = sl.fcnNodePricing(...
                    sl.VirtualDataCenters.Price, sl.VirtualDataCenters.Load);
                sl.node_reconfig_cost = DynamicSlice.THETA*sl.node_reconfig_cost;
                %                 else
                %                     sl.edge_reconfig_cost = link_price / a;
                %                     sl.node_reconfig_cost = node_price / a;
                %                 end
            end
        end
    end
    
    methods(Static)
        %%%
        % subclass can override this method.
        % |slice|: type of DynamicAccessSlice
        function slice_opt = update_slice_option(slice, slice_opt)
            switch slice.Options.FlowPattern
                case FlowPattern.RandomInterDataCenter
                    slice_opt.NodeSet = slice.getDCNI;
                case FlowPattern.RandomInterBaseStation
                    slice_opt.NodeSet = slice.getBSNI;
                case FlowPattern.RandomDataCenter2BaseStation
                    slice_opt.BSNodeSet = slice.getBSNI;
                    slice_opt.DCNodeSet = slice.getDCNI;
                otherwise
                    error('error: cannot handle the flow pattern <%s>.', ...
                        slice.Options.FlowPattern.char);
            end
        end
    end
    methods(Access=protected)          
        %%%
        % *Create new flows*
        % Creating new flows in the slice could guarantee no extra node or link would be
        % needed. If we enable new flows from new locations, we should create the flow
        % in the network.
        % |ft|: return flow table entries.
        function ft = createflow(this, slice, numflow)
            % map virtual network to physical network
            A = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
            C = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
            for i = 1:slice.NumberVirtualLinks
                h = slice.Topology.Head(i);
                t = slice.Topology.Tail(i);
                ph = slice.VirtualNodes{h, 'PhysicalNode'};
                pt = slice.VirtualNodes{t, 'PhysicalNode'};
                A(ph, pt) = slice.Topology.Adjacent(h,t); %#ok<SPRIX>
                C(ph, pt) = slice.Topology.Capacity(h,t); %#ok<SPRIX>
            end
            graph = DirectedGraph(A, C);
            slice_opt.Pattern = slice.Options.FlowPattern;
            slice_opt.DelayConstraint = slice.Options.DelayConstraint;
            slice_opt.MiddleNodes = slice.getDCNI;
            switch slice.Options.FlowPattern
                case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
                    slice_opt.NodeSet = slice.VirtualNodes.PhysicalNode;
                case FlowPattern.RandomInterDataCenter
                    slice_opt.NodeSet = slice.getDCNI;
                case FlowPattern.RandomInterBaseStation
                    slice_opt.NodeSet = slice.getBSNI;
                case FlowPattern.RandomDataCenter2BaseStation
                    slice_opt.BSNodeSet = slice.getBSNI;
                    slice_opt.DCNodeSet = slice.getDCNI;
                otherwise
                    error('error: cannot handle the flow pattern <%s>.', ...
                        slice.Options.FlowPattern.char);
            end
            if nargin <= 2
                slice_opt.NumberFlows = 1;
            else
                slice_opt.NumberFlows = numflow;
            end
            slice_opt.NumberPaths = slice.Options.NumberPaths;
            slice_opt.method = 'dynamic-slicing';
            b_vailid_flow = false;
            while ~b_vailid_flow
                try
                    b_vailid_flow = true;
                    ft = this.generateFlowTable(graph, slice_opt);
                catch ME
                    disp(ME)
                    b_vailid_flow = false;
                end
            end
            %%%
            % Update slice information as creating slice.
            %% TODO
            % When new nodes/edges should be added.
            ft.Properties.VariableNames = ...
                {'Source', 'Target', 'Rate', 'Delay', 'Paths'};
            for k = 1:height(ft)
                path_list = ft{k,{'Paths'}};
                for p = 1:path_list.Width
                    path = path_list.paths{p};
                    path.node_list = slice.PhyscialNodeMap{path.node_list,'VirtualNode'};
                    path.id = this.path_identifier_generator.next;
                end
            end
        end
    end
        
end

