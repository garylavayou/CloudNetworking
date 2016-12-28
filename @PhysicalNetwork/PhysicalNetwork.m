%% Physical Network

%%
classdef PhysicalNetwork < handle
    
    %% Properties

    properties (SetAccess = private)
        %%%
        % * *Topology*: including a NodeTable |Nodes| and an EdgeTable |Edges|.
        %
        % |Nodes|: the fields in node table include _Name_, _Location_, _Capacity_, _StaticCost_,
        % _Load_, _Price_.
        %
        % |Edges|: the fields in the edge table include _EndNodes_, _Weight_, _Capacity_, _Index_,
        % _Load_, _Price_.
        Topology;
        graph;
        slices;
        VNFTable;            % Meta data of virtual network function
        
        LinkOption;          % Options for link properties
        NumberNodes;         % Number of physical nodes
        NumberLinks;         % Number of physical links
        NumberSlices;        % Number of network slices
        %%%
        % * *VNFTable*: including following fields.
        %
        % _ProcessEfficiency_: The coefficient for converting service rate to processing resource
        %  requirement, i.e. $ProcessLoad = ServiceRate \times ProcessEfficiency$;
        NumberVNFs;          % Number of VNF types
        NumberPaths;         % Number of candidate paths of all slices
        NumberFlows;         % Number of flows of all slices
    end
    properties
        delta;          % the proportion of link when computing the network ultilization.
        slice_template;
    end
    properties (Access = private)
        lambda;
    end
    methods
        %% Constructor
        function this = PhysicalNetwork(node_opt, link_opt, VNF_opt, options)
            % Store graph information
            %     PhysicalNetwork(node_opt, link_opt, VNF_opt, options)
            this.Topology = PhysicalNetwork.loadNetworkData(node_opt, link_opt);
            this.graph = DirectedGraph(this.Topology);
            %%%
            % *Edge Index Mapping*:
            % In the edge table, link is indexed by rows, which is different from the
            % default index scheme of matlab matrix. So we add a column-index to the Edge
            % Table.
            [s,t] = this.Topology.findedge;
            idx = this.graph.IndexEdge(s,t);
            this.Topology.Edges.Index = idx;
            this.Topology.Edges.Properties.VariableDescriptions{4} = ...
                'Index the edges by column.';
            this.setLinkField('Load', 0);
            this.setNodeField('Load', 0);
            
            % Initialize VNF Specification
            switch VNF_opt.Model
                case VNFIntegrateModel.AllInOne
                    if VNF_opt.StaticCostOption == NodeStaticCostOption.Random
                        if ~isfield(VNF_opt, 'static_cost_range') ||...
                                isempty(VNF_opt.static_cost_range)
                            scr = [0.1 0.3];
                        else
                            scr = VNF_opt.static_cost_range;
                        end
                        rng(VNF_opt.RandomSeed(1));
                        node_capacity = this.Topology.Nodes.Capacity;
                        node_uc = this.Topology.Nodes.UnitCost;
                        avg_static_cost = mean(node_capacity.*node_uc)...
                            .*(rand([this.NumberNodes, 1])*(scr(2)-scr(1))+scr(1));
                        this.Topology.Nodes.StaticCost = ...
                            round(avg_static_cost, -2);
                    elseif VNF_opt.StaticCostOption == NodeStaticCostOption.NetworkSpecified
                        % TODO
                    end
                case VNFIntegrateModel.SameTypeInOne
                    % TODO
                case VNFIntegrateModel.Separated
                    % TODO
            end
            
            this.VNFTable = table;
            if isfield(VNF_opt, 'ProcessEfficiency')
                this.VNFTable.ProcessEfficiency = VNF_opt.ProcessEfficiency;
            else
                rng(VNF_opt.RandomSeed(2));
                this.VNFTable.ProcessEfficiency = 0.5 + rand([VNF_opt.Number, 1]);
            end
            
            % Add slice data
            this.NumberSlices = 0;
            %             if nargin >=5
            %                 this.setLinkField('Beta', beta{1});
            %                 this.setNodeField('Beta', beta{2});
            %             else
            %                 this.setLinkField('Beta', 1*this.getLinkField('Capacity'));
            %                 this.setNodeField('Beta', 1*this.getNodeField('Capacity'));
            %             end
            if nargin<=4
                options.delta = 0.5;
                %                 options.beta.node = 1;          % NOTE: this is not used temporarily
                %                 options.beta.link = 1;
            else
                if ~isfield(options, 'delta') || isempty(options.delta)
                    options.delta = 0.5;
                end
            end
            this.delta = options.delta;
        end
    end

%% Methods
    methods (Static)
        %%%
        % * *loadNetworkData* |static| : generate graph data.
        %
        %       graph_data = LoadNetworkData(link_opt, node_opt)
        %
        graph_data = loadNetworkData(node_opt, link_opt);
        %%%
        % * *LinkDelay* |static| : convert bandwidth to link delay.
        %
        %       dt = PhysicalNetwork.LinkDelay(delay_opt, bandwidth)
        %
        % |delay_opt|: enumeration type of LinkDelayOption.
        %
        % |bandwidth|: bandwidth with unit of |Mbps|.
        %
        % |dt|: delay with unit of |ms|.
        %
        function dt = LinkDelay(delay_opt, bandwidth)
            if delay_opt == LinkDelayOption.BandwidthPropotion
                dt = 0.001*bandwidth;
            elseif delay_opt == LinkDelayOption.BandwidthInverse
                dt = 100./bandwidth;
            else
                error('the delay option (%s) cannot be handled.', delay_opt.char);
            end
        end
        
    end
    
    % property access functions
    methods
        function opt = get.LinkOption(this)
            opt = this.Topology.Edges.Properties.UserData{1};
        end
        function n = get.NumberNodes(this)
            n = this.Topology.numnodes;
        end
        function m = get.NumberLinks(this)
            m = this.Topology.numedges;
        end
        function n = get.NumberVNFs(this)
            n = height(this.VNFTable);
        end
        function n = get.NumberPaths(this)
            n = 0;
            for i = 1:this.NumberSlices
                n = n + this.slices{i}.NumberPaths;
            end
        end
        function n = get.NumberFlows(this)
            n = 0;
            for i = 1:this.NumberSlices
                n = n + this.slices{i}.NumberFlows;
            end
        end
    end
    
    methods
        
        %%%
        % * *AllocateFlowId* : Allocate flow identifier.
        %
        %      AllocateFlowId(phy_network, start_slice)
        %
        % |start_slice|: the slice index or the handle of the slice.
        function AllocateFlowId(this, start_slice)
            % Allocate flow identifier
            if nargin <= 1
                slice_id = 1;
            else
                if isa(start_slice, 'Slice')
                    slice_id = this.findSlice(start_slice);
                else
                    slice_id = start_slice;
                end
            end
            if isempty(slice_id) || slice_id == 1
                flow_id = 0;
            else
                flow_id = this.slices{slice_id-1}.FlowTable.Identifier(end);
            end
            for s = slice_id:this.NumberSlices
                slice = this.slices{s};
                slice.FlowTable.Identifier = flow_id + (1:height(slice.FlowTable))';
                flow_id = flow_id + height(slice.FlowTable);
            end
        end
        
        %%%
        % * *AllocatePathId* : Allocate path identifier
        %
        %      AllocatePathId(phy_network)
        function AllocatePathId(this, start_slice)
            if nargin <= 1
                slice_id = 1;
            else
                if isa(start_slice, 'Slice')
                    slice_id = this.findSlice(start_slice);
                else
                    slice_id = start_slice;
                end
            end
            if isempty(slice_id) || slice_id == 1
                path_id = 0;
            else
                path_list = this.slices{slice_id-1}.FlowTable.Paths(end);
                path_id = path_list.paths{end}.id;
            end
            for s = slice_id:this.NumberSlices
                for j = 1:height(this.slices{s}.FlowTable)
                    path_list = this.slices{s}.FlowTable.Paths(j).paths;
                    for k = 1:length(path_list)
                        path_id = path_id + 1;
                        path_list{k}.id = path_id;
                    end
                end
            end
        end    
        %%%
        % * *LinkId*: link index of links in the edge table
        %
        %       idx = LinkId(this, s, t)
        %
        % |s|: source nodes of links;
        %
        % |t|: tail nodes of links;
        function idx = LinkId(this, s, t)
            % LinkId  get the index of links by the head and tail nodes' id.
            switch nargin
                case 1
                    idx = this.graph.IndexEdge;
                case 2
                    idx = this.graph.IndexEdge(s);
                otherwise
                    idx = this.graph.IndexEdge(s,t);
            end
        end
        %%%
        % * *setLinkField*: set the value for a field in the Edge Table. Link price
        % corresponds to column indexed links, so this method is used when a
        % column-indexed data value is provided.
        %
        %      setLinkField(this, name, value)
        %
        % |name|: a character array represents the field name.
        %
        % |value|: a column vector stores values to be set for the target field.
        function setLinkField(this, name, value)
            if strcmp(name,'Index')
                error('Index cannot be set from the outside.');
            end
            if isscalar(value)
                this.Topology.Edges{:,{name}} = value;
            else
                this.Topology.Edges{:,{name}} = value(this.Topology.Edges.Index);
            end
        end
        %%%
        % * *getLinkField*: get the value from a field in the Edge Table. See also setLinkField.
        %
        %      value = getLinkField(this, name, link_id)
        %
        %    link field corresponds to column indexed links.
        % |link_id| can be integers or logical numbers.
        function value = getLinkField(this, name, link_id)
            if nargin < 2 || isempty(name)
                error('input arguments are not enough (name, link_id).');
            elseif nargin < 3
                link_id = 1:this.NumberLinks;
            end
            switch name
                case 'ResidualCapacity'
                    value(this.Topology.Edges.Index,1) = ...
                        this.Topology.Edges{:,{'Capacity'}} - ...
                        this.Topology.Edges{:,{'Load'}};
                otherwise
                    value(this.Topology.Edges.Index,1) = this.Topology.Edges{:,{name}};
            end
            value = value(link_id);
        end
        
        function setNodeField(this, name, value)
            this.Topology.Nodes{:,{name}} = value;
        end
        
        function value = getNodeField(this, name, node_id)
            if nargin < 2 || isempty(name)
                error('input arguments are not enough (name, node_id).');
            elseif nargin < 3
                node_id = 1:this.NumberNodes;
            end
            switch name
                case 'ResidualCapacity'
                    value = this.Topology.Nodes{node_id,{'Capacity'}} - ...
                        this.Topology.Nodes{node_id,{'Load'}};
                otherwise
                    value = this.Topology.Nodes{node_id, {name}};
            end
        end
        
        function c = getLinkCost(this, link_load)
            % compute the total link cost.
            if nargin == 1
                c = dot(this.getLinkField('Load'), this.getLinkField('UnitCost'));
            else
                c = dot(link_load, this.getLinkField('UnitCost'));
            end
        end
        
        function c = getNodeCost(this, node_load)
            % compute the total node cost.
            if nargin == 1
                c = dot(this.Topology.Nodes.Load, this.getNodeField('UnitCost'));
            else
                c = dot(node_load, this.getNodeField('UnitCost'));
            end
        end
        
        function c = unitStaticNodeCost(this)
            c = mean(this.getNodeField('StaticCost'));
        end
        
        function c = getStaticCost(this, node_load, link_load)
            if nargin <= 2
                link_load = this.getLinkField('Load');
            end
            if nargin <=1
                node_load = this.getNodeField('Load');
            end
            epsilon = this.unitStaticNodeCost;
            theta = this.utilizationRatio(node_load, link_load);
            c = epsilon*((this.NumberNodes-1)*theta+1);
        end
        
        %%% 
        % *Network Operation Cost*: 
        % two methods to calculate network cost.
        % # Calculate with the approximate model, where the static node cost is computed
        % by the approximate formula.
        % # Calculate with the accurate model, where the static node cost is computed by
        % the solution of VNF deployment.
        function c = getNetworkCost(this, node_load, link_load, model)
            if nargin <=1 || isempty(node_load)
                node_load = this.getNodeField('Load');
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.getLinkField('Load');
            end
            if nargin <= 3
                model = 'Approximate';
            end
            
            if strcmp(model, 'Approximate')
                c = this.getNodeCost(node_load) + this.getLinkCost(link_load) +...
                    this.getStaticCost(node_load, link_load);
            elseif strcmp(model, 'Accurate')
                c = this.getNodeCost(node_load) + this.getLinkCost(link_load);
                b_deployed = node_load > 0;
                c = c + sum(this.getNodeField('StaticCost', b_deployed));
            end
        end
        
        function V = totalNodeCapacity(this)
            V = sum(this.getNodeField('Capacity'));
        end
        
        function C = totalLinkCapacity(this)
            C = sum(this.getLinkField('Capacity'));
        end
        
        function theta = utilizationRatio(this, node_load, link_load)
            if nargin == 1
                node_load = this.getNodeField('Load');
                link_load = this.getLinkField('Load');
            end
            theta_v = sum(node_load)/this.totalNodeCapacity;
            theta_l = sum(link_load)/this.totalLinkCapacity;
            theta = this.delta*theta_v + (1-this.delta)*theta_l;
        end
        
        % Remove the slice with identifier |id|.
        function sl = RemoveSlice(this, id)
            for s = 1:this.NumberSlices
                if this.slices{s}.Identifier == id
                    sl = this.slices{s};
                    break;
                end
            end
            this.NumberSlices = this.NumberSlices - 1;
            this.slices(s) = [];
            this.AllocateFlowId(s);
            this.AllocatePathId(s);
        end
        
        %%%
        % findSlice
        % Find the given slice's Index.
        % Find the index of slices with Type |key|, return a row vector of index.
        function sid = findSlice(this, key, ~)
            if isa(key, 'Slice')
                for s = 1:this.NumberSlices
                    if key == this.slices{s}
                        sid = s;
                        return;
                    end
                end
                sid = [];
            else
                sid = false(1,this.NumberSlices);
                for s = 1:this.NumberSlices
                    if key == this.slices{s}.Type
                        sid(s) = true;
                    end
                end
                sid = find(sid);
            end
        end
        
        % type_index is a scalar.
        function [p,r] = statSlice(this, type_index, profit)
            s_index = this.findSlice(type_index);
            if isempty(s_index)
                p = {0, 0, 0};
                if nargout >= 2
                    r = {0, 0, 0};
                end
            else
                p = {mean(profit(s_index)), max(profit(s_index)), min(profit(s_index))};
                if nargout >= 2
                    sum_rate = 0;
                    num_flow = 0;
                    max_rate = 0;
                    min_rate = inf;
                    for s = s_index     % s_index is a row vector
                        sum_rate = sum(this.slices{s}.FlowTable.Rate) + sum_rate;
                        num_flow = this.slices{s}.NumberFlows + num_flow;
                        max_rate = max(max_rate, max(this.slices{s}.FlowTable.Rate));
                        min_rate = min(min_rate, min(this.slices{s}.FlowTable.Rate));
                    end
                    avg_rate = sum_rate/num_flow;
                    r = {avg_rate, max_rate, min_rate};
                end
            end
        end
    end
    methods
        [output] = optimizeResourcePrice(this, init_price, options);
        [output, single_slice] = singleSliceOptimization(this , options);
        [output] = resourcePartitionOptimization(this, slice_weight, options);
        [output] = optimizeNetSocialWelfare1( this, options );
        output = StaticSlicing(this, slice, options);
        welfare = optimizeNetSocialWelfare2( this, options );
        welfare = optimizeNetSocialWelfare2a( this );
        welfare = optimizeNetSocialWelfare3( this );
        welfare = optimizeNetSocialWelfare4( this );
        %%%
        % * *AddSlice* : Add Slice to Substrate Network
        %
        %      AddSlice(phy_network, slice_opt)
        %
        % |slice_opt|:  option for the added slice;
        AddSlice(this, slice_opt);
        %%%
        % * *plot* : Visualize Substrate Network and Network Slices
        %
        %      plot(phy_network)
        plot(this);
    end
    methods(Access=private)
        function saveStates(this, lambda)
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.Variables.x = sl.x_path;
                sl.Variables.z = sl.z_npf;
                sl.VirtualNodes.Load = sl.getNodeLoad;
                sl.VirtualLinks.Load = sl.getLinkLoad;
                sl.FlowTable.Rate = sl.getFlowRate;
                sl.setPathBandwidth;
                this.lambda = lambda;
            end
        end
        function clearStates(this)
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.VirtualNodes.Load = zeros(sl.NumberVirtualNodes,1);
                sl.VirtualLinks.Load = zeros(sl.NumberVirtualLinks,1);
                sl.FlowTable.Rate = zeros(sl.NumberFlows,1);
                sl.setPathBandwidth(zeros(sl.NumberPaths,1));
                sl.Variables.x = [];
                sl.Variables.z = [];
                this.lambda = [];
            end
        end
    end
    
    properties (Constant)
        EdgeColor = [ % Predefined set of color data, used to specify the edge color when visualize network.
            0     0     1;    %blue
            0     0.5   0;
            0     0.7   0.7;
            0.078 0.169 0.549;
            0     0.447 0.741;
            0.494 0.184 0.557;
            0.467 0.675 0.188;
            0.302 0.745 0.933;
            1     0.6   0.78;
            ];
        NodeColor = [ % Predefined set of color data, used to specify the node color when visualize network.
            1     0     0;  %red
            1     0     1;  %magenta
            0.75  0.75  0;
            0.871 0.5   0;
            0.851 0.325 0.098;
            0.25  0.25  0.25];
    end
    
end

