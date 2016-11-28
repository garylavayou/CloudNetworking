%% Physical Network

%%
classdef PhysicalNetwork < handle
    
    properties (Constant)
        AF0 = 1;      % 1 unit of processing resource per 1Mbps data
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
    
    properties
        
        Topology;
        slices;
        VNFTable;            % Meta data of virtual network function
        Delta;               % the proportion of link when computing the network ultilization.
        
        LinkOption;          % Options for link properties
        NumberNodes;         % Number of physical nodes
        NumberLinks;         % Number of physical links
        NumberSlices;        % Number of network slices
        NumberVNFs;          % Number of VNF types
        NumberPaths;         % Number of candidate paths of all slices
    end
    properties (Access = private)
        graph;
        % mapping link index to edge table index.
        %         edge_table_index;        
    end
    
    methods 
        function this = PhysicalNetwork(graph_data, VNF_data, slice_opt, delta, beta)
            % Store graph information
            %     PhysicalNetwork(graph_data, VNF_opt, slice_opt)
            
            this.graph = DirectedGraph(graph_data);
            this.Topology = graph_data;
            %% Edge Index Mapping
            % In the edge table, link is indexed by rows, which is different from the
            % default index scheme of matlab matrix. So we add a column-index to the Edge
            % Table. 
            [s,t] = graph_data.findedge;
            idx = this.graph.IndexEdge(s,t);
            this.Topology.Edges.Index = idx;
            this.Topology.Edges.Properties.VariableDescriptions{4} = ...
                'Index the edges by column.';
%             [~, this.edge_table_index] = sort(this.Topology.Edges.Index,'ascending');

            % Initialize VNF Specification
            if isfield(VNF_data, 'StaticCost')
                this.Topology.Nodes.StaticCost = VNF_data.StaticCost;
            else
                average_node_capacity = min(this.Topology.Nodes.Capacity);
                switch VNF_data.Model
                    case VNFIntegrateModel.AllInOne
                        rng(VNF_data.RandomSeed(1));
                        this.Topology.Nodes.StaticCost = ...
                            round(average_node_capacity./randi([10 20], [this.NumberNodes, 1]), -2);
                    case VNFIntegrateModel.SameTypeInOne
                        
                    case VNFIntegrateModel.Separated
                        
                end
            end
            this.VNFTable = table;
            if isfield(VNF_data, 'ProcessEfficiency')
                this.VNFTable.ProcessEfficiency = VNF_data.ProcessEfficiency;
            else
                rng(VNF_data.RandomSeed(2));
                this.VNFTable.ProcessEfficiency = 0.5 + rand([VNF_data.Number, 1]);
            end
            
            % Add slice data
            this.NumberSlices = 0;
            this.slices = cell(length(slice_opt),1);
            for i = 1:length(slice_opt)
                this.AddSlice(slice_opt(i));
                this.slices{i}.Parent = this;
            end
            
            this.AllocateFlowId;
            this.AllocatePathId;
            for i = 1:length(slice_opt)
                this.slices{i}.initializeState;
            end
            
            if nargin >=4 && ~isempty(delta)
                this.Delta = delta;
            else
                this.Delta = 0.5;
            end
            if nargin >=5
                this.setLinkField('Beta', beta{1});
                this.setNodeField('Beta', beta{2});
            else
                this.setLinkField('Beta', 1*this.getLinkField('Capacity'));
                this.setNodeField('Beta', 1*this.getNodeField('Capacity'));
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
    end
    
    methods 
        
        function AllocateFlowId(this)
            % Allocate flow identifier
            flow_id = 0;
            for i = 1:this.NumberSlices
                slice = this.slices{i};
                slice.FlowTable.Identifier = flow_id + (1:height(slice.FlowTable));
                flow_id = flow_id + height(slice.FlowTable);
            end
        end
        
        function AllocatePathId(this)
            % Allocate path identifier
            path_id = 0;
            for i = 1:this.NumberSlices
                for j = 1:height(this.slices{i}.FlowTable)
                    path_list = this.slices{i}.FlowTable.Paths(j).paths;
                    for k = 1:length(path_list)
                        path_id = path_id + 1;
                        path_list{k}.id = path_id;
                    end
                end
            end
        end
        
        function idx = LinkId(this, s, t)
            % LinkId  get the index of links by the head and tail nodes' id.
            if nargin == 1
                idx = this.graph.IndexEdge;
            elseif nargin == 2
                idx = this.graph.IndexEdge(s);
            else
                idx = this.graph.IndexEdge(s,t);
            end
        end
    
        function setLinkField(this, name, value)
            % setLinkPrice
            %    the link index mapping between DirectedGraph and digraph is as 
            %    follows.
            if strcmp(name,'Index')
                error('Index cannot be set from the outside.');
            end
            this.Topology.Edges{:,{name}} = value(this.Topology.Edges.Index);
        end
        
        function value = getLinkField(this, name, link_id)
            % getLinkPrice
            %    link price corresponds to column indexed links.
            if nargin == 2
                % v(x,1) makes it column vector.
                value(this.Topology.Edges.Index,1) = this.Topology.Edges{:,{name}};
            elseif nargin >= 3
                value(this.Topology.Edges.Index,1) = this.Topology.Edges{:,{name}};
                value = value(link_id);
            else
                error('input arguments are not enough.');
            end
        end
        
        function setNodeField(this, name, value)
            this.Topology.Nodes{:,{name}} = value;
        end
        
        function value = getNodeField(this, name, node_id)
            if nargin == 2
                value = this.Topology.Nodes{:,{name}};
            elseif nargin >= 3
                value = this.Topology.Nodes{node_id, {name}};
            else
                error('input arguments are not enough.');
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
        
        function c = staticNodeCost(this)
            c = mean(this.getNodeField('StaticCost'));
        end
        
        function V = totalNodeCapacity(this)
            V = sum(this.getNodeField('Capacity'));
        end
        
        function C = totalLinkCapacity(this)
            C = sum(this.getLinkField('Capacity'));
        end
        
        function theta = networkUtilization(this, node_load, link_load)
            if nargin == 1
                node_load = this.getNodeField('Load');
                link_load = this.getLinkField('Load');
            end
            theta_v = sum(node_load)/sum(this.getNodeField('Capacity'));
            theta_l = sum(link_load)/sum(this.getLinkField('Capacity'));
            theta = this.Delta*theta_v + (1-this.Delta)*theta_l;
        end
        
        optimizeResourcePrice(this, init_link_price, init_node_price);
        [single_slice, rate, utility] = singleSliceOptimization( this );
        welfare = optimizeNetSocialWelfare( this );
        welfare = optimizeNetSocialWelfare2( this );
        welfare = optimizeNetSocialWelfare3( this );
        welfare = optimizeNetSocialWelfare4( this );
        
        AddSlice(this, slice_opt);       
        plot(this);
        
        %         function RemoveSlice(this)
        %         end
        % generate slice data
        %         slice_data = LoadSliceData(this, model);
        
    end
    methods(Access=private)
        function saveStates(this)
            for i = 1:this.NumberSlices
                this.slices{i}.Variables.x = this.slices{i}.x_path;
                this.slices{i}.Variables.z = this.slices{i}.z_npf;
                this.slices{i}.VirtualNodes.Load = this.slices{i}.getNodeLoad;
                this.slices{i}.VirtualLinks.Load = this.slices{i}.getLinkLoad;
                this.slices{i}.FlowTable.Rate = this.slices{i}.getFlowRate;
                this.slices{i}.setPathBandwidth;
            end
        end
        function clearStates(this)
            for i = 1:this.NumberSlices
                this.slices{i}.VirtualNodes.Load = zeros(this.slices{i}.NumberVirtualNodes,1);
                this.slices{i}.VirtualLinks.Load = zeros(this.slices{i}.NumberVirtualLinks,1);
                this.slices{i}.FlowTable.Rate = zeros(this.slices{i}.NumberFlows,1);
                this.slices{i}.setPathBandwidth(zeros(this.slices{i}.NumberPaths,1));
                this.slices{i}.Variables.x = [];
                this.slices{i}.Variables.z = [];
            end
        end
    end
    methods (Static)
        graph_data = LoadNetworkData(model, link_opt, node_opt);
        
        function dt = LinkDelay(delay_opt, bandwidth)
            % link bandwidth-delay convert function
            %     dt = PhysicalNetwork.LinkDelay(delay_opt, bandwidth) convert bandwidth
            %     to delay.
            if delay_opt == LinkDelayOption.BandwidthPropotion
                dt = 0.001*bandwidth;
            elseif delay_opt == LinkDelayOption.BandwidthInverse
                dt = 100./bandwidth;
            else
                error('the delay option (%s) cannot be handled.', delay_opt.char);
            end
        end
        
        [profit, grad] = fcnNetWelfare(x_vars, S);
        hess = fcnHessian(x_vars, ~, S);
    end
    
end

%% Properties
% * *Topology*: including a NodeTable |Nodes| and an EdgeTable |Edges|.
%
% |Nodes|: the fields in node table include _Name_, _Location_, _Capacity_, _StaticCost_,
% _Load_, _Price_.
%
% |Edges|: the fields in the edge table include _EndNodes_, _Weight_, _Capacity_, _Index_,
% _Load_, _Price_.  
%
% * *NumberNodes*
% * *NumberLinks* |get|
% * *NumberSlices*
% * *VNFTable*: including following fields.
%
% _ProcessEfficiency_: The coefficient for converting service rate to processing resource
%  requirement, i.e. $ProcessLoad = ServiceRate \times ProcessEfficiency$;

%% Methods
% * *LoadNetworkData* |static| : generate graph data.
% 
%       graph_data = LoadNetworkData(model, link_opt, node_opt)
%
% * *LoadVNFData* |static| : Generate virtual network function data.
% 
%       VNF_data = LoadVNFData(this, number_VNF, VNF_model)
%
% * *LinkDelay* |static| : convert bandwidth to link delay.
% 
%       dt = LinkDelay(delay_opt, bandwidth)
%
% |delay_opt|: enumeration type of LinkDelayOption.
%
% |bandwidth|: bandwidth with unit of |Mbps|.
%
% |dt|: delay with unit of |ms|.
%
% * *LinkId*: link index of links in the edge table
%
%       idx = LinkId(this, s, t)
%
% |s|: source nodes of links;
%
% |t|: tail nodes of links;
%
% * *setLinkField*: set the value for a field in the Edge Table. Link price corresponds to
% column indexed links, so this method is used when a column-indexed data value is
% provided.
%
%      setLinkField(this, name, value)
%
% |name|: a character array represents the field name.
%
% |value|: a column vector stores values to be set for the target field.
%
% * *getLinkField*: get the value from a field in the Edge Table. See also setLinkField.
%
%      value = getLinkField(this, name)
%
% * *AddSlice* : Add Slice to Substrate Network
%
%      AddSlice(phy_network, slice_opt)
% 
% |slice_opt|:  option for the added slice;
%
% * *AllocateFlowId* : Allocate flow identifier.
%
%      AllocateFlowId(phy_network)
%
% * *AllocatePathId* : Allocate path identifier
%
%      AllocatePathId(phy_network)
%
% * *plot* : Visualize Substrate Network and Network Slices
%
%      plot(phy_network)