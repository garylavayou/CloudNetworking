%% Physical Network
% network resource abstraction and allocation description; No resource allocation method
% defined, which can be realized by subclasses, see also <CloudNetwork>. 
classdef PhysicalNetwork < matlab.mixin.Copyable 
    
    %% Properties
    properties (Dependent)
        NumberNodes;         % Number of physical nodes
        NumberLinks;         % Number of physical links
        NumberSlices;        % Number of network slices
        NumberDataCenters;   % Number of data centers
        %%%
        % * *VNFTable*: including following fields.
        %
        % _ProcessEfficiency_: The coefficient for converting service rate to processing resource
        %  requirement, i.e. $ProcessLoad = ServiceRate \times ProcessEfficiency$;
        NumberVNFs;          % Number of VNF types
        NumberPaths;         % Number of candidate paths of all slices
        NumberFlows;         % Number of flows of all slices

        LinkOptions;         % Options for link properties
        NodeOptions;         % Options for node properties
    end
    
    properties (GetAccess = public)
        %%%
        % * *Topology*: including a NodeTable |Nodes| and an EdgeTable |Edges|.
        %
        % |Nodes|: the fields in node table include _Name_, _Location_, _Capacity_,
        % _StaticCost_, _Load_, _Price_.
        %
        % |Edges|: the fields in the edge table include _EndNodes_, _Weight_, _Capacity_,
        % _Index_, _Load_, _Price_. 
        Topology;            % Class <digraph> object
        graph;               % Class <DirectedGraph> object
        DataCenters;         % the forwarding node mapping of data center.
                             % data_center(dc_id) returns the physical node id;
        VNFTable = table;    % Meta data of virtual network function
        slices;              % a list of <Slice> objects        
    end
    
    properties
        slice_template;
    end
    
    properties (Access = {?PhysicalNetwork, ?Slice})
        % |options| will be visited by <CloudNetwork> and <DynamicNetwork>
        options;
    end    

    properties (Access = protected)
        path_identifier_generator;
        flow_identifier_generator;
    end
    properties (Access = {?PhysicalNetwork,?VirtualNetwork})
        AggregateLinkUsage;
        AggregateNodeUsage;
    end
    
    methods
        %% Constructor
        %   PhysicalNetwork(node_opt, link_opt, VNF_opt, options)
        function this = PhysicalNetwork(varargin)
            if isempty(varargin)
                return;
            elseif isa(varargin{1}, 'PhysicalNetwork')
                this = varargin{1}.copy;
                return;
            elseif length(varargin)<3
                return;
            end
            node_opt = varargin{1};
            link_opt = varargin{2};
            VNF_opt = varargin{3};
            % Store graph information
            this.Topology = PhysicalNetwork.loadNetworkData(node_opt, link_opt);
            dc_node_index = find(this.Topology.Nodes.Capacity>0);
            this.DataCenters = table(dc_node_index, 'VariableName', {'NodeIndex'});
            c = 1;
            while c<=width(this.Topology.Nodes)
                name = this.Topology.Nodes.Properties.VariableNames{c};
                switch name
                    case {'Name', 'Location'}
                        c = c + 1;
                    otherwise
                        this.DataCenters{:,{name}} = ...
                            this.Topology.Nodes{dc_node_index, {name}};
                        this.Topology.Nodes(:,{name}) = [];
                end
            end
            this.Topology.Nodes.DataCenter = zeros(this.NumberNodes,1);
            this.Topology.Nodes{dc_node_index, 'DataCenter'} = ...
                (1:height(this.DataCenters))';            
            this.graph = DirectedGraph(this.Topology);
            %%%
            % *Edge Index Mapping*:
            % In the edge table, link is indexed by rows, which is different from the
            % default index scheme of matlab matrix. So we add a column-index to the Edge
            % Table.
            [s,t] = this.Topology.findedge;
            idx = this.graph.IndexEdge(s,t);
            this.Topology.Edges.Index = idx;
            this.Topology.Edges.Properties.VariableDescriptions{5} = ...
                'Index the edges by column in the adjacent matrix.';
            this.setLinkField('Load', 0);
            this.setDataCenterField('Load', 0);
            
            % Initialize VNF Specification
            this.initializeVNF(VNF_opt);   
            
            this.flow_identifier_generator = SerialNumber(1, [], true);
            this.path_identifier_generator = SerialNumber(1, [], true);
            
            this.AggregateLinkUsage = zeros(this.NumberLinks,1);
            this.AggregateNodeUsage = zeros(this.NumberNodes,1);
        end
        
        function delete(this)
            delete(this.graph);
            for i = 1:length(this.slices)
                delete(this.slices{i});
            end
        end
    end
    
    methods (Access = protected)
        % Child classes can overload these functions, which can be called from the
        % superclass, without consideration the specific class of the object. 
        function this = copyElement(pn)
            % Make a shallow copy of all properties
            this = copyElement@matlab.mixin.Copyable(pn);
            %% Deep Copy Issue
            % *slice.Parent*: update this link to the copyed network object.
            this.graph = pn.graph.copy;
            for i = 1:pn.NumberSlices
                this.slices{i} = pn.slices{i}.copy;
                this.slices{i}.Parent = this;
            end
        end        
    end
    
    %% Methods
    %%%
    % Subclass can override these functions
    methods (Abstract, Access = protected)
        %%%
        % Called by <AddSlice>.
        sl = createslice(this, slice_opt, varargin);
        %%%
        % Called by <AddSlice>.
        slice_opt = preAddingSlice(this, slice_opt);
    end
    methods (Access = protected)
        %%%
        % Called by _Constructor_;
        function initializeVNF(this, VNF_opt)
            if isfield(VNF_opt, 'ProcessEfficiency')
                this.VNFTable.ProcessEfficiency = VNF_opt.ProcessEfficiency;
            else
                if isfield(VNF_opt, 'RandomSeed')
                    rng(VNF_opt.RandomSeed(1));
                else
                    rng(floor(now));
                    warning('random seed is not sepecifed for VNF, set as %d', floor(now));
                end
                this.VNFTable.ProcessEfficiency = 0.5 + rand([VNF_opt.Number, 1]);
            end
        end
        %%%
        % called by generateFlowTable.
        function info = updateDemandInfo(this, slice_opt)
            switch slice_opt.FlowPattern
                case FlowPattern.RandomSingleFlow
                    info.NumberFlows = 1;
                case FlowPattern.RandomMultiFlow
                    info.NumberFlows = this.NumberNodes*(this.NumberNodes-1);
                otherwise
                    error('error: cannot handle the flow pattern <%s>.', ...
                        slice_opt.FlowPattern.char);
            end
            if isfield(slice_opt, 'NodeSet')
                info.NumberNodes = length(slice_opt.NodeSet);
                info.NodeSet = slice_opt.NodeSet;
            else
                info.NumberNodes = this.Topology.numnodes;
                info.NodeSet = 1:info.NumberNodes;
            end
        end
        %%%
        % called by generateFlowTable.
        function options = updatePathConstraints(this, slice_opt)
            options = getstructfields(slice_opt, ...
                {'DelayConstraint'}, 'ignore');
            options.DelayModel = this.LinkOptions.DelayModel;
        end
        function end_points = generateEndPoints(this, info, slice_opt) %#ok<INUSL,INUSD>
            end_points = info.NodeSet(unique_randi(info.NumberNodes, 2, 'stable'));
        end
        %%%
        % Called by access AddSlice
        [flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt );      
        %%%
        % Called by _AddSlice_, return value is used by _generateFlowTable_;
        function graph = residualgraph(this)
            graph = this.graph;
        end
        
    end
    
    methods (Static)
        %%%
        % * *loadNetworkData* |static| : generate graph data.
        %
        %       graph_data = LoadNetworkData(link_opt, node_opt)
        %
        %% TODO: override it for subclass.
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
    methods (Static, Access=protected)
        function flag = assert_path_list(~, path_list)
            if isempty(path_list)
                flag = -1;
            else
                flag = 0;
            end
        end
    end
    % property access functions
    methods
        function opt = get.LinkOptions(this)
            opt = this.Topology.Edges.Properties.UserData{1};
        end
        
        function opt = get.NodeOptions(this)
            opt = this.Topology.Nodes.Properties.UserData{1};
        end
        
        function n = get.NumberNodes(this)
            n = this.Topology.numnodes;
        end
        
        function n = get.NumberSlices(this)
            n = length(this.slices);
        end
        
        function n = get.NumberDataCenters(this)
            n = height(this.DataCenters);
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
        
        setOptions(this, opt_name, opt_value);
    end

    methods    
        %%%
        % statistics of slices.
        function type_count = CountSlices(this)
            type_list = [this.slice_template.Type];
            type_count = zeros(1,length(type_list));
            type = zeros(1,this.NumberSlices);
            for i = 1:this.NumberSlices
                type(i) = this.slices{i}.Type;
            end
            [~, tid] = ismember(type, type_list);
            for t = tid
                assert(t~=0, 'error: unknown slice type.');
                type_count(t) = type_count(t) + 1;
            end
        end
        %%%
        % * *LinkId*: link index of links in the edge table
        %
        %      [head, tail] = LinkId(this, idx); 
        %      idx = LinkId(this, s, t)
        %
        % |s|: source nodes of links;
        %
        % |t|: tail nodes of links;
        %
        % |idx|: the numeric index of links.
        %
        % See Also _DirectedGraph.IndexEdge_.
        function [argout_1, argout_2] = LinkId(this, argin_1, argin_2)
            switch nargin
                case 1
                    [argout_1, argout_2] = this.graph.IndexEdge;
                case 2
                    [argout_1, argout_2] = this.graph.IndexEdge(argin_1);
                otherwise
                    argout_1 = this.graph.IndexEdge(argin_1,argin_2);
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
        % * *getLinkField*: get the value from a field in the Edge Table. See also
        % setLinkField. 
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
        
        function setNodeField(this, name, value, node_index)
            if nargin < 3 || isempty(name) || isempty(value)
                error('input arguments are not enough (name, value).');
            end
            if nargin == 3
                node_index = 1:this.NumberNodes;
            end
            this.Topology.Nodes{node_index,{name}} = value;
            
        end
        
        function setDataCenterField(this, name, value, dc_index)
            if nargin < 3 || isempty(name) || isempty(value)
                error('input arguments are not enough (name, value).');
            end
            if nargin == 3
                dc_index = 1:this.NumberDataCenters;
            end
            if strcmp(name, 'ResidualCapacity')
                error('error: ResidualCapacity cannot be set.');
            else
                this.DataCenters{dc_index,{name}} = value;
            end            
        end
        
        % * *node_id*: physical node index.
        function value = getNodeField(this, name, node_id)
            if nargin < 2 || isempty(name)
                error('input arguments are not enough (name, node_id).');
            end
            if nargin < 3
                node_id = 1:this.NumberNodes;
            end
            switch name
                case {'Name', 'Location', 'DataCenter'}
                    value = this.Topology.Nodes{node_id, {name}};
                otherwise
                    value = zeros(this.NumberNodes, 1);
                    value(this.DataCenters.NodeIndex) = this.getDataCenterField(name);
                    value = value(node_id);
            end
        end
        
        % * *dc_id*: data center index.
        function value = getDataCenterField(this, name, dc_id)
            if nargin < 2 || isempty(name)
                error('input arguments are not enough (name, node_id).');
            end
            if nargin < 3
                dc_id = 1:this.NumberDataCenters;
            end
            if strcmp('ResidualCapacity', name)
                value = this.DataCenters{dc_id,{'Capacity'}} - ...
                    this.DataCenters{dc_id,{'Load'}};
            else
                value = this.DataCenters{dc_id, {name}};
            end
        end
        
        function V = totalNodeCapacity(this)
            V = sum(this.getDataCenterField('Capacity'));
        end
        
        function C = totalLinkCapacity(this)
            C = sum(this.getLinkField('Capacity'));
        end
                
		function [r_mean, r_max, r_min, r_std] = nodeUtilization(this)
			node_load = this.getDataCenterField('Load');
			node_capacity = this.getDataCenterField('Capacity');
            % 			node_index = node_load > 1;
            %%%
            % Another method: ratio = sum(node_load)/sum(node_capacity);
			ratio = node_load ./ node_capacity;
			r_mean = mean(ratio);
            %%%
            % The range of node utilization may large, since the load of nodes depends on
            % the flow's location, the node's cost, and our objective is not to balancing
            % the node load. 
            if nargout >= 2
                r_max = max(ratio);
            end
            if nargout >= 3
                r_min = min(ratio);
            end
            if nargout >= 4
                r_std = std(ratio);
            end
		end
		
        function [r_mean, r_max, r_min, r_std] = linkUtilization(this)
            link_load = this.getLinkField('Load');
            link_capacity = this.getLinkField('Capacity');
            % link_index = link_load > 1;
            ratio = link_load ./ link_capacity;
            r_mean = mean(ratio);
            if nargout >= 2
                r_max = max(ratio);
            end
            if nargout >= 3
                r_min = min(ratio);
            end
            if nargout >= 4
                r_std = std(ratio);
            end
        end

        %%% 
        % * *RemoveSlice*:
        % Remove the slice with identifier |id|.
        % if no slice with identifier |id|, this method do not perform any operation.
        % |b_update| should be 1, if using static slicing method.
        function sl = RemoveSlice(this, arg1)
            if isnumeric(arg1)
                id = arg1;
                sid = this.findSlice(id, 'Identifier');
            elseif isa(arg1, 'Slice')
                sl = arg1;
                sid = this.findSlice(sl);
            end
%             if nargin <= 2
%                 b_update = false;
%             end
            if ~isempty(sid)
                sl = this.slices{sid};
                link_id = sl.VirtualLinks.PhysicalLink;
                this.AggregateLinkUsage(link_id) = this.AggregateLinkUsage(link_id) - 1;
                node_id = sl.VirtualNodes.PhysicalNode;
                this.AggregateNodeUsage(node_id) = this.AggregateNodeUsage(node_id) - 1;
%                 if b_update
%                     %%%
%                     % Update the load of the substrate network.
%                     % for calculate the residual capacity.
%                     node_load = zeros(this.NumberNodes, 1);
%                     node_load(sl.VirtualNodes.PhysicalNode) = sl.VirtualNodes.Load;
%                     this.setDataCenterField('Load', node_load - node_load);
%                     link_load = zeros(this.NumberLinks, 1);
%                     link_load(sl.VirtualLinks.PhysicalLink) = sl.VirtualLinks.Load;
%                     this.setLinkField('Load', link_load - link_load);
%                 end
                this.slices(sid) = [];
                %
                % Since the flow/path id might be used by other associate entities, the
                % cost to reallocate flow id is large. And reallocation is not necessary,
                % since the space of identifier is large enough.
                %    this.AllocateFlowId(sid);
                %    this.AllocatePathId(sid);
            else
                sl = [];
            end
        end
        
        %%%
        % * *findSlice*: 
        %   sid = findSlice(this, slice)
        %   sid = findSlice(this, key, field)
        % (1) Find the given slice's Index.
        % (2) Find the index of slices with Type |key|, return a row vector of index.
        function sid = findSlice(this, key, field)
            sid = [];
            if isa(key, 'Slice')
                for s = 1:this.NumberSlices
                    if key == this.slices{s}
                        sid = s;
                        return;
                    end
                end
            else
                if nargin <= 2
                    field = 'Type';
                end
                switch field
                    case 'Type'
                        sid = false(1,this.NumberSlices);
                        for s = 1:this.NumberSlices
                            if key == this.slices{s}.Type
                                sid(s) = true;
                            end
                        end
                        sid = find(sid);
                    case 'Identifier'
                        for s = 1:this.NumberSlices
                            if key == this.slices{s}.Identifier
                                sid = s;
                                break;
                            end
                        end
                    otherwise
                        warning('undefined key type [%s].', key);
                end
            end
        end
        
    end
    methods
        %%%
        % * *AddSlice* : Add Slice to Substrate Network
        %
        %      AddSlice(phy_network, slice_opt)
        %
        % |slice_opt|:  option for the added slice;
        sl = AddSlice(this, slice_opt, varargin);
        %%%
        % * *plot* : Visualize Substrate Network and Network Slices
        %
        %      plot(phy_network)
        plot(this, b_undirect);
    end
    methods(Access=protected)
        function allocatepathid(this, slice)
            for j = 1:height(slice.FlowTable)
                path_list = slice.FlowTable.Paths(j).paths;
                for k = 1:length(path_list)
                    path_list{k}.id = this.path_identifier_generator.next;
                end
            end
        end
                
        %% getNetworkLoad
        % Network load equals to the sums of occupied capacity from all slices. See also
        % <Slice.getLinkCapacity> and <Slice.getNodeCapacity>.
        %
        % If the 2nd argument is provided, calculate load from the set of |sub_slices|,
        % otherwise, calculate from all slices. 
        % If |option| is not provided, directly copy from 'Capacity' field of each slice.
        % Otheriwse if |option='sum'|, use the temporary variables |x| and |z| of
        % each slice, to calculate temporary capacity of each slice.  
        % Makesure the temporary variables is up-to-date , when calling this method.
        function [node_load, link_load] = getNetworkLoad(this, sub_slices, option)
            if nargin <= 1 || isempty(sub_slices)
                sub_slices = this.slices;
            elseif ~iscell(sub_slices)
                sub_slices = {sub_slices};
            end
            if nargin <= 2
                option = 'copy';
            end
            node_load = zeros(this.NumberDataCenters, 1);
            link_load = zeros(this.NumberLinks, 1);
            for i = 1:length(sub_slices)
                sl = sub_slices{i};
                link_id = sl.VirtualLinks.PhysicalLink;
                dc_id = sl.getDCPI;
                if strcmpi(option, 'copy')
                    node_load(dc_id) = node_load(dc_id) + sl.getNodeCapacity;
                    link_load(link_id) = link_load(link_id) + sl.getLinkCapacity;
                else  % 'sum'
                    node_load(dc_id) = node_load(dc_id) + sl.getNodeCapacity(sl.temp_vars.z);
                    link_load(link_id) = link_load(link_id)+ sl.getLinkCapacity(sl.temp_vars.x);                    
                end
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
    
    methods  % Deprecated
        %%%
        % * *AllocateFlowId* : Allocate flow identifier.
        %
        %      AllocateFlowId(phy_network, start_slice)
        %
        % |start_slice|: the slice index or the handle of the slice. If this argument is
        % not provided, we allocate flow id from the first slice.
        %         function id = AllocateFlowId(this, start_slice)
        %             % Allocate flow identifier
        %             if nargin <= 1
        %                 slice_id = 1;
        %             else
        %                 if isa(start_slice, 'Slice')
        %                     slice_id = this.findSlice(start_slice);
        %                 else
        %                     slice_id = start_slice;
        %                 end
        %             end
        %             if isempty(slice_id) || slice_id == 1
        %                 flow_id = uint64(0);
        %             else
        %                 flow_id = this.slices{slice_id-1}.FlowTable.Identifier(end);
        %             end
        %             for s = slice_id:this.NumberSlices
        %                 slice = this.slices{s};
        %                 % Since identifier space is large enough (64-bit), no need to worry that
        %                 % the identifier will duplicate.
        %                 slice.FlowTable.Identifier = flow_id + (1:height(slice.FlowTable))';
        %                 flow_id = flow_id + height(slice.FlowTable);
        %             end
        %         end
        
        %%%
        % * *AllocatePathId* : Allocate path identifier
        %
        %      AllocatePathId(phy_network)
        %         function AllocatePathId(this, start_slice)
        %             if nargin <= 1
        %                 slice_id = 1;
        %             else
        %                 if isa(start_slice, 'Slice')
        %                     slice_id = this.findSlice(start_slice);
        %                 else
        %                     slice_id = start_slice;
        %                 end
        %             end
        %             if isempty(slice_id) || slice_id == 1
        %                 path_id = uint64(0);
        %             else
        %                 path_list = this.slices{slice_id-1}.FlowTable.Paths(end);
        %                 path_id = path_list.paths{end}.id;
        %             end
        %             for s = slice_id:this.NumberSlices
        %                 for j = 1:height(this.slices{s}.FlowTable)
        %                     path_list = this.slices{s}.FlowTable.Paths(j).paths;
        %                     for k = 1:length(path_list)
        %                         path_id = path_id + 1;
        %                         path_list{k}.id = path_id;
        %                     end
        %                 end
        %             end
        %         end
    end
end