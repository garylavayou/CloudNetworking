%% Physical Network
% network resource abstraction and allocation description; No resource allocation method
% defined, which can be realized by subclasses, see also <CloudNetwork>.
classdef PhysicalNetwork < INetwork
	properties
		DataCenters;         % the forwarding node mapping of data center.
		% data_center(dc_id) returns the physical node id;
		VNFTable = table;    % Meta data of virtual network function
		%%%
		% * *VNFTable*: including following fields.
		%
		% _ProcessEfficiency_: The coefficient for converting service rate to processing resource
		%  requirement, i.e. $ProcessLoad = ServiceRate \times ProcessEfficiency$;
		slices;              % a list of <Slice> objects
		%%%
		% * *topo*: including a NodeTable |Nodes| and an EdgeTable |Edges|.
		%
		% |Nodes|: the fields in node table include _Name_, _Location_, _Capacity_,
		% _StaticCost_, _Load_, _Price_.
		%
		% |Edges|: the fields in the edge table include _EndNodes_, _Weight_, _Capacity_,
		% _Index_, _Load_, _Price_.
		AggregateLinkUsage;
		AggregateNodeUsage;
		slice_template;
	end
	
	%% Properties
	properties (Dependent)
		NumberSlices;        % Number of network slices
		NumberDataCenters;   % Number of data centers
		NumberFlows;         % Number of flows of all slices
		NumberVNFs;          % Number of VNF types
		
		LinkOptions;         % Options for link properties
		NodeOptions;         % Options for node properties
	end
	
	properties (Access = protected)
		flow_identifier_generator;
		path_identifier_generator;
	end
	
	methods
		%% Constructor
		%   PhysicalNetwork(node_opt, link_opt, VNF_opt, options)
		% TODO: alternativelt input arguments by struct.
		function this = PhysicalNetwork(varargin)
			if nargin >= 2
				node_opt = varargin{1};
				link_opt = varargin{2};
				netdata.topo = PhysicalNetwork.loadNetworkData(node_opt, link_opt);
				args = {netdata};
			else
				% args = cell(0);
				return;
			end
			this@INetwork(args{:});
			% 					if isempty(varargin)
			% 						return;
			% 					end
			
			VNF_opt = varargin{3};
			this.Nodes = netdata.topo.Nodes;
			this.Links = netdata.topo.Links;
			dc_node_index = find(this.Nodes.Capacity>0);
			this.DataCenters = table(dc_node_index, 'VariableName', {'Node'});
			c = 1;
			while c<=width(this.Nodes)
				name = this.Nodes.Properties.VariableNames{c};
				switch name
					case {'Name', 'Location'}
						c = c + 1;
					otherwise
						this.DataCenters{:,{name}} = this.Nodes{dc_node_index, {name}};
						this.Nodes(:,{name}) = [];
				end
			end
			this.Nodes.DataCenter = zeros(this.NumberNodes,1);
			this.Nodes{dc_node_index, 'DataCenter'} = (1:height(this.DataCenters))';
			this.writeLink('Load', 0);
			this.writeDataCenter('Load', 0);
			
			% Initialize VNF Specification
			this.initializeVNF(VNF_opt);
			
			this.flow_identifier_generator = SerialNumber(1, [], true);
			this.path_identifier_generator = SerialNumber(1, [], true);
			
			this.AggregateLinkUsage = zeros(this.NumberLinks,1);
			this.AggregateNodeUsage = zeros(this.NumberNodes,1);
		end
		
		function delete(this)
			for i = 1:length(this.slices)
				delete(this.slices{i});
			end
			delete@INetwork(this);
		end
	end
	
	methods (Access = protected)
		function newobj = copyElement(this)
			%% Make a shallow copy of all properties
			newobj = copyElement@INetwork(this);
			%% Deep Copy
			% *slice.Parent*: update this link to the copyed network object.
			for i = 1:this.NumberSlices
				newobj.slices{i} = this.slices{i}.copy;
				newobj.slices{i}.Parent = newobj;
			end
		end
	end
	
	%% Methods
	methods
		setOptions(this, opt_name, opt_value);
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
		%%
		% called by generateFlowTable.
		function info = updateDemandInfo(this, slice_opt)
			switch slice_opt.FlowPattern
				case FlowPattern.RandomSingleFlow
					info.NumberFlows = 1;
				case FlowPattern.RandomMultiFlow
					info.NumberFlows = this.NumberNodes*(this.NumberNodes-1);
				otherwise
					error('error: [%s] unidentified flow pattern <%s>.', ...
						calledby(0), slice_opt.FlowPattern.char);
			end
			if isfield(slice_opt, 'NodeSet')
				info.NumberNodes = length(slice_opt.NodeSet);
				info.NodeSet = slice_opt.NodeSet;
			else
				info.NumberNodes = this.NumberNodes;
				info.NodeSet = 1:info.NumberNodes;
			end
		end
		
		%%%
		% called by generateFlowTable.
		function options = updatePathConstraints(this, slice_opt)
			options = getstructfields(slice_opt, {'DelayConstraint'}, 'ignore');
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
			opt = this.topo.Edges.Properties.UserData{1};
		end
		
		function opt = get.NodeOptions(this)
			opt = this.topo.Nodes.Properties.UserData{1};
		end
		
		function n = get.NumberSlices(this)
			n = length(this.slices);
		end
		
		function n = get.NumberDataCenters(this)
			n = height(this.DataCenters);
		end
				
		function n = get.NumberVNFs(this)
			n = height(this.VNFTable);
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
					argout_1 = this.Links.EndNodes(:,1);
					argout_2 = this.Links.EndNodes(:,2);
				case 2
					argout_1 = this.Links.EndNodes(argin_1,1);
					argout_2 = this.Links.EndNodes(argin_1,2);
				otherwise
					argout_1 = this.graph.IndexEdge(argin_1,argin_2);
			end
		end
		%%%
		% * *writeLink*: set the value for a field in the Edge Table. 
		%
		%      writeLink(this, name, value)
		%
		% |name|: a character array represents the field name.
		%
		% |value|: a column vector stores values to be set for the target field.
		function writeLink(this, name, value)
			this.Links{:,{name}} = value;
		end		
		%%%
		% * *readLink*: get the value from a field in the Edge Table. See also
		% writeLink.
		%
		%      value = readLink(this, name, link_id)
		%
		%    link field corresponds to column indexed links.
		% |link_id| can be integers or logical numbers.
		function value = readLink(this, name, link_id)
			if nargin < 2 || isempty(name)
				error('input arguments are not enough (name, link_id).');
			end
			if nargin <= 2
				link_id = 1:this.NumberLinks;
			end
			switch name
				case 'ResidualCapacity'
					value = this.Links{:,{'Capacity'}} - this.Links{:,{'Load'}};
				otherwise
					value = readLink@INetwork(this, name, link_id);
			end
			if nargin >= 3
				value = value(link_id);
			end
		end
		
		function writeNode(this, name, value, node_index)
			if nargin < 3 || isempty(name) || isempty(value)
				error('input arguments are not enough (name, value).');
			end
			if nargin == 3
				node_index = 1:this.NumberNodes;
			end
			if strcmp(name, 'DataCenter')
				warning('does not support change data center configurations.');
			elseif contains(name, this.Nodes.Properties.VariableNames)
				this.Nodes{node_index,{name}} = value;
			end			
		end
		
		function writeDataCenter(this, name, value, dc_index)
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
		function value = readNode(this, name, node_id)
			if nargin < 2 || isempty(name)
				error('input arguments are not enough (name, node_id).');
			end
			if nargin < 3
				node_id = 1:this.NumberNodes;
			end
			if contains(name, this.Nodes.Properties.VariableNames)
					value = readNode@INetwork(this, name, node_id);
			elseif contains(name, this.DataCenters.Properties.VariableNames)
					value = zeros(this.NumberNodes, 1);
					value(this.DataCenters.Node) = this.readDataCenter(name);
					value = value(node_id);
					warning('[%s] read information from data center nodes.', calledby);
			else
				error('error:[%s] ''%s'' does not exist.', calledby, name);
			end
		end
		
		% * *dc_id*: data center index.
		function value = readDataCenter(this, name, dc_id)
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
				link_id = sl.Links.PhysicalLink;
				this.AggregateLinkUsage(link_id) = this.AggregateLinkUsage(link_id) - 1;
				node_id = sl.Nodes.PhysicalNode;
				this.AggregateNodeUsage(node_id) = this.AggregateNodeUsage(node_id) - 1;
				%                 if b_update
				%                     %%%
				%                     % Update the load of the substrate network.
				%                     % for calculate the residual capacity.
				%                     node_load = zeros(this.NumberNodes, 1);
				%                     node_load(sl.VirtualNodes.PhysicalNode) = sl.VirtualNodes.Load;
				%                     this.writeDataCenter('Load', node_load - node_load);
				%                     link_load = zeros(this.NumberLinks, 1);
				%                     link_load(sl.VirtualLinks.PhysicalLink) = sl.VirtualLinks.Load;
				%                     this.writeLink('Load', link_load - link_load);
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
	
	methods(Access=protected)
		function [cap, load] = readCapacityLoad(this)
			cap.node = this.readDataCenter('Capacity');
			cap.link = this.readLink('Capacity');
			load.node = this.readDataCenter('Load');
			load.link = this.readLink('Load');
		end
		%% getNetworkLoad
		% Network load equals to the sums of occupied capacity from all slices. See also
		% <Slice.getLinkCapacity>, <Slice.getNodeCapacity> and
		% <DynamicSlice.getLinkCapacity>, <DynamicSlice.getNodeCapacity>.
		%
		% |slices|: If the argument is provided, calculate load from the
		%						set of |slices|, otherwise, calculate from all slices.
		% |options|
		%			_Stage_
		%					'final': directly copy from 'Capacity' field of each slice.
		%					'temp': use the temporary variables of each slice, to
		%							calculate temporary capacity of each slice. Makesure
		%							the temporary variables is up-to-date, when calling
		%							this method.
		function load = getNetworkLoad(this, slices, options)
			defaultopts = struct('Stage', 'final');
			if nargin <= 2
				options = defaultopts;
			else
				options = structupdate(defaultopts, options);
			end
			
			if nargin <= 1 || isempty(slices)
				slices = this.slices;
			elseif ~iscell(slices)
				slices = {slices};
			end
			
			load.Node = zeros(this.NumberDataCenters, 1);
			load.Link = zeros(this.NumberLinks, 1);
			for i = 1:length(slices)
				sl = slices{i};
				link_id = sl.VirtualLinks.PhysicalLink;
				dc_id = sl.getDCPI;
				if sl.isFinal() || strcmpi(options.State, 'final')
					load.Node(dc_id) = load.Node(dc_id) + sl.getNodeCapacity;
					load.Link(link_id) = load.Link(link_id) + sl.getLinkCapacity;
				elseif strcmpi(options.State, 'temp')
					load.Node(dc_id) = load.Node(dc_id) + sl.getNodeCapacity(false);
					load.Link(link_id) = load.Link(link_id)+ sl.getLinkCapacity(false);
				else
					error('error:[%s] Invalid value for options ''Stage''=''%s''.', options.Stage);
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
