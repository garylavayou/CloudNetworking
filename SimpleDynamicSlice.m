classdef SimpleDynamicSlice < DynamicSlice
	%DynamicSlice Event-driven to dynamic configure slice.
	
	properties(Dependent)
		UnitReconfigureCost;
	end
	
	%% Constructor
	methods
		function this = SimpleDynamicSlice(slice_data)
      if nargin == 0
        args = {};
			else
				args = {slice_data};
      end
			this@DynamicSlice(args{:});
			this.op = SimpleDynamicSliceOptimizer(this, args{:});
		end
	end
	
	%% Public Methods
	methods		
		function op = getOptimizer(this, options)
			if nargin == 1
				this.op = SimpleDynamicSliceOptimizer(this);
			else
				this.op = SimpleDynamicSliceOptimizer(this, options);
			end
			op = this.op;
		end
		
		function finalize(this, prices)
			finalize@DynamicSlice(this, prices);
			if isfield(this.op.temp_vars,'v') && ~isempty(this.op.temp_vars.v)
				if this.NumberFlows > 0
					% Override the capacity setting in the super class.
					this.ServiceNodes.Capacity = full(sum(reshape(this.VNFCapacity, ...
						this.NumberServiceNodes,this.NumberVNFs),2));
				else
					warning('zero flow, Variables.v not initialized.');
				end
			elseif this.NumberFlows > 0
				this.op.Variables.v = this.getVNFCapacity;
			else
				warning('zero flow, Variables.v not initialized.');
			end
      if isfield(this.op.temp_vars,'c') && ~isempty(this.op.temp_vars.c)
        this.Links.Capacity = full(this.op.temp_vars.c);
      end
      
    end
		
		%%%
		% *Create new flows*
		% Creating new flows in the slice could guarantee no extra node or link would be
		% needed. If we enable new flows from new locations, we should create the flow
		% in the network.
		%         function ft = createflow(this, slice)
		%             graph = slice.Topology;
		%             switch slice.options.FlowPattern
		%                 case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
		%                     ft = this.generateFlowTable(graph,
		%                 case FlowPattern.RandomInterDataCenter
		%                 case FlowPattern.RandomInterBaseStation
		%                 case FlowPattern.RandomDataCenter2BaseStation
		%                 otherwise
		%                     error('error: cannot handle the flow pattern <%s>.', ...
		%                         slice.options.FlowPattern.char);
		%             end
		%         end
		
	end
	
	methods
		%% Get link and node capacity
		% <getLinkCapacity>
		% <getNodeCapacity>
		% Due to resource reservation or reconfiguration cost constraint, the link/node
		% capcity might be larger than the demand.
		function c = getLinkCapacity(this, isfinal)
			if nargin == 1 || isfinal
				c = this.Links.Capacity;
			else
				if this.invoke_method == 0
					c = getLinkCapacity@SimpleSlice(this, isfinal);
				else
					c = this.op.temp_vars.c;
				end
			end
		end
		
		function c = getNodeCapacity(this, isfinal) % isfinal=true by default
			if nargin == 1 || isfinal
				c = this.ServiceNodes.Capacity;
			else
				if this.invoke_method == 0
					c = getNodeCapacity@SimpleSlice(this, isfinal);
				else
					% since we do not reconfigure VNF capacity through fast slice reconfiguration,
					% the sum of VNF capacity equals to the node capacity.
					c = sum(reshape(this.temp_vars.v, ...
						this.NumberServiceNodes,this.NumberVNFs),2);
				end
			end
		end
		
	end
	
	methods
		%%%
		% At present this method inherit the superclass. See also <Slice.checkFeasible>.
		%         function b = checkFeasible(this, vars, options)
		%             if nargin <= 1
		%                 vars = [];
		%             else
		%                 vars = vars(1:this.NumberVariables);
		%             end
		%             if nargin <= 2
		%                 options = struct;
		%             end
		%             b = checkFeasible@Slice(this, vars, options);
		%             %             if b
		%             %                 % TODO add other check conditions.
		%             %                 switch options.Action
		%             %                     case 'add'
		%             %                     case 'remove'
		%             %                 end
		%             %             end
		%         end
	end
	
	methods		
		function stat = get_reconfig_stat(this, stat_names)
			if nargin == 1
				stat_names = {'All'};
			elseif ischar(stat_names)
				stat_names = {stat_names};
			end
			options = getstructfields(this.Parent.options, ...
				{'DiffNonzeroTolerance', 'NonzeroTolerance'});
			s = this.diffstate(true);
			
			stat = table;
			tol_vec = options.DiffNonzeroTolerance;
			for i = 1:length(stat_names)
				%% Time
				if contains(stat_names{i},{'All', 'Time'},'IgnoreCase',true)
					stat.Time = this.time.Current;
				end
				if contains(stat_names{i},{'All', 'Utilization'},'IgnoreCase',true)
					stat.Utilization = this.utilizationRatio();
				end
				%% Reconfiguration Cost
				if contains(stat_names{i},{'All', 'Cost'},'IgnoreCase',true)
					stat.Cost = this.op.get_reconfig_cost('const');
				end
				if contains(stat_names{i},{'All', 'LinearCost'},'IgnoreCase',true)
					stat.LinearCost = this.op.get_reconfig_cost('linear', true);
				end
				%% Number of Variables
				if contains(stat_names{i},{'All', 'NumberVariables'},'IgnoreCase',true)
					stat.NumberVariables = this.op.get_reconfig_state(s, 'NumberVariables');
				end
				if contains(stat_names{i},{'All', 'NumberReconfigVariables'},'IgnoreCase',true)
					stat.NumberReconfigVariables = this.op.get_reconfig_state(s, 'NumberReconfigVariables');
				end
				%% Number of Flows
				if contains(stat_names{i},{'All', 'NumberFlows'},'IgnoreCase',true)
					% For convenience of comparison, we store the number of flows including the
					% removed one.
					stat.NumberFlows = max(this.NumberFlows, height(this.old_state.flow_table));
				end
				if contains(stat_names{i},{'All', 'NumberReconfigFlows'},'IgnoreCase',true)
					stat.NumberReconfigFlows = this.op.get_reconfig_state(s, 'NumberReconfigFlows');
				end
				if contains(stat_names{i},{'All', 'ResourceCost'},'IgnoreCase',true)
					options = structmerge(options, ...
						getstructfields(this.options, 'PricingPolicy', 'default', {'quadratic'}));
					stat.ResourceCost = this.getCost(options.PricingPolicy);
				end
				if contains(stat_names{i},{'All', 'FairIndex'},'IgnoreCase',true)
					stat.FairIndex = (sum(this.FlowTable.Rate))^2/...
						(this.NumberFlows*sum(this.FlowTable.Rate.^2));
				end
				if contains(stat_names{i},{'All', 'Interval'},'IgnoreCase',true)
					if this.b_dim
						stat.Interval = this.time.DimensionIntervalModified;
					else
						stat.Interval = this.time.ConfigureInterval;
					end
				end
			end
		end
			
	end
	
	methods (Access=protected)
		function s = save_state(this)
			this.old_state.flow_table = this.FlowTable;
			this.old_state.link_load = this.Links.Load;
			this.old_state.link_capacity = this.Links.Capacity;
			this.old_state.node_load = this.ServiceNodes.Load;
			this.old_state.node_capacity = this.ServiceNodes.Capacity;
      this.old_state.number_paths = this.NumberPaths;
      this.old_state.path_owner = this.path_owner;
			
			sop = this.op.save_state();
      
      this.old_state = structmerge(this.old_state, sop);
			if nargout >= 1
				s = this.old_state;
			end
    end
    
		function restore_state(this, s)
      if nargin <=1 
        s = this.old_state;
        this.op.restore_state();
      end
			this.FlowTable = s.flow_table;
      % As the removing/adding flow operations do not change the original flow data, we
      % shallow copy the flow table.
			% As path is shallow copyed, its local identifier might have been changed when removing
			% flows. So in the restore stage, the ID is recaculated. 
			pid = 1;
			this.path_owner = zeros(size(s.path_owner));
			for fid = fidx
				for p = 1:this.FlowTable{fid, 'Paths'}.Width
					this.FlowTable{fid, 'Paths'}{p}.local_id = pid;
					this.path_owner(pid) = fid;
					pid = pid + 1;
				end
			end
			assert(this.path_owner, s.path_owner);
			this.Links.Load = s.link_load;
			this.Links.Capacity = s.link_capacity;
			this.ServiceNodes.Load = s.node_load;
			this.ServiceNodes.Capacity = s.node_capacity;
    end
		
	end
	
	methods (Access = protected)
		%% Copy
		function newobj = copyElement(this)
      if this.isShallowCopyable
        newobj = copyElement@DynamicSlice(this);
        newobj.isShallowCopyable = false;
        newobj = copyElement@EventSender(newobj);
        newobj.isShallowCopyable = true;
      else
        newobj = this;
      end
			%% Reset the listener of the new instance
			% We should reconfigure the listeners by using AddListeners
			% outside.
			% see <EventSender>, <RepeatSliceReconfiguration>.
			newobj.ClearListener();
		end
		
		function updateFlowTable(this, action, varargin)
			switch action
				case 'remove'
					%% Update flow table and local path id.
					% After removing flows, re-allocate local identifier to subsequent flows (flow
					% ID) and paths (path ID).
					flow_index = varargin{1};
					Nf = height(this.old_state.flow_table);
					pid = this.FlowTable{flow_index(1), 'Paths'}(1).local_id-1;
					fidt = flow_index(1):(Nf-length(flow_index));
					for fid = fidt
						path_list = this.FlowTable{fid, 'Paths'};
						for j = 1:path_list.Width
							pid = pid + 1;
							path_list{j}.local_id = pid;
							this.path_owner(pid) = fid;
						end
					end
			end
		end
		
		function release_resource_description(this)
			%% Release Unused Resource Description
			% release unused resource description: datacenters, nodes, links;
			%   a datacenter is unused when its capcacity is zero;
			%   a node is unused only when all adjecent links are unused;
			%   a link is unused when its capcity is zero and no paths use it;
			%
			this.save_state;
			b_removed_dcs = this.ServiceNodes.Capacity <= eps;
			b_removed_links = this.Links.Capacity <= eps;
			for i = 1:this.NumberFlows
				pathlist = this.FlowTable{i, 'Paths'};
				for j = 1:pathlist.Width
					path = pathlist{j};
					h = path.node_list(1:(end-1));
					t = path.node_list(2:end);
					eid = this.graph.IndexEdge(h,t);
					b_removed_links(eid) = false;
				end
			end
			b_removed_nodes = this.graph.Remove([], b_removed_links);
			
			%% Update variables
			% see also <OnRemoveFlow>.
			this.net_changes.NodeIndex = b_removed_nodes;
			this.net_changes.DCIndex = b_removed_dcs;
			this.net_changes.LinkIndex = b_removed_links;
			this.op.identify_change(false(this.NumberPaths,1));
			this.op.I_dc_path(b_removed_dcs, :) = [];
			this.op.I_edge_path(b_removed_links, :) = [];
			this.op.Variables.z(this.changed_index.z) = [];
			this.op.Variables.v(this.changed_index.v) = [];
			%             this.vnf_reconfig_cost(this.changed_index.v) = [];
			
			%% Update resources
			this.PhysicalNodeMap{:,'Node'} = 0;
			this.Nodes(b_removed_nodes,:) = [];
			this.PhysicalNodeMap{this.Nodes.PhysicalNode,'Node'} = ...
				(1:this.NumberNodes)';
			this.PhysicalLinkMap{:,'Link'} = 0;
			this.Links(b_removed_links,:) = [];
			this.PhysicalLinkMap{this.Links.PhysicalLink,'Link'} = ...
				(1:this.NumberLinks)';
			dc_node_index = this.ServiceNodes.Node;
			this.ServiceNodes(b_removed_dcs, :) = [];
			this.Nodes{dc_node_index(b_removed_dcs), 'ServiceNode'} = 0;
			dc_node_index(b_removed_dcs) = [];
			this.Nodes{dc_node_index, 'ServiceNode'} = 1:this.NumberServiceNodes;
		end
		
	end

end
