classdef SimpleDynamicSlice < DynamicSlice
	%DynamicSlice Event-driven to dynamic configure slice.
	
	% 	properties(Dependent)
	% 		UnitReconfigureCost;
	% 	end
	
	%% Constructor
	methods
		function this = SimpleDynamicSlice(slice_data)
      if nargin == 0
        args = {};
			else
				args = {slice_data};
      end
			this@DynamicSlice(args{:});
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
			op.initializeState();
		end
		
		function finalize(this, prices)
			finalize@DynamicSlice(this, prices);
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

		%% Get link and node capacity
		% <getLinkCapacity>
		% <getNodeCapacity>
		% Due to resource reservation or reconfiguration cost constraint, the link/node
		% capcity might be larger than the demand.
		function c = getLinkCapacity(this, isfinal)
			if nargin == 1 || isfinal
				c = this.Links.Capacity;
			else
				if this.op.invoke_method == 0
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
				if this.Optimizer.invoke_method == 0
					c = getNodeCapacity@SimpleSlice(this, isfinal);
				else
					% since we do not reconfigure VNF capacity through fast slice reconfiguration,
					% the sum of VNF capacity equals to the node capacity.
					c = sum(reshape(this.op.temp_vars.v, ...
						this.NumberServiceNodes,this.NumberVNFs),2);
				end
			end
		end
		
		%%%
		% At present this method inherit the superclass. See also <Slice.checkFeasible>.
		%         function b = checkFeasible(this, vars, options)
		%             if nargin <= 1
		%                 vars = [];
		%             else
		%                 vars = vars(1:num_vars);
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
	
	methods (Access = protected)
		function release_resource_description(this)
			error('error: not implemented!');
			%% Release Unused Resource Description
			% release unused resource description: datacenters, nodes, links;
			%   A link is unused when its capcity is zero and no candidate paths use it;
			%   A node is unused only when all adjecent links are unused; If a node is unused
			%   and it is colocated with a DC, the DC must also be unused.
			%   A datacenter is unused when its capcacity is zero, and the co-located node is
			%   unused. 
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
			% 'b_removed_nodes', 'b_removed_links', 'b_removed_dcs' helps to identify the
			% variables that should be removed.
			b_removed_dcs = b_removed_dcs & (b_removed_nodes(this.ServiceNodes.VirtualNode));
			this.net_changes.NodeIndex = find(b_removed_nodes);
			this.net_changes.DCIndex = find(b_removed_dcs);
			this.net_changes.LinkIndex = find(b_removed_links);
			this.op.identify_change(false(this.NumberPaths,1));
			this.op.I_dc_path(b_removed_dcs, :) = [];
			this.op.I_edge_path(b_removed_links, :) = [];
			this.op.Variables.z(this.changed_index.z) = [];
			this.op.Variables.v(this.changed_index.v) = [];
			%             this.vnf_reconfig_cost(this.changed_index.v) = [];
			
			%% Update resources
			% Should not change the relative indexing order.
			this.Links(b_removed_links,:) = [];
			this.Nodes(b_removed_nodes,:) = [];
			this.ServiceNodes(b_removed_dcs, :) = [];
			dc_node_index = find(this.Nodes.SeviceNode~=0);
			this.Nodes{dc_node_index, 'ServiceNode'} = 1:this.NumberServiceNodes;
			this.ServiceNodes.VirtualNode = dc_node_index;
			% TO BE REMOVED: this.PhysicalNodeMap{:,'Node'} = 0;
			% TO BE REMOVED: this.PhysicalNodeMap{this.Nodes.PhysicalNode,'Node'} = (1:this.NumberNodes)';
			% TO BE REMOVED: this.PhysicalLinkMap{:,'Link'} = 0;
			% TO BE REMOVED: this.PhysicalLinkMap{this.Links.PhysicalLink,'Link'} = (1:this.NumberLinks)';
			warning('Debug Required!');
		end
		
	end
end
