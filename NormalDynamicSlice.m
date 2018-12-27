classdef NormalDynamicSlice < DynamicSlice	
	%% Constructor
	methods
		function this = NormalDynamicSlice(slice_data)
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
				this.op = NormalDynamicSliceOptimizer(this);
			else
				this.op = NormalDynamicSliceOptimizer(this, options);
			end
			op = this.op;
			op.initializeState();
		end
		
		function finalize(this, prices)
			finalize@DynamicSlice(this, prices);
		end
		
		function c = getLinkCapacity(this, isfinal)
			if nargin == 1 || isfinal
				c = this.Links.Capacity;
			else
				if this.Optimizer.invoke_method == 0
					c = getLinkCapacity@NormalSlice(this, isfinal);
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
					c = getNodeCapacity@NormalSlice(this, isfinal);
				else
					% since we do not reconfigure VNF capacity through fast slice reconfiguration,
					% the sum of VNF capacity equals to the node capacity.
					c = sum(reshape(this.op.temp_vars.v, ...
						this.NumberServiceNodes,this.NumberVNFs),2);
				end
			end
		end
		
	end
	
	methods (Access = protected)
		function release_resource_description(this)
			error('error: not implemented!');
			%% Release Unused Resource Description
			% Release unused resource description: datacenters, nodes, links;
			%   A link is unused when its capcity is zero and no (candidate paths of) flow uses
			%   it; 
			%   A node is unused only when all adjecent links are unused; If a node is unused
			%   and it is co-located with a DC, the DC must also be unused.
			%   A datacenter is unused when its capcacity is zero, and the co-located node is
			%   unused. 
			this.save_state;
			b_removed_dcs = this.ServiceNodes.Capacity <= eps;
			b_removed_links = this.Links.Capacity <= eps;
			%%
			if ~isempty(this.op.I_flow_path)
				for i = 1:this.NumberFlows
					b_removed_links(this.op.I_flow_edge(i,:)) = false;
				end
			end % otherwise existing flows potentially use all edges.
			b_removed_nodes = this.graph.Remove([], b_removed_links);
			%% Update variables
			% see also <OnRemoveFlow>.
			% A node can be removed only when the co-located DC should be removed.
			if ~isempty(intersect(this.ServiceNodes{~b_removed_dcs, 'VirtualNode'},...
					find(b_removed_nodes)))
				error('error: datacenter is still alive while the colocated node is removed.');
			end
			if ~isempty(intersect(this.ServiceNodes{b_removed_dcs, 'VirtualNode'},...
					find(~b_removed_nodes)))
				warning('datacenter to be removed before the co-located node.');
			end
			b_removed_dcs = b_removed_dcs & (b_removed_nodes(this.ServiceNodes.VirtualNode));
			this.net_changes.NodeIndex = find(b_removed_nodes);
			this.net_changes.DCIndex = find(b_removed_dcs);
			this.net_changes.LinkIndex = find(b_removed_links);
			this.op.identify_change([]);
			this.op.I_edge_path(b_removed_links, :) = [];
			this.op.Variables.z(this.changed_index.z) = [];
			this.op.Variables.v(this.changed_index.v) = [];
			%             this.vnf_reconfig_cost(this.changed_index.v) = [];
			
			%% Update resources
			this.Links(b_removed_links,:) = [];
			this.Nodes(b_removed_nodes,:) = [];
			this.ServiceNodes(b_removed_dcs, :) = [];
			dc_node_index = find(this.Nodes.ServiceNode~=0);
			this.Nodes{dc_node_index, 'ServiceNode'} = 1:this.NumberServiceNodes;
			this.ServiceNodes.VirtualNode = dc_node_index;

			this.op.OnUpdateResources('release');
			warning('Debug Required!');
		end
	end
		
end

