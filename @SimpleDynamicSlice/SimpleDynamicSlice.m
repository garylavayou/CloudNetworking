classdef SimpleDynamicSlice < DynamicSlice
	%DynamicSlice Event-driven to dynamic configure slice.
	
	properties(Dependent)
		UnitReconfigureCost;
  end
	
	methods
		function this = SimpleDynamicSlice(slice_data)
      if nargin == 0
        return;
      end
			this@DynamicSlice(slice_data);
      this.op = SimpleDynamicSliceOptimizer(this);
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
					c = this.temp_vars.c;
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
					stat.Cost = this.get_reconfig_cost('const');
				end
				if contains(stat_names{i},{'All', 'LinearCost'},'IgnoreCase',true)
					stat.LinearCost = this.get_reconfig_cost('linear', true);
				end
				%% Number of Variables
				if contains(stat_names{i},{'All', 'NumberVariables'},'IgnoreCase',true)
					stat.NumberVariables = ...
						nnz(s.tI_edge_path)+nnz(s.tI_dc_path)*this.NumberVNFs;
					if isfield(s, 'mid_v')
						stat.NumberVariables = stat.NumberVariables + nnz(s.mid_v);
					end
				end
				%% Number of reconfigurations between current and previous solution
				% NOTE: resource allocation and release will not take place at the same
				% time.
				%   Number of edge variables, VNF assignment variables
				%   the reconfiguration of VNF instances.
				if contains(stat_names{i},{'All', 'NumberReconfigVariables'},'IgnoreCase',true)
					paths_length = sum(s.tI_edge_path,1);
					stat.NumberReconfigVariables = nnz(s.diff_z_norm>tol_vec) ...
						+ dot(paths_length,s.diff_x_norm>tol_vec);
					if isfield(s, 'diff_v_norm')
						stat.NumberReconfigVariables = stat.NumberReconfigVariables + ...
							nnz(s.diff_v_norm>tol_vec);
					end
				end
				%% Number of Flows
				if contains(stat_names{i},{'All', 'NumberFlows'},'IgnoreCase',true)
					% For convenience of comparison, we store the number of flows including the
					% removed one.
					stat.NumberFlows = max(this.NumberFlows, height(this.old_state.flow_table));
				end
				%% Number of Reconfigured Flows
				if contains(stat_names{i},{'All', 'NumberReconfigFlows'},'IgnoreCase',true)
					np = numel(s.diff_x_norm);
					nv = this.NumberVNFs;
					nn = numel(s.diff_z_norm)/(nv*np);
					diff_path = (s.diff_x_norm>tol_vec)' + ...
						sum(sum(reshape(full(s.diff_z_norm>tol_vec),nn,np,nv),3),1);
					if length(this.op.path_owner) <= length(this.old_state.path_owner)
						stat.NumberReconfigFlows = ...
							numel(unique(this.old_state.path_owner(diff_path~=0)));
					else
						stat.NumberReconfigFlows = ...
							numel(unique(this.op.path_owner(diff_path~=0)));
					end
				end
				if contains(stat_names{i},{'All', 'ResourceCost'},'IgnoreCase',true)
					options = structmerge(options, ...
						getstructfields(this.options, 'PricingPolicy', 'default', 'quadratic'));
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
      for fid = fidx
        for p = 1:this.FlowTable{fid, 'Paths'}.Width
          this.FlowTable{fid, 'Paths'}.paths{p}.local_id = pid;
          pid = pid + 1;
        end
      end
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
		
		%% fast slice reconfiguration when flow arriving and depaturing
		% * *flow*: flow table entries.
		%
		% *Update state*: when single flow arrive or depart, update the state record,
		% including: |I_dc_path, I_edge_path, I_flow_path, path_owner| and |As|.
		% Since only add/remove entry for only one flow, this operation is more efficient
		% than <initializeState>.
		function fidx = OnAddingFlow(this, flows)
			% temporarily allocated flow id.
			global DEBUG total_iter_num; %#ok<NUSED>
			num_exist_flows = height(this.FlowTable);
			num_new_flows = height(flows);
			num_exist_paths = this.NumberPaths;
			fidx = (1:num_new_flows) + num_exist_flows;
			if num_exist_flows == 0
				flows{:,'Identifier'} = 1:num_new_flows;
			else
				flows{:,'Identifier'} = (1:num_new_flows) + this.FlowTable{end,'Identifier'}; % temporary identifier
			end
			if num_exist_flows == 0
				pid = 0;
			else
				pid = this.FlowTable{end, 'Paths'}.paths{end}.local_id;
			end
			
			this.save_state; % if fail to adding flow, recover flow table
			
			this.FlowTable = [this.FlowTable; flows];
			changed_path_index((num_exist_paths+1):this.NumberPaths) = true;
			this.op.identify_change(changed_path_index);
			
			for fid = fidx
				for p = 1:this.FlowTable{fid, 'Paths'}.Width
					pid = pid + 1;
					this.op.I_flow_path(fid, pid) = 1;
					this.op.path_owner(pid) = fid;
					path = this.FlowTable{fid, 'Paths'}.paths{p};
					path.local_id = pid;
					for k = 1:(path.Length-1)
						e = path.Link(k);
						eid = this.graph.IndexEdge(e(1),e(2));
						this.op.I_edge_path(eid, pid) = 1;
						dc_index = this.Nodes{e(1),'ServiceNode'};
						if dc_index~=0
							this.op.I_dc_path(dc_index, pid) = 1;
						end
					end
					dc_index = this.Nodes{e(2),'ServiceNode'}; % last node
					if dc_index~=0
						this.op.I_dc_path(dc_index, pid) = 1;
					end
				end
			end
			
			%% Previous State
			% In _onAddingFlow_ and _OnRemovingFlow_, the slice topology does not change.
			% The data centers and NVF type do not changes, so there is no new/deleted VNF
			% instance variables.
			%
			% Record the state changes after flow is added to the flow table to maintain
			% the full set of variables, different from that in _<OnRemovingFlow>_.
			this.op.old_variables.x = zeros(this.NumberPaths,1);
			this.op.old_variables.z = zeros(this.num_varz,1);
			this.op.old_variables.x(~this.changed_index.x) = this.op.old_state.variables.x;
			this.op.old_variables.z(~this.changed_index.z) = this.op.old_state.variables.z;
			this.op.topts.old_variables_x = this.op.old_variables.x;
			this.op.topts.old_variables_z = this.op.old_variables.z;
			this.op.update_reconfig_costinfo('add');
			this.op.sh_options.action = 'add';
			[ef, ~] = this.op.executeMethod('add');
			%% Failure handling
			% Zero-rate: now we leave it unhandled, the zero-rate flow will stay in
			% the slice. It may be allocated with resource at later stage. On the
			% other hand, we can remove the zero-rate flows from the slice (reject).
			if ef < 0
				fidx = [];
				this.restore_state();
			end
		end
		
		function ef = OnRemovingFlow(this, fid)
			global DEBUG; %#ok<NUSED>
			fidx = this.FlowTable.Identifier == fid;
			if isempty(find(fidx==true,1))
				warning('flow (identifier %d) not in slice (identifer %d).',...
					fid, this.Identifier);
			end
			this.save_state; % if fail to adding flow, recover slice data
			changed_path_index = false(this.NumberPaths,1);
			flow_index = find(fidx);
			for fid = flow_index
				for p = 1:this.FlowTable{fid, 'Paths'}.Width
					pid = this.FlowTable{fid, 'Paths'}.paths{p}.local_id;
					changed_path_index(pid) = true;
				end
			end
			%% Previous State
			% |old_variables|(to be renamed) is used for all methods to calculate the
			% reconfiguration cost.
			%
			% When removing flows, the deleted variables' value is treated as zeros, its
			% difference is constant (0-x0). Reconfiguration cost corresponding to the
			% previous state, should be updated before performing optimization.
			% *Link reconfigure cost* is a vector, and we know edge-path incident matrix,
			% so we can calculate the |x_reconfig_cost| for all paths.
			%
			% Record changes before flows are removed from flow table, _identify_change_
			% is called.
			this.op.identify_change(changed_path_index);
			this.op.old_variables.x = this.op.old_state.variables.x;
			this.op.old_variables.z = this.op.old_state.variables.z;
			this.op.topts.old_variables_x = this.op.old_variables.x(~this.changed_index.x);
			this.op.topts.old_variables_z = this.op.old_variables.z(~this.changed_index.z);
			this.update_reconfig_costinfo('remove');
			% Before performing optimization, there is not situation that VNF instances
			% will be removed. So there is no need to copy the |old_variable.v| like x/z.
			
			%% Update flow table and local path id.
			% After removing flows, re-allocate local identifier to subsequent flows (flow
			% ID) and paths (path ID).
			pid = this.FlowTable{flow_index(1), 'Paths'}.paths{1}.local_id-1;
			fidt = flow_index(1):(this.NumberFlows-length(flow_index));
			this.FlowTable(fidx,:) = [];
			this.op.path_owner(changed_path_index) = [];
			for fid = fidt
				path_list = this.FlowTable{fid,{'Paths'}};
				for j = 1:path_list.Width
					pid = pid + 1;
					path_list.paths{j}.local_id = pid;
					this.op.path_owner(pid) = fid;
				end
			end
			%% Update incident matrix
			% remove the path/flow-related rows/columns.
			this.op.I_edge_path(:, changed_path_index) = [];
			this.op.I_dc_path(:, changed_path_index) = [];
			this.op.I_flow_path(:, changed_path_index) = [];
			this.op.I_flow_path(fidx, :) = [];
			
			this.op.sh_options.action = 'remove';
			[ef,~] = this.op.executeMethod('remove');
			if ef < 0
				this.restore_state();
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
					path = pathlist.paths{j};
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
	
	methods(Static, Access = protected)
		[profit, grad]= fcnProfitReconfigureSlicing(vars, slice, options);
		[profit, grad]= fcnProfitReserveSlicing(vars, slice, options);
		[profit, grad] = fcnSocialWelfare(vars, slice, options);
		hs = hessSlicing(var_x, lambda, slice, options);
		hs = hessInitialSlicing(vars, lambda, slice, options);
		%%
		% The resource consumption cost is constant, since the slice will not
		% release/request resouces to/from the substrate network, and the resource prices
		% during the reconfiguration does not change. Therefore, the resource consumption
		% cost it not counted in the objective function. The objective function contains
		% the user utility and the reconfiguration cost.
		%
		% options:
		%   # |num_varx|,|num_varz|,|num_varv|: the number of main variables. For compact
		%     mode (bCompact=true), |num_varz| is less than the original number of variables |z|.
		%   # |bFinal|: set to 'true' to output the original objective function value.
		%
		% NOTE: hession matrix of the objective is only has non-zero values for x
		% (super-linear). See also <fcnHessian>.
		function [profit, grad] = fcnFastConfigProfit(vars, this, options)
			num_vars = length(vars);
			num_basic_vars = options.num_varx + options.num_varz;
			% |tx|,|tz| and |tv| are auxilliary variables to transform L1 norm to
			% linear constraints and objective part. The number of auxiliary variables
			% eqauls to the number of the main variables.
			num_var_tx = length(this.topts.x_reconfig_cost);
			num_var_tz = length(this.topts.z_reconfig_cost);
			if isfield(options, 'num_varv')  % for _fastReconfigure2_, including |v|, and |tv|
				var_v_index = num_basic_vars + (1:options.num_varv);
				num_basic_vars = num_basic_vars + options.num_varv;
				num_var_tv = length(this.topts.vnf_reconfig_cost);
				var_tv_index = num_basic_vars + num_var_tx + num_var_tz + (1:num_var_tv);
				var_tv = vars(var_tv_index);
			end
			var_x = vars(1:options.num_varx);
			var_tx_index = num_basic_vars+(1:num_var_tx);
			var_tx = vars(var_tx_index);
			var_tz_index = num_basic_vars+num_var_tx+(1:num_var_tz);
			var_tz = vars(var_tz_index);
			flow_rate = this.getFlowRate(var_x);
			%% objective value
			profit = -this.weight*sum(fcnUtility(flow_rate));
			%%%
			% calculate reconfiguration cost by |tx| and |tz| (L1 approximation),
			% which is similar to calculating by |x-x0| and |z-z0|.
			profit = profit + dot(var_tx, this.topts.x_reconfig_cost) + ...
				dot(var_tz, this.topts.z_reconfig_cost);
			if isfield(options, 'num_varv') % for _fastConfigure2_
				profit = profit + dot(var_tv, this.topts.vnf_reconfig_cost);
			end
			% If there is only one output argument, return the real profit (positive)
			if isfield(options, 'bFinal') && options.bFinal
				profit = -profit;
			else
				%% Gradient
				% The partial derivatives are computed by dividing the variables into four
				% parts, i.e., $x,z,t_x,t_z$. Since |z| does not appear in objective
				% function, the corrsponding derivatives is zero.
				% The partial derivatives of x
				grad = spalloc(num_vars, 1, num_vars-options.num_varz);
				for p = 1:options.num_varx
					i = this.path_owner(p);
					grad(p) = -this.weight/(1+this.I_flow_path(i,:)*var_x); %#ok<SPRIX>
				end
				%%%
				% The partial derivatives of t_x is the vector |x_reconfig_cost|;
				% The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
				% since the node reconfiguration cost is only depend on node.
				grad(var_tx_index) = this.topts.x_reconfig_cost;
				grad(var_tz_index) = this.topts.z_reconfig_cost;
				if isfield(options, 'num_varv')
					grad(var_v_index) = this.topts.vnf_reconfig_cost;
				end
			end
		end
		
		%%
		% See also <Slice.fcnHessian>.
		function h = hessReconfigure(vars, lambda, this, options) %#ok<INUSL>
			var_x = vars(1:options.num_varx);
			num_vars = length(vars);
			h = spalloc(num_vars, num_vars, options.num_varx^2);   % non-zero elements less than | options.num_varx^2|
			for p = 1:options.num_varx
				i = this.path_owner(p);
				h(p,1:options.num_varx) = this.weight *...
					this.I_flow_path(i,:)/(1+(this.I_flow_path(i,:)*var_x))^2; %#ok<SPRIX>
			end
		end
	end
end
