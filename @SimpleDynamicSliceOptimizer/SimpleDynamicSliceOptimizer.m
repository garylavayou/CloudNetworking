classdef SimpleDynamicSliceOptimizer < SimpleSliceOptimizer & IDynamicSliceOptimizer
	
  %% Constructor
  methods
    function this = SimpleDynamicSliceOptimizer(slice, options)
			if nargin <= 1
				args = {};
			else
				args = {options};
			end
      this@SimpleSliceOptimizer(slice, args{:});
      this@IDynamicSliceOptimizer(slice, args{:});
    end
  end
  
  %% Public Methods
  methods
		%% Initialize Problem
		% 
		function initializeParallel(this, procedure, options)
			if nargin <= 2
				options = Dictionary();
			else
				options = Dictionary(options);
			end
			initializeParallel@SimpleSliceOptimizer(this, procedure, options);
			initializeParallel@IDynamicSliceOptimizer(this, procedure, options);
		end
		
		function setProblem(this, varargin) 
			setProblem@SimpleSliceOptimizer(this, varargin{:});
			setProblem@IDynamicSliceOptimizer(this, varargin{:});
		end
		
		% Need to be public, for subclasses override it to provider an interface to network.
		[exitflag,fidx] = executeMethod(this, action);
				
		%%
		% *Update state*: when single flow arrive or depart, update the state record,
		% including: |I_dc_path, I_edge_path, I_flow_path, path_owner| and |As|.
		% Since only add/remove entry for only one flow, this operation is more efficient
		% than <initializeState>.
		% 
		% See also <DynamicSlice>.<onAddingFlow>.
		function onAddingFlow(this, fidx)
			assert(isrow(fidx), 'fidx should be an row vector.');
			slice = this.hs;
			num_exist_paths = this.old_state.number_paths;
			num_exist_flows = height(slice.old_state.flow_table);

			changed_path_index((num_exist_paths+1):slice.NumberPaths) = true;
			this.identify_change(changed_path_index);

			if num_exist_flows == 0
				pid = 0;
			else
				paths = slice.FlowTable{num_exist_flows, 'Paths'};
				pid = paths{paths.Width}.local_id;
			end
			slice.updateFlowTable('add', fidx, struct('pid', pid));
			for fid = fidx
				for p = 1:slice.FlowTable{fid, 'Paths'}.Width
					pid = pid + 1;
					this.I_flow_path(fid, pid) = 1;
					path = slice.FlowTable{fid, 'Paths'}{p};
					for k = 1:path.Length
						e = path.Link(k);
						eid = slice.graph.IndexEdge(e(1),e(2));
						this.I_edge_path(eid, pid) = 1;
						dc_index = slice.Nodes{e(1),'ServiceNode'};
						if dc_index~=0
							this.I_dc_path(dc_index, pid) = 1;
						end
					end
					dc_index = slice.Nodes{e(2),'ServiceNode'}; % last node
					if dc_index~=0
						this.I_dc_path(dc_index, pid) = 1;
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
			Np = slice.NumberPaths;
			this.old_variables.x = zeros(Np,1);
			this.old_variables.z = zeros(slice.NumberVNFs*Np*slice.NumberServiceNodes,1);
			this.old_variables.x(~this.changed_index.x) = this.old_state.variables.x;
			this.old_variables.z(~this.changed_index.z) = this.old_state.variables.z;
			this.topts.old_variables = getstructfields(this.old_variables, {'x', 'z'});
			this.update_reconfig_costinfo('add');
			this.sh_options.action = 'add';
		end
		
		function changes = onRemovingFlow(this, flow_index)
			slice = this.hs;
			changes.path_index = false(this.old_state.number_paths,1);
			for fid = flow_index
				for p = 1:slice.FlowTable{fid, 'Paths'}.Width
					pid = slice.FlowTable{fid, 'Paths'}{p}.local_id;
					changes.path_index(pid) = true;
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
			this.identify_change(changes.path_index);
			this.old_variables = getstructfields(this.old_state.variables, {'x', 'z'});
			this.topts.old_variables = Dictionary(...
				'x', this.old_variables.x(~this.changed_index.x), ...
				'z', this.old_variables.z(~this.changed_index.z));
			this.update_reconfig_costinfo('remove');
			% Before performing optimization, there is not situation that VNF instances
			% will be removed. So there is no need to copy the |old_variable.v| like x/z.
			
			%% Update incident matrix
			% remove the path/flow-related rows/columns.
			this.I_edge_path(:, changes.path_index) = [];
			this.I_dc_path(:, changes.path_index) = [];
			this.I_flow_path(:, changes.path_index) = [];
			this.I_flow_path(flow_index, :) = [];
			
			this.sh_options.action = 'remove';
			changes.Nf = height(slice.FlowTable);
			changes.pid = slice.FlowTable{flow_index(1), 'Paths'}{1}.local_id;
			%%
			% return the changes to slice to update the flow table
		end
		
		%% TODO: only update the colums arriving/removing
		% since the resource changes after each slice configuration, the matrix needs be
		% computed each time. Therefore, we define the access function to evalue it.
		% This decoupt the computation of incident matrix and As_res, in
		% <initializeState>.
		%         function As = getAs_res(this)
		%         end
 
    function update_options(this, options)
      update_options@IDynamicSliceOptimizer(this, options);
      update_options@SimpleSliceOptimizer(this, options);
    end
    
    function s = save_state(this)
      this.old_state.vnf_capacity = this.getVNFCapacity();
      this.old_state.I_dc_path = this.I_dc_path;
      this.old_state.I_edge_path = this.I_edge_path;
      this.old_state.I_flow_path = this.I_flow_path;
      this.old_state.variables = Dictionary(this.Variables);
      this.old_state.x_reconfig_cost = this.x_reconfig_cost;
      this.old_state.z_reconfig_cost = this.z_reconfig_cost;
      this.old_state.vnf_reconfig_cost = this.vnf_reconfig_cost;
      this.old_state.number_paths = this.hs.NumberPaths;
      % this.old_state.As_res = this.As_res;
      if nargout >= 1
        s = this.old_state;
      end
    end
    
    function restore_state(this, s)
      if nargin <= 1
        s = this.old_state;
      end
      this.I_dc_path = s.I_dc_path;
      this.I_edge_path = s.I_edge_path;
      this.I_flow_path = s.I_flow_path;
      this.Variables = Dictionary(s.Variables);
      this.x_reconfig_cost = s.x_reconfig_cost;
      this.z_reconfig_cost = s.z_reconfig_cost;
      this.vnf_reconfig_cost = s.vnf_reconfig_cost;
      this.As_res = s.As_res;
    end
    
		function stat = get_reconfig_stat(this, type)
			s = this.diffstate(true);
			tol_vec = this.options.DiffNonzeroTolerance;
			slice = this.hs;
			switch type
				case 'NumVariables'
					stat = ...
						nnz(s.tI_edge_path)+nnz(s.tI_dc_path)*slice.NumberVNFs;
					if isfield(s, 'mid_v')
						stat = stat + nnz(s.mid_v);
					end
				case 'ReVariables'
					%% Number of reconfigurations between current and previous solution
					% NOTE: resource allocation and release will not take place at the same
					% time.
					%   Number of edge variables, VNF assignment variables
					%   the reconfiguration of VNF instances.
					paths_length = sum(s.tI_edge_path,1);
					stat = nnz(s.diff_z_norm>tol_vec) + dot(paths_length,s.diff_x_norm>tol_vec);
					if isfield(s, 'diff_v_norm')
						stat = stat + nnz(s.diff_v_norm>tol_vec);
					end
				case 'ReFlows'
					%% Number of Reconfigured Flows
					Np = numel(s.diff_x_norm);
					Nvnf = slice.NumberVNFs;
					Nn = numel(s.diff_z_norm)/(Nvnf*Np);
					diff_path = (s.diff_x_norm>tol_vec)' + ...
						sum(sum(reshape(full(s.diff_z_norm>tol_vec),Nn,Np,Nvnf),3),1);
					if length(slice.path_owner) <= length(slice.old_state.path_owner)
						stat = numel(unique(slice.old_state.path_owner(diff_path~=0)));
					else
						stat = numel(unique(slice.path_owner(diff_path~=0)));
					end
			end
		end
		
    [utility, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts);
    
    %% optimalFlowRate
		%			[output, profit, cost] = optimalFlowRate(this, options)
    % |new_opts|:
    % * *FixedCost*: used when slice's resource amount and resource prices are fixed,
    % resource cost is fixed. When this options is specified, the base method
    % <Slice.optimalFlowRate> returns profit without resource cost.
    % See also <fcnSocialWelfare>,<optimize>.
    function varargout = optimalFlowRate(this, options)
			if nargin <= 1
				options = Dictionary();
			else
				options = Dictionary(options);
			end
			if isempty(this.pardata)
				this.initializeParallel('optimalFlowRate', options);
			end
			pardata = this.pardata;
			if pardata.NumberFlows == 0
				[varargout{1:nargout}] = this.handle_zero_flow(options);
			else
				if isfield(options, 'CostModel') && strcmpi(options.CostModel, 'fixcost')
					theta0 = this.GLOBAL_OPTIONS.theta0;
					pardata.NodeCapacity = theta0*pardata.NodeCapacity;
					pardata.LinkCapacity = theta0*pardata.LinkCapacity;
					[profit, cost, varargout{3:nargout}] = ...
						optimalFlowRate@SimpleSliceOptimizer(this, options);
					cost = cost + this.get_reconfig_cost('const');
					profit = profit - cost;
					varargout{1} = {profit};
					if nargout >= 2
						varargout{2} = cost;
					end
				else
					[varargout{1:nargout}] = optimalFlowRate@SimpleSliceOptimizer(this, options);
				end
			end
		end
	
		%% TODO postProcessing with reconfiguration examining
		% Fist check constraint violation, then post-process the reconfiguration.
		%   flow processing requirement;
		%   link capacity constraint;
		%   VNF capacity constraint;
		% Resolve constraint violation by down-scaling, which might induce more
		% reconfiguration.
		%
		% As a result, we currently ignore the violation, only focus on
		% how to remove unacessary reconfigurations.
		function [b, vars] = postProcessing(this)
			%             % [TODO] May need post processing for VNF capacity constraint;
			%% post process reconfiguration
			postProcessing@SimpleSliceOptimizer(this);
			if isfield(this.temp_vars, 'v')
				tol_zero = this.options.NonzeroTolerance;
				var_v = this.temp_vars.v;
				var_v(var_v<tol_zero*max(var_v)) = 0;
				this.Variables.v = var_v;
			end
			%% Another simple method
			% directly recover those little changes, which might cause constraint violation
			tol_vec = this.options.DiffNonzeroTolerance;
			if ~this.b_derive_vnf
				s = this.diffstate(true);
				if length(s.diff_x) > length(this.Variables.x)
					s.diff_x_norm = s.diff_x_norm(~this.changed_index.x);
					s.diff_z_norm = s.diff_z_norm(~this.changed_index.z);
					old_x = this.old_variables.x(~this.changed_index.x);
					old_z = this.old_variables.z(~this.changed_index.z);
				else
					old_x = this.old_variables.x;
					old_z = this.old_variables.z;
				end
				b_diss_x = s.diff_x_norm < tol_vec;
				this.Variables.x(b_diss_x) = old_x(b_diss_x);
				b_diss_z =  s.diff_z_norm < tol_vec;
				this.Variables.z(b_diss_z) = old_z(b_diss_z);
				if isfield(this.temp_vars, 'v')
					if length(s.diff_v) > length(var_v)
						s.diff_v_norm = s.diff_v_norm(~this.changed_index.v);
						old_v = this.old_variables.v(~this.changed_index.v);
					else
						old_v = this.old_variables.v;
					end
					b_diss_v =  s.diff_v_norm < tol_vec;
					var_v(b_diss_v) = old_v(b_diss_v);
					this.Variables.v = var_v;
				end
				this.diff_state = struct([]);
				b = true; 
				if nargout>= 2
					vars = this.Variables.data;
				end
				return;
				%{
                find(this.As_res * [this.Variables.x; this.Variables.z]>0,1)
                res = this.I_edge_path*this.Variables.x-this.Links.Capacity;
                find(res>0)
                disp(res)
                this.Hdiag*this.Variables.z-this.Variables.v>0
				%}
			end
			if ~this.b_derive_vnf
				%% do processing to discard minor changes.
				% When performing slice dimensioning, this only be performed after the
				% final iteration. The intermediate iteration only use the linear
				% approximated cost.
				s = this.diffstate(true);
				if length(s.diff_x) > length(this.Variables.x)
					s.diff_x = s.diff_x(~this.changed_index.x);
					s.diff_x_norm = s.diff_x_norm(~this.changed_index.x);
					s.diff_z = s.diff_z(~this.changed_index.z);
					s.diff_z_norm = s.diff_z_norm(~this.changed_index.z);
					old_x = this.old_variables.x(~this.changed_index.x);
					old_z = this.old_variables.z(~this.changed_index.z);
				else
					old_x = this.old_variables.x;
					old_z = this.old_variables.z;
				end
				restore_x = this.Variables.x;
				restore_z = this.Variables.z;
				if isfield(this.temp_vars, 'v')
					restore_v = this.temp_vars.v;
					if length(s.diff_v) > length(restore_v)
						s.diff_v = s.diff_v(~this.changed_index.v);
						s.diff_v_norm = s.diff_v_norm(~this.changed_index.v);
						s.mid_v = s.mid_v(~this.changed_index.v);
						old_v = this.old_variables.v(~this.changed_index.v);
					else
						old_v = this.old_variables.v;
					end
					b_diss_v = (s.diff_v<0) & (s.diff_v_norm<tol_vec);
					restore_v(b_diss_v) = old_v(b_diss_v);
					s.diff_v_norm(b_diss_v) = 0;
					s.diff_v(b_diss_v) = 0;
				else
					restore_v = this.Variables.v;
				end
				%                 b_diss_z = (s.diff_z_norm>0) && (s.diff_z_norm<tol_vec) && s.diff_z>0;
				%                 if length(this.old_variables.x) <= length(this.temp_vars.x)
				%                     this.Variables.x(b_diss_x) = this.old_variables.x(b_diss_x);
				%                 else
				%
				%                 end
				% separately process the two part of flows, firstly process the path's
				% that might release some resource could improve the possibility to accept
				% the second part of flows, which need a little more resources.
				%%
				% first case: |x| does not change or the increase amount less than
				% |tol_vec|, try to restore these changes.
				b_diss_x1 = s.diff_x>=0;   % x might not change while z might still change
				slice = this.hs;
				Ndc = slice.NumberServiceNodes;
				Np = slice.NumberPaths;
				Nvnf = slice.NumberVNFs;
				af = slice.Parent.VNFTable{slice.VNFList, 'ProcessEfficiency'};
				%% 1-1
				for p = (find(b_diss_x1))'
					if s.diff_x_norm(p)<tol_vec
						% If x has significant change, usually z also has significant
						% change. But we cannot claim that all z components has
						% significant changes.
						restore_x(p) = old_x(p);
						s.diff_x_norm(p) = 0;
						s.diff_x(p) = 0;
					end
					check_z(1);
				end
				%% 1-2
				delta_vnf = restore_v - this.Hdiag*restore_z;  % residual VNF capacity
				for p = (find(b_diss_x1))'
					check_z(2);
				end
				%% 2-1
				% second case: |x|'s decrease amount less than |tol_vec|, try to restore
				% these changes.
				b_diss_x2 = s.diff_x<0;
				% Links.Capacity has been updated in <finalize>.
				if isfield(this.temp_vars, 'c')
					link_capacity = this.temp_vars.c;
				else
					link_capacity = slice.Links.Capacity;
				end
				delta_link = link_capacity - this.I_edge_path*restore_x; % residual link capacity
				for p = (find(b_diss_x2))'
					if s.diff_x_norm(p)<tol_vec
						temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
						if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
							temp_x = restore_x(p);
							restore_x(p) = old_x(p); % restore: up-scaling
							tf = check_z(1);
							if tf == false % up-scaling restore failed, recover current value.
								restore_x(p) = temp_x;
							else
								delta_link = temp_delta_link;
								s.diff_x_norm(p) = 0;
								s.diff_x(p) = 0;
							end
						end
					else
						check_z(1);
					end
				end
				%% 2-2
				for p = (find(b_diss_x2))'
					if s.diff_x_norm(p)<tol_vec
						temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
						if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
							temp_x = restore_x(p);
							restore_x(p) = old_x(p); % restore: up-scaling
							tf = check_z(2);
							if tf == false % up-scaling restore failed, recover current value.
								restore_x(p) = temp_x;
							else
								delta_link = temp_delta_link;
								s.diff_x_norm(p) = 0;
								s.diff_x(p) = 0;
							end
						end
					else
						check_z(2);
					end
				end
				%% 1-3
				for p = (find(b_diss_x1))'
					check_z(3);
				end
				%% 2-3
				for p = (find(b_diss_x2))'
					if s.diff_x_norm(p)<tol_vec
						temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
						if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
							temp_x = restore_x(p);
							restore_x(p) = old_x(p); % restore: up-scaling
							tf = check_z(3);
							if tf == false % up-scaling restore failed, recover current value.
								restore_x(p) = temp_x;
							else
								delta_link = temp_delta_link;
								s.diff_x_norm(p) = 0;
								s.diff_x(p) = 0;
							end
						end
					else
						check_z(3);
					end
				end
				%% reover VNF capacity
				% must be performed after the flow assignement hast been done. Otherwise,
				% there is no chance to restore VNF capacity. Due to the reconfiguration
				% cost constraint, there is no redundant reconfiguration.
				% Instead, after the flow has been restored, there might be redundancy of
				% VNF capacity.
				if isfield(this.temp_vars, 'v')
					b_diss_v = (s.diff_v>0) & (s.diff_v_norm<tol_vec);
					idv = this.Hdiag*restore_z<=old_v & b_diss_v; % if the old value satisfy the capacity constraint
					restore_v(idv) = old_v(idv);
					s.diff_v_norm(idv) = 0;
					s.diff_v(idv) = 0;
					
					%                     delta_dc = full(sum(reshape(this.Variables.v, Ndc, Nvnf),2)) -...
					%                         full(sum(reshape(restore_v, Ndc, Nvnf),2)); % residual node capacity
					%                     % all nodes can be check in parallel, with VNF varying
					%                     idv = 1:Ndc;
					%                     for v = 1:Nvnf
					%                         tv = restore_v(idv);
					%                         old_tv = old_v(idv);
					%                         idtv = (s.diff_v_norm(idv)>0) & (s.diff_v_norm(idv)<tol_vec) & s.diff_v(idv)<0;
					%                         tv(idtv) = old_tv(idtv);
					%                         temp_delta_dc = delta_dc + restore_v(idv) - tv;
					%                         fx = temp_delta_dc >= 0;
					%                         restore_v(idv(fx)) = tv(idv(fx));
					%                         delta_dc(fx) = temp_delta_dc(fx);
					%                         s.diff_v_norm(idv(fx)) = 0;
					%                         s.diff_v(idv(fx)) = 0;
					%                         idv = idv + Ndc;
					%                     end
					this.Variables.v = restore_v;
				end
				this.Variables.x = restore_x;
				this.Variables.z = restore_z;
				this.diff_state = struct([]);
			end
			b = true; 
			if nargout >= 2
				vars = this.Variables.data;
			end
			this.setPathBandwidth();
			%% nest function: check z
			function tf = check_z(t)
				tf = true;
				b_reset = false(Nvnf,1);
				temp_z = zeros(Ndc, Nvnf);
				temp_sz = zeros(Ndc, Nvnf);
				temp_szn = zeros(Ndc, Nvnf);
				delta_sum = zeros(Nvnf,1);
				if t > 1
					temp_delta_vnf = zeros(Ndc, Nvnf);
				end
				
				idz = (1:Ndc) + (p-1)*Ndc;
				idn = 1:Ndc;
				for k = 1:Nvnf
					tz = restore_z(idz);
					temp_z(:,k) = tz;
					temp_sz(:,k) = s.diff_z(idz);
					temp_szn(:,k) = s.diff_z_norm(idz);
					if t > 1
						temp_delta_vnf(:,k) = delta_vnf(idn);
					end
					old_tz = old_z(idz);
					switch t
						case 1
							% try restore z(:,p,f) that should be restored by down-scaling;
							% make space for other requests.
							idtz = (s.diff_z(idz)>0) & (s.diff_z_norm(idz)<tol_vec);
						case 2
							%%% try restore z(:,p,f) that includes both up-sclaing and down-scaling;
							idtz = (s.diff_z_norm(idz)>0) & (s.diff_z_norm(idz)<tol_vec);
						case 3
							%%% try restore z(:,p,f) that only includes up-scaling;
							idtz = (s.diff_z(idz)<0) & (s.diff_z_norm(idz)<tol_vec);
					end
					if isempty(find(idtz,1))
						if nargout == 1 && dot(this.I_dc_path(:,p), tz) < af(k)*restore_x(p)
							% when x increase, but z keeps unchange, so we need to check the processing
							% constraint.
							tf = false;
						end
						idz = idz + Np*Ndc;
						idn = idn + Ndc;
						continue;
					end
					tz(idtz) = old_tz(idtz);  % recover: down-scaling | down/up-scaling | up-scaling
					if t == 2 || t == 3
						delta_tz = restore_z(idz) - tz; % current - past
						t_delta_vnf = delta_vnf(idn)+this.I_dc_path(:,p).*delta_tz;
					end
					if t == 1
						if dot(this.I_dc_path(:,p), tz) >= af(k)*restore_x(p) % check processing constraint
							restore_z(idz) = tz;
							s.diff_z_norm(idz(idtz)) = 0;
							s.diff_z(idz(idtz)) = 0;
							delta_sum(k) = 0;
						else
							b_reset(k) = true;
						end
					end
					if t == 2
						if dot(this.I_dc_path(:,p), tz) >= af(k)*restore_x(p) &&...
								isempty(find(t_delta_vnf<0,1))
							%% check both processing constraint and VNF capacity constraint
							restore_z(idz) = tz; % accept old value or keep current value.
							s.diff_z_norm(idz(idtz)) = 0;
							s.diff_z(idz(idtz)) = 0;
							delta_vnf(idn) = t_delta_vnf; % update residual capacity
							delta_sum(k) = sum(delta_tz);
						else
							b_reset(k) = true;
						end
					end
					if t == 3
						if (nargout == 0 && isempty(find(t_delta_vnf<0,1))) || (nargout == 1 && isempty(find(t_delta_vnf<0,1)) && dot(this.I_dc_path(:,p), tz) >= af(k)*restore_x(p))
							restore_z(idz) = tz; % accept old value or keep current value.
							s.diff_z_norm(idz(idtz)) = 0;
							s.diff_z(idz(idtz)) = 0;
							delta_vnf(idn) = t_delta_vnf; % update residual capacity
							delta_sum(k) = sum(delta_tz);
						else
							b_reset(k) = true;
						end
					end
					idz = idz + Np*Ndc;
					idn = idn + Ndc;
				end
				if ~isempty(find(b_reset,1))
					tf = false;
					idz = (1:Ndc) + (p-1)*Ndc;
					idn = 1:Ndc;
					for k = 1:Nvnf
						if ~b_reset(k) && delta_sum(k) < 0
							restore_z(idz) = temp_z(:,k);
							s.diff_z(idz) = temp_sz(:,k);
							s.diff_z_norm(idz) = temp_szn(:,k);
							if t > 1
								delta_vnf(idn) = temp_delta_vnf(:,k);
							end
						end
						idz = idz + Np*Ndc;
						idn = idn + Ndc;
					end
				end
			end
    end
  end
  
  methods (Access=protected)
    [profit,cost] = fastReconfigure2(this, action, options);
    [profit,cost] = fastReconfigure(this, action, options);
		% |parameters|: include fields x0, As, bs, Aeq, beq, lbs, ub, lb;
		% |options|: include fields minopts, CostModel, PricingPolicy (if
		%       CostModel='fixcost'), num_orig_vars;
		function [x, fval] = optimize(this, options)
			if isfield(options, 'CostModel') && strcmpi(options.CostModel, 'fixcost')
				prbm = this.problem;
				[xs, fval, exitflag, output] = ...
					fmincon(@(x)SimpleDynamicSliceOptimizer.fcnSocialWelfare(x, this, options), ...
					prbm.x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, prbm.ub, [], ...
					prbm.minopts);
				this.interpretExitflag(exitflag, output.message);
				if options.bCompact
					x = zeros(prbm.num_orig_vars, 1);
					x(this.I_active_variables) = xs;
				else
					x = xs;
				end
				constr_tol = getstructfields(prbm.options, {'ConstraintTolerance'}, ...
					'default-ignore', this.options.ConstraintTolerance);
				assert(this.checkFeasible(x, constr_tol), 'error: infeasible solution.');
			else
				[x, fval] = optimize@SimpleSliceOptimizer(this, rmstructfields(options, 'CostModel'));
			end
		end
 		
		%% Convert the optimizatio results to temporary variables
    % temp variables might include: x,z,v,tx,tz,tv,c (see also
    % <Dynamic.priceOptimalFlowRate>);
    function convert(this, x, ~)
      Np = this.hs.NumberPaths;
      this.temp_vars.x = x(1:Np);
			num_varxz = sum(this.num_vars(1:2));
      this.temp_vars.z = x((Np+1):num_varxz);
      offset = num_varxz;
			num_varv = this.num_vars(3);
      this.temp_vars.v = x(offset+(1:num_varv));
      offset = offset + num_varv;
      if this.invoke_method >= 2
        this.temp_vars.tx = x(offset+(1:Np));
        this.temp_vars.tz = x(offset+((Np+1):num_varxz));
        offset = offset + num_varxz;
        this.temp_vars.tv = x(offset+(1:num_varv));
        offset = offset + num_varv;
      end
      this.temp_vars.c = x(offset+(1:this.hs.NumberLinks));
      %             if nargin >= 3
      %                 this.Variables = getstructfields(this.temp_vars, {'x','z','v','c'}, 'ignore');
      %             end
		end
    
		%% TODO
		% optimize only on the subgraph, which carries the flows that intersection with
		% the arriving/departuring flow. this can significantly reduce the problem scale
		% when there is lots of flow, and the numbe of hops of the flow is small.
		% 		function temp_vars = get_temp_variables(this, bfull)
		% 			if nargin >= 2 && bfull
		% 				temp_vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v; ...
		% 					this.temp_vars.tx; this.temp_vars.tz; this.temp_vars.tv];
		% 			else
		% 				temp_vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v];
		% 			end
		% 		end
		
    %% Find the index of variables for the newly added/removed flow
    % assume the newly added/removed flow index is u, the the index of
    % corrsponding flow and VNF allocation variables is
    % [x_u, z_npf], for all p in p_u.
    %
    % NOTE: This function should be called after adding new flow, or before removing
    % departing flow.
    function identify_change(this, changed_path_index)
      global DEBUG; %#ok<NUSED>
			if ~islogical(changed_path_index)
				pidx = changed_path_index;
				changed_path_index = logical(sparse(1, this.hs.NumberPaths));
				changed_path_index(pidx) = true;
			end
      Nsn = this.hs.NumberServiceNodes;
      Nvnf = this.hs.NumberVNFs;
      this.changed_index.x = logical(sparse(changed_path_index));
      base_z_index = logical(sparse(Nsn, this.hs.NumberPaths));
      base_z_index(:,changed_path_index) = true;
      if ~isempty(fieldnames(this.hs.net_changes))
        base_z_index(this.hs.net_changes.DCIndex,:) = true;
      end
      base_z_index = sparse(base_z_index);
      this.changed_index.z = repmat(base_z_index(:), Nvnf, 1);
      if ~isempty(fieldnames(this.hs.net_changes))
        this.changed_index.v = sparse(Nsn, Nvnf);
        this.changed_index.v(this.hs.net_changes.DCIndex,:) = true;
        this.changed_index.v = logical(this.changed_index.v(:));
      end
		end
    
    % function update_reconfig_costvinfo(this)  SEE <IDynamicSliceOptimizer>.
    
		%% calculate the difference of slice state
		function ds = diffstate(this, isfinal)
			ds = diffstate@IDynamicSliceOptimizer(this, isfinal);
			
			% At the beginging, |old_variables| equals to |Variables|.
			% |old_variables| is still valid after adding/removing flows.
			% |old_variables| have more elements than |this.Variables.x| when removing
			% flows.
			if length(this.old_variables.x) <= length(this.temp_vars.x)  % only ==
				ds.tI_dc_path = this.I_dc_path;
				% adding flow, |path| will increase, and |edge| might(might not) increase;
				ds.tI_edge_path = this.I_edge_path;
			else
				ds.tI_dc_path = this.old_state.I_dc_path;
				ds.tI_edge_path = this.old_state.I_edge_path;
			end
			if isfinal
				this.diff_state = ds;
			end
		end

	end  % end protected methods
	
	methods (Access = {?IDynamicSliceOptimizer, ?DynamicNetwork})
		% update reconfiguration cost with scaler.
		function update_reconfig_costinfo(this, action, bDim, ~)
			if isempty(action)
				action = this.sh_options.action;
			end
			%% Period re-dimensioing
			% The number of virtual nodes/links is not change in the slice, as well as the
			% number of VNF instances.
			% When performing re-dimensioning, the reconfiguration cost is larger than
			% that of fast reconfiguration, so we need to update the reconfiguration cost.
			slice = this.hs;
			if nargin >= 3 && bDim
				slice.Parent.updateRedimensionCost(this.hs);
				this.update_reconfig_costvinfo();
			end
			if strcmpi(action, 'add')
				this.x_reconfig_cost = (this.I_edge_path)' * slice.Links.ReconfigCost;
				this.z_reconfig_cost = repmat(slice.ServiceNodes.ReconfigCost, ...
					slice.NumberPaths*slice.NumberVNFs, 1);
				if nargin >= 4   % the 4th argument passed in from <init_reconfig_info>
					return;
				end
				%% Ignore the reconfiguration cost when adding flow
				% The newly added flow's reconfiguration cost set to zero.
				% This helps the new flow be admitted into the slice, when there is
				% resource while the reconfiguration cost of the new flow is too high.
				% This can be explained as: the reconfiguration cost of the first time can
				% be treated as that it is distrbuted in the whole life time of the flow.
				% As long as the flow's lifetime is long enough, the first time cost could
				% be ignored.
				%
				% In the output, we still count the reconfiguration cost for the new flow,
				% see also <diffstate> and <get_reconfig_cost>.
				this.x_reconfig_cost(this.changed_index.x) = 0;
				this.z_reconfig_cost(this.changed_index.z) = 0;
				if this.ENABLE_DYNAMIC_NORMALIZER
					beta =this.getbeta();
					this.topts.x_reconfig_cost = beta.x*this.x_reconfig_cost;
					this.topts.z_reconfig_cost = beta.z*this.z_reconfig_cost;
				else
					this.topts.x_reconfig_cost = this.options.ReconfigScaler*this.x_reconfig_cost;
					this.topts.z_reconfig_cost = this.options.ReconfigScaler*this.z_reconfig_cost;
				end
			else
				% We might call this method after the flow has been removed (see
				% <priceOptimalFlowRate>). Therefore, we should use the record of flow
				% information.
				this.x_reconfig_cost = (this.old_state.I_edge_path)' * slice.Links.ReconfigCost;
				old_num_paths = this.old_state.number_paths;
				this.z_reconfig_cost = ...
					repmat(slice.ServiceNodes.ReconfigCost, old_num_paths*slice.NumberVNFs, 1);
				%%
				% The removed flow's reconfiguration cost is not considered in the
				% optimization, which is constant.
				if this.ENABLE_DYNAMIC_NORMALIZER
					beta =this.getbeta();
					this.topts.x_reconfig_cost = beta.x*this.x_reconfig_cost(~this.changed_index.x);
					this.topts.z_reconfig_cost = beta.z*this.z_reconfig_cost(~this.changed_index.z);
				else
					this.topts.x_reconfig_cost = ...
						this.options.ReconfigScaler*this.x_reconfig_cost(~this.changed_index.x);
					this.topts.z_reconfig_cost = ...
						this.options.ReconfigScaler*this.z_reconfig_cost(~this.changed_index.z);
				end
			end
		end
		
	end
	
	methods (Static, Access = protected)
		[profit, grad] = fcnProfitReconfigureSlicing(vars, slice, options);
		[profit, grad] = fcnProfitReserveSlicing(vars, slice, options);
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
			numvar = this.problem.num_vars;
			num_basic_vars = sum(numvar(1:2));
			% |tx|,|tz| and |tv| are auxilliary variables to transform L1 norm to
			% linear constraints and objective part. The number of auxiliary variables
			% eqauls to the number of the main variables.
			num_var_tx = length(this.topts.x_reconfig_cost);
			num_var_tz = length(this.topts.z_reconfig_cost);
			if length(numvar) >= 3  % for _fastReconfigure2_, including |v|, and |tv|
				var_v_index = num_basic_vars + (1:num_vars(3));
				num_basic_vars = num_basic_vars + num_vars(3);
				num_var_tv = length(this.topts.vnf_reconfig_cost);
				var_tv_index = num_basic_vars + num_var_tx + num_var_tz + (1:num_var_tv);
				var_tv = vars(var_tv_index);
			end
			var_x = vars(1:numvar(1));
			var_tx_index = num_basic_vars+(1:num_var_tx);
			var_tx = vars(var_tx_index);
			var_tz_index = num_basic_vars+num_var_tx+(1:num_var_tz);
			var_tz = vars(var_tz_index);
			flow_rate = this.getFlowRate(var_x);
			%% objective value
			slice = this.hs;
			profit = -slice.Weight*sum(fcnUtility(flow_rate));
			%%%
			% calculate reconfiguration cost by |tx| and |tz| (L1 approximation),
			% which is similar to calculating by |x-x0| and |z-z0|.
			profit = profit + dot(var_tx, this.topts.x_reconfig_cost) + ...
				dot(var_tz, this.topts.z_reconfig_cost);
			if length(numvar) >= 3 % for _fastConfigure2_
				profit = profit + dot(var_tv, this.topts.vnf_reconfig_cost);
			end
			% If there is only one output argument, return the real profit (positive)
			if options.bFinal
				profit = -profit;
			else
				%% Gradient
				% The partial derivatives are computed by dividing the variables into four
				% parts, i.e., $x,z,t_x,t_z$. Since |z| does not appear in objective
				% function, the corrsponding derivatives is zero.
				% The partial derivatives of x
				grad = spalloc(num_vars, 1, num_vars-numvar(2));
				for p = 1:numvar(1)
					i = slice.path_owner(p);
					grad(p) = -slice.Weight/(1+this.I_flow_path(i,:)*var_x); %#ok<SPRIX>
				end
				%%%
				% The partial derivatives of t_x is the vector |x_reconfig_cost|;
				% The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
				% since the node reconfiguration cost is only depend on node.
				grad(var_tx_index) = this.topts.x_reconfig_cost;
				grad(var_tz_index) = this.topts.z_reconfig_cost;
				if length(numvar) >= 3
					grad(var_v_index) = this.topts.vnf_reconfig_cost;
				end
			end
		end
		
		%%
		% See also <SliceOptimzier.fcnHessian>.
		function h = hessReconfigure(vars, lambda, this) %#ok<INUSL>
			num_varx = this.problem.num_vars(1);
			var_x = vars(1:num_varx);
			num_vars = length(vars);
			h = spalloc(num_vars, num_vars, num_varx^2);   % non-zero elements less than | options.num_varx^2|
			slice = this.hs;
			for p = 1:num_varx
				i = slice.path_owner(p);
				h(p,1:num_varx) = slice.Weight *...
					this.I_flow_path(i,:)/(1+(this.I_flow_path(i,:)*var_x))^2; %#ok<SPRIX>
			end
		end
	end
end