classdef NormalDynamicSliceOptimizer < NormalSliceOptimizer & IDynamicSliceOptimizer

	methods
		function this = NormalDynamicSliceOptimizer(slice, options)
			if nargin <= 1
				args = {};
			else
				args = {options};
			end
			this@NormalSliceOptimizer(slice, args{:});
			this@IDynamicSliceOptimizer(slice, args{:});
		end
	end
	
	methods
		[utility, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts);
		
		function initializeParallel(this, procedure, options)
			initializeParallel@NormalSliceOptimizer(this, procedure, options);
			initializeParallel@IDynamicSliceOptimizer(this, procedure, options);
		end
		
		function setProblem(this, varargin) 
			setProblem@NormalSliceOptimizer(this, varargin{:});
			setProblem@IDynamicSliceOptimizer(this, varargin{:});
		end
		
		%% Can move to IDynamicSliceOptimizer
		% Need to be public, for subclasses override it to provider an interface to network.
		[exitflag,fidx] = executeMethod(this, action);
		
		%% OnAddlingFlow
		% When single flow arrive, incrementally update the state variables including:
		%			flow_section_table;
		%			I_flow_node, I_flow_edge, I_flow_path, I_edge_path, I_flow_edge_ex, I_flow_node_ex;
		%			As_flow, As_proc, As_procz, As_load
		function onAddingFlow(this, fidx)
			%{
			if ~isempty
				this.identify_change([]);
				this.onUpateResources('request');
			%}
			assert(isrow(fidx), 'fidx should be an row vector.');
			slice = this.hs;
			Nf = slice.NumberFlows;
			num_segs = slice.NumberVNFs+1;
			Nvf = this.NumberFlowSections;
			Nsn = slice.NumberServiceNodes;
			Nvnf = slice.NumberVNFs;
			nnf = numel(fidx);   % number of new flows
			nnfs = num_segs*nnf;  % number of new flow sections
			Nl = slice.NumberLinks;
			Nn = slice.NumberNodes;
			
			this.identify_change(fidx);
			
			this.updateFlowSectionTable('add', fidx);
			if ~isempty(this.I_flow_path)
				% As the stored paths in the flow table are not the same as the candidate path, we
				% need to used 'I_flow_path' to calculate the id of candidate paths. Therefore,
				% the pid may not be continuous to existing flow ID.
				%
				% See also <SimpleDynamicSliceOptimizer>.<onAddingFlow>,
				% <DynamicSlice>.<updateFlowTable>, <SliceOptimizer>.<InitializeState>, 
				% <NormalSliceOptimizer>.<InitializeState>. Different from
				% SimpleDynamicSliceOptimizer, we do not need the data of 'I_dc_path'.
				pid = size(this.I_flow_path,2);
				slice.updateFlowTable('add', fidx, struct('pid', pid));
				for fid = fidx
					for p = 1:slice.FlowTable{fid, 'Paths'}.Width
						pid = pid + 1;
						this.I_flow_path(fid, pid) = true;
						path = slice.FlowTable{fid, 'Paths'}{p};
						for ei = 1:path.Length
							e = path.Link(ei);
							eid = slice.graph.IndexEdge(e(1),e(2));
							this.I_edge_path(eid, pid) = true;
						end
					end
				end
				this.I_flow_edge = double(this.I_flow_path)*double(transpose(this.I_edge_path));
				% this.I_flow_edge(fidx,:) = double(this.I_flow_path(fidx,:)* ...
				%			double(transpose(this.I_edge_path); assert(isequal(this.I_flow_edge(1:end-1,:), this.old_state.I_flow_edge))
				for fid = fidx
					pl = slice.FlowTable{fid, 'Paths'};	% PathList
					for j = 1:pl.Width
						this.I_flow_node(fid, pl{j}.node_list) = true;
					end
				end
			else
				this.I_flow_edge(fidx,:) = true(nnf, Nl);
				this.I_flow_node(fidx,:) = true(nnf, Nn);
			end
			Nal = this.NumberAugmentedLinks;
			Nan = this.NumberAugmentedNodes;
			this.I_flow_edge_ex(fidx, 1:Nl) = this.I_flow_edge(fidx,:);
			this.I_flow_node_ex(fidx, 1:Nn) = this.I_flow_node(fidx,:);
			for fid = fidx
				for i = (Nl+1):Nal
					[src,dest] = this.augraph.IndexEdge(i);
					if src<=Nn && this.I_flow_node(fid, src)
						this.I_flow_edge_ex(fid, i) = 1;
						this.I_flow_node_ex(fid, dest) = 1;
					end
					if dest<=Nn && this.I_flow_node(fid, dest)
						this.I_flow_edge_ex(fid, i) = 1;
						this.I_flow_node_ex(fid, src) = 1;
					end
				end
			end
			
			%% Under the condition that the topology of the slice is not changed
			% Otherwise the whole state matrices need to be updated
			aug_nid = (Nn+1):Nan;
			aug_eid = (Nl+1):Nal;
			A = this.augraph.GetIncidentMatrix;		% There might exist single direction links in the virtual topology.
			As = block_diag(A, num_segs);
			As((num_segs-1)*Nan+aug_nid, :) = 0;   % (1)
			Ap = block_diag(A, num_segs);	
			for nidx = aug_nid
				eout = find(A(nidx,:)==-1);
				for j = 1:(num_segs-1)
					As((j-1)*Nan+nidx, (j-1)*Nal+eout) = 0;
					As((j-1)*Nan+nidx,     j*Nal+eout) = -1;		% TODO:  -1 => -α : rate change
					Ap((j-1)*Nan+nidx, (j-1)*Nal+eout) = 0;
				end
			end
			% Extend the coefficent to all flows
			As_flow_temp = block_diag(As, nnf);
			As_proc_temp = block_diag(Ap, nnf);			
			As_procz_temp = speye(Nsn*Nvnf*nnf);
			% Specify the demand for the flow reservation constraints
			Ids_temp = spalloc(size(As_flow_temp,1), nnf, 2*nnf);
			nidx_offset = 0;
			for i = 1:numel(fidx)
				src = slice.FlowTable{fidx(i), 'Source'};
				Ids_temp(nidx_offset+src, i) = -1;			%#ok<SPRIX> % 2 only outgoing traffic at source
				nidx_offset = nidx_offset + (num_segs-1)*Nan;
				dest = slice.FlowTable{fidx(i), 'Target'};
				Ids_temp(nidx_offset+dest, i) = 1;			%#ok<SPRIX> % 3
				nidx_offset = nidx_offset + Nan;
			end
			this.Ids = blkdiag(this.Ids, Ids_temp);
			
			this.num_vars = [...
				Nvf * Nal;...    % number of edge variables: x(f,e)
				Nvnf*Nf*Nsn;...  % number of node variables: z(n,v,f)
				Nf ...					 % number of flows: r(f)
				];
			I_effect_varx = reshape(repmat(this.I_flow_edge_ex(fidx,:)', num_segs, 1), nnfs*Nal, 1);
			I_effect_rows = reshape(repmat(this.I_flow_node_ex(fidx,:)', num_segs, 1), Nan*nnfs, 1);
			I_unused_varx = logical(sparse(Nal*num_segs,1)); 
			idx = find(sum(As(1:Nn,aug_eid),1)==1); % 2(c) 1st segment has no traffic out from dc
			I_unused_varx(Nl+idx) = true; 			
			idx = find(sum(As((1:Nn)+Nvnf*Nan,aug_eid+Nvnf*Nal),1)==-1);% 2(d) last segment has no traffic into dc
			I_unused_varx(Nvnf*Nal+Nl+idx) = true; 
			I_unused_varx = repmat(I_unused_varx, nnf, 1);
			I_effect_varx = I_effect_varx & ~I_unused_varx;
			nidx_offset = 0;
			eidx_offset = 0;
			for fid = fidx
				src = slice.FlowTable{fid, 'Source'};
				eidx_off = find(A(src, 1:Nal)==1);		% no incoming flow at source
				I_effect_varx(eidx_offset+eidx_off) = false;
				nidx_offset = nidx_offset + (num_segs-1)*Nan;
				eidx_offset = eidx_offset + (num_segs-1)*Nal;
				dest = slice.FlowTable{fid, 'Target'};
				eidx_off = find(A(dest, 1:Nal)==-1);  % no outgoing flow at target
				I_effect_varx(eidx_offset+eidx_off) = false;
				nidx_offset = nidx_offset + Nan;
				eidx_offset = eidx_offset + Nal;
			end
			%% Fake vars
			xoffset = 0;
			I_fake_edgevars = logical(sparse(nnfs*Nal,1));
			for i =1:nnf
				for j = 1:num_segs
					I_fake_edgevars(xoffset+aug_eid) = true;
					xoffset = xoffset + Nal;
				end
			end
			b_filter_dc = logical(sparse(Nan*nnfs,1));
			nidx_offset = 0;
			for k = 1:nnf
				for j = 1:Nvnf
					b_filter_dc(nidx_offset+aug_nid) = true; 
					nidx_offset = nidx_offset + Nan;
				end
				nidx_offset = nidx_offset + Nan; % last segment: no processing function required
			end
			As_proc_temp = As_proc_temp(b_filter_dc,:);
			
			As_flow_temp(~I_effect_rows, :) = 0;
			As_flow_temp(:, ~I_effect_varx) = 0;
			As_proc_temp(~I_effect_rows(b_filter_dc), :) = 0;
			As_proc_temp(:, ~I_effect_varx) = 0;
			As_procz_temp(~I_effect_rows(b_filter_dc),:) = 0;
			
			this.As_flow = blkdiag(this.As_flow, As_flow_temp);
			this.As_proc = blkdiag(this.As_proc, As_proc_temp);
			this.As_procz = blkdiag(this.As_procz, As_procz_temp);
			this.I_fake_edgevars = [this.I_fake_edgevars; I_fake_edgevars];
			this.As_load = [repmat(speye(Nl, Nal),1, Nvf), sparse(Nl, Nsn*Nvnf*Nf);...
				sparse(Nsn, Nal*Nvf), repmat(speye(Nsn), 1, Nvnf*Nf)]; % RE-ALOOCATE
				
			%% [TODO]: sumarize the above procedure as method
			% function partialUpdate(this, action, fidx)

			%% Previous State
			% In _onAddingFlow_ and _OnRemovingFlow_, the slice topology does not change.
			% The data centers and NVF type do not changes, so there is no new/deleted VNF
			% instance variables.
			%
			% Record the state changes after flow is added to the flow table to maintain
			% the full set of variables, different from that in _<OnRemovingFlow>_.
			this.old_variables.x = sparse(Nal*Nvf,1);
			this.old_variables.z = sparse(Nvnf*Nsn*Nf,1);
			this.old_variables.r = sparse(Nf,1);
			this.old_variables.x(~this.changed_index.x) = this.old_state.variables.x; % sparse matrices with full form
			this.old_variables.z(~this.changed_index.z) = this.old_state.variables.z;
			this.old_variables.r(~this.changed_index.r) = this.old_state.variables.r;
			this.topts.old_variables = getstructfields(this.old_variables, {'x', 'z', 'r'});
			this.update_reconfig_costinfo('add');
			this.sh_options.action = 'add';
		end
		
		function changes = onRemovingFlow(this, fidx)
			slice = this.hs;
			%% Previous State
			% When removing flows, the deleted variables' value is treated as zeros, its
			% difference is constant (0-x0). Reconfiguration cost corresponding to the
			% previous state, should be updated before performing optimization.
			%
			% Record changes before flows are removed from flow table, _identify_change_
			% is called.
			this.identify_change(fidx);
			this.old_variables = structupdate(this.old_variables, this.old_state.variables, ...
				{'x', 'z', 'r'});
			this.topts.old_variables = Dictionary( ...
				'x', this.old_variables.x(~this.changed_index.x), ...
				'z', this.old_variables.z(~this.changed_index.z), ...
				'r', this.old_variables.r(~this.changed_index.r));
			this.update_reconfig_costinfo('remove');
			% Before performing optimization, there is not situation that VNF instances
			% will be removed. So there is no need to copy the |old_variable.v| like x/z.
			
			%% Update State Matrices
			this.updateFlowSectionTable('remove', fidx);
			if ~isempty(this.I_flow_path)
				this.I_flow_path(fidx, :) = [];
				idx_rm_paths = sum(this.I_flow_path,1)==0;
				this.I_flow_path(:, idx_rm_paths) = [];
				this.I_edge_path(:, idx_rm_paths) = [];
			end
			Nl = slice.NumberLinks;
			Nal = this.NumberAugmentedLinks;
			Nan = this.NumberAugmentedNodes;
			Nsn = slice.NumberServiceNodes;
			Nvnf = slice.NumberVNFs;
			num_segs = Nvnf + 1;
			rm_mat = sparse(false(size(this.I_flow_node_ex)));
			rm_mat(fidx,:) = true;
			rm_rows = reshape(repmat(rm_mat', num_segs, 1), Nan*num_segs*slice.NumberFlows, 1); 
			this.As_flow(:, this.changed_index.x) = [];  
			this.As_flow(rm_rows, :) = [];
			this.Ids(:, fidx) = [];
			this.Ids(rm_rows, :) = [];
			this.As_proc(:, this.changed_index.x) = [];
			this.As_proc(this.changed_index.z, :) = [];
			this.As_procz(this.changed_index.z, :) = [];
			this.As_procz(:, this.changed_index.z) = [];
			this.I_fake_edgevars(this.changed_index.x) = [];
			
			this.I_flow_edge(fidx, :) = [];
			this.I_flow_node(fidx, :) = [];
			this.I_flow_edge_ex(fidx, :) = [];
			this.I_flow_node_ex(fidx, :) = [];
			Nf = slice.NumberFlows - numel(fidx);  % flow has not yet been removed from flow table
			Nvf = (Nvnf+1)*Nf;		% so this.NumberFlowSections cannot be used.
			this.num_vars = [Nvf * Nal;  Nvnf*Nf*Nsn; Nf ];
			this.As_load = [repmat(speye(Nl, Nal),1, Nvf), sparse(Nl, Nsn*Nvnf*Nf);...
				sparse(Nsn, Nal*Nvf), repmat(speye(Nsn), 1, Nvnf*Nf)]; % RE-ALOOCATE
			this.sh_options.action = 'remove';
			changes.Nf = height(slice.FlowTable);
			changes.pid = slice.FlowTable{fidx(1), 'Paths'}{1}.local_id;
			%%
			% return the changes to slice to update the flow table
		end
		
		function onUpdateResources(this, action)
			switch action
				case 'release'
					%% Re-build the augmented graph
					% method-1: build from the slice.graph
					%		
					% method-2: update this.augraph
					%		1. identify the link/nodes that should be removed (which has been done by
					%		   <idnetify_change>);
					%		2. trim the augraph (see also <DirectedGraph>.<Remove>);
					%			 The node/edge removal/add must be carefully handled.
					%
					% The two methods should have the same results, as the node/link order is not
					% changed. See also <DynamicNetwork>.<createflow>
					error('error: not implemented!');
				case 'request'
					%% The state matrices and the variables (x,z,v) should be augmented
					% and the orginal ones are embedded in.
					error('error: not implemented!');
				otherwise
					error('error: not support!');
			end
			
		end
		
		%% optimalFlowRate
		%			[output, profit] = optimalFlowRate(this, options)
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
						optimalFlowRate@NormalSliceOptimizer(this, options);
					cost = cost + this.get_reconfig_cost('const');
					profit = profit - cost;
					varargout{1} = profit;
					if nargout >= 2
						varargout{2} = cost;
					end
				else
					[varargout{1:nargout}] = optimalFlowRate@NormalSliceOptimizer(this, options);
				end
			end
		end
		
		function update_options(this, options)
			update_options@IDynamicSliceOptimizer(this, options);
			update_options@NormalSliceOptimizer(this, options);
		end
		
		%%
		function s = save_state(this)
			this.old_state.augraph = this.augraph.copy();
			this.old_state.flow_section_table = this.flow_section_table;
			this.old_state.I_flow_node = this.I_flow_node;
			this.old_state.I_flow_edge = this.I_flow_edge;
			this.old_state.I_flow_node_ex = this.I_flow_node_ex;
			this.old_state.I_flow_edge_ex = this.I_flow_edge_ex;
			this.old_state.vnf_capacity = this.getVNFCapacity();
			this.old_state.I_flow_path = this.I_flow_path;
			this.old_state.I_edge_path = this.I_edge_path;
			this.old_state.As_flow = this.As_flow;
			this.old_state.As_proc = this.As_proc;
			this.old_state.As_procz = this.As_procz;
			this.old_state.Ids = this.Ids;
			this.old_state.As_load = this.As_load;
			this.old_state.I_fake_edgevars = this.I_fake_edgevars;
			this.old_state.variables = Dictionary(this.Variables);
			this.old_state.I_active_variables = this.I_active_variables;
			this.old_state.x_reconfig_cost = this.x_reconfig_cost;
			this.old_state.z_reconfig_cost = this.z_reconfig_cost;
			this.old_state.vnf_reconfig_cost = this.vnf_reconfig_cost;
			if nargout >= 1
				s = this.old_state;
			end
		end
		
		function restore_state(this, s)
      if nargin <= 1
        s = this.old_state;
      end
			this.augraph = s.old_state.augraph;
			this.flow_section_table = s.flow_section_table;
			this.I_flow_node = s.I_flow_node;
			this.I_flow_edge = s.I_flow_edge;
			this.I_flow_node_ex = s.I_flow_node_ex;
			this.I_flow_edge_ex = s.I_flow_edge_ex;
			this.I_flow_path = s.I_flow_path;
			this.I_edge_path = s.I_edge_path;
			this.As_flow = s.As_flow;
			this.As_proc = s.As_proc;
			this.As_procz = s.As_procz;
			this.Ids = s.Ids;
			this.As_load = s.As_load;
			this.I_fake_edgevars = s.I_fake_edgevars;
			this.Variables = s.variables;
			this.I_active_variables = s.I_active_variables;
			this.x_reconfig_cost = s.x_reconfig_cost;
			this.z_reconfig_cost = s.z_reconfig_cost;
			this.vnf_reconfig_cost = s.vnf_reconfig_cost;
		end
    
		function [b, vars] = postProcessing(this)
			%% post process reconfiguration
			% the super class calls <setPathBandwidth>
			postProcessing@NormalSliceOptimizer(this); % xzr
			if isfield(this.temp_vars, 'v')
				tol_zero = this.options.NonzeroTolerance;
				var_v = this.temp_vars.v;
				var_v(var_v<tol_zero*max(var_v)) = 0;
				this.Variables.v = var_v;
			end
			%%% c
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
				b_diss_z = s.diff_z_norm < tol_vec;
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
			end
			b = true;
			if nargout >= 2
				vars = this.Variables.data;
			end
		end
		
		function stat = get_reconfig_stat(this, type)
			s = this.diffstate(true);
			tol_vec = this.options.DiffNonzeroTolerance;
			slice = this.hs;
			switch type
				case 'NumVariables'
					if ~isfield(this.problem, 'I_active_xz')
						error('error: create ''I_active_xz'' when initialize problem.');
					end
					if ~isfield(this.problem, 'I_fake_vars')
						error('error: create ''I_fake_vars'' when initialize problem.');
					end
					stat = nnz(this.problem.I_active_xz(~this.problem.I_fake_vars));
					if isfield(s, 'mid_v')
						stat = stat + nnz(s.mid_v);
					end
				case 'ReVariables'
					%% Number of reconfigurations between current and previous solution
					% NOTE: resource allocation and release will not take place at the same time.
					%   Number of edge variables, VNF assignment variables the reconfiguration of
					%   VNF instances. 
					s.diff_x_norm(this.I_fake_edgevars) = 0;
					stat = nnz(s.diff_x_norm>tol_vec) + nnz(s.diff_z_norm>tol_vec); % inactive variables is 0
					if isfield(s, 'diff_v_norm')
						stat = stat + nnz(s.diff_v_norm>tol_vec);
					end
				case 'ReFlows'   % Number of Reconfigured Flows
					Nal = this.NumberAugmentedLinks;
					Nvnf = slice.NumberVNFs;
					num_segs = Nvnf + 1;
					Nf = max(slice.NumberFlows, height(slice.old_state.flow_table));
					Nsn = slice.NumberServiceNodes;
					s.diff_x_norm(this.I_fake_edgevars) = 0;
					idx_flow = sum(reshape(s.diff_x_norm>tol_vec, Nal*num_segs, Nf),1) & ...
						sum(reshape(s.diff_z_norm>tol_vec, Nsn*Nvnf, Nf),1);
					stat = nnz(idx_flow);
			end
		end
		
	end % end of Public Methods
	
	methods (Access = protected)
    [profit,cost] = fastReconfigure(this, action, options);
		% Call this function after flow is added to slice's FlowTable or before flow is
		% removed from slice's flow table.
		function updateFlowSectionTable(this, action, fidx)
			if nargin <= 1
				action = 'add';
			end
			if strcmpi(action, 'add')
				if nargin <= 2
					updateFlowSectionTable@NormalSliceOptimizer(this);
				else
					updateFlowSectionTable@NormalSliceOptimizer(this, fidx);
				end
			else
				assert(isrow(fidx), 'fidx should be an row vector.');
				num_segs = this.hs.NumberVNFs + 1;
				for fid = fidx
					this.flow_section_table{(fid-1)*num_segs+(1:num_segs), 'FlowIndex'} = 0;
				end
				this.flow_section_table(this.flow_section_table{:,'FlowIndex'}==0,:) = [];
				this.flow_section_table{:, 'FlowIndex'} = ...
					repelem((1:height(this.flow_section_table)/num_segs)', num_segs, 1);
			end
		end
		
    % Override <NormalSliceOptimizer>.<optimize>.
		% |options|: include fields minopts, CostModel, PricingPolicy (if
		%       CostModel='fixcost'), num_orig_vars;
		function [x, fval] = optimize(this, options)
			if isfield(options, 'CostModel') && strcmpi(options.CostModel, 'fixcost')
				prbm = this.problem;
				[xs, fval, exitflag, output] = ...
					fmincon(@(x)NormalDynamicSliceOptimizer.fcnSocialWelfare(x, this, options), ...
					prbm.x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, prbm.ub, [], prbm.minopts);
				if exitflag < 0
					warning([strtok(output.message, newline), 'Retry.']);
					x0 = [this.topts.old_variables.x;
						this.topts.old_variables.z;
						this.topts.old_variables.r;
						];
					if isfield(options, 'unit')
						x0 = x0(this.I_active_variables)/options.unit;
					end
					[xs, fval, exitflag, output] = ...
						fmincon(@(x)NormalDynamicSliceOptimizer.fcnSocialWelfare(x, this, options), ...
						x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, prbm.ub, [], prbm.minopts);
				end
				this.interpretExitflag(exitflag, output.message);
				if options.bCompact
					x = sparse(sum(this.num_vars), 1);
					x(this.I_active_variables) = xs;
				else
					x = xs;
				end
				constr_tol = getstructfields(prbm.minopts, {'ConstraintTolerance'}, ...
					'default-ignore', this.options.ConstraintTolerance);
				assert(this.checkFeasible(x, constr_tol), 'error: infeasible solution.');
			else
				[x, fval] = optimize@NormalSliceOptimizer(this, rmstructfields(options, 'CostModel'));
			end
    end
		
    function postl1normalizer(this)
      [~, cost] = this.get_reconfig_cost('const');
      [~, linear_cost] = this.get_reconfig_cost('linear', true);
      field_names = {'const', 'linear'};
      b = this.getbeta();
      if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
        arr_cost = [cost, linear_cost];
        fn = fieldnames(arr_cost(2));
        for j = 1:length(fn)
          arr_cost(2).(fn{j}) = arr_cost(2).(fn{j})/b.(fn{j});
        end
        for i = 1:2
          if length(this.raw_cost.(field_names{i})) < this.NUM_MEAN_BETA
            this.raw_cost.(field_names{i})(end+1) = ...
              getstructfields(arr_cost(i), {'x','z'});
          else
            this.raw_cost.(field_names{i}) = ...
              [this.raw_cost.(field_names{i})(2:end), ...
              getstructfields(arr_cost(i), {'x','z'})];
          end
        end
        if isfield(arr_cost, 'v')
          if arr_cost(1).v < 10^-3 || arr_cost(2).v < 10^-3
            return;
          end
          if this.raw_costv.bInit
            this.raw_costv.const = arr_cost(1).v;
            this.raw_costv.linear = arr_cost(2).v;
            this.raw_costv.bInit = false;
            return;
          end
          for i = 1:2
            if length(this.raw_costv.(field_names{i})) < this.NUM_MEAN_BETA
              this.raw_costv.(field_names{i})(end+1) = arr_cost(i).v;
            else
              this.raw_costv.(field_names{i}) = ...
                [this.raw_costv.(field_names{i})(2:end), arr_cost(i).v];
            end
          end
        end
      elseif strcmpi(this.GET_BETA_METHOD, 'OverallAverage')
        
      else
        if this.b_dim
          field_names = {'x','z','v'};
        else
          field_names = {'x','z'};
        end
        for i = 1:length(field_names)
          if isfield(cost, field_names{i}) && isfield(linear_cost, field_names{i})
            if cost.(field_names{i})<10^-6
              continue;
            else
              % If dierctly set to f, the real (const) cost and the estimated (linear)
              % cost are equal, providing that the optimal solution does not change.
              % However, if |f>1|, the new solution result in less changes (with drastic
              % value changes), which drastically decreases the const cost, leading to the
              % const cost and linear cost imbalance again, so we guess |sqrt(f)| is a
              % relatively good factor. Similar condition applies to |f<1|.
              f = cost.(field_names{i})/linear_cost.(field_names{i});
              beta.(field_names{i}) = sqrt(f);
            end
            if length(this.raw_beta.(field_names{i})) < this.NUM_MEAN_BETA
              this.raw_beta.(field_names{i})(end+1) = beta.(field_names{i});
            else
              this.raw_beta.(field_names{i}) = ...
                [this.raw_beta.(field_names{i})(2:end), beta.(field_names{i})];
            end
          end
        end
      end
    end
    
    function b = getbeta(this)
      field_names = {'x','z','v'};
      if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
        for i=1:2
          b.(field_names{i}) = sum([this.raw_cost.const.(field_names{i})].*[this.raw_cost.linear.(field_names{i})])./ ...
            sum([this.raw_cost.linear.(field_names{i})].^2);
        end
        b.v = sum(this.raw_costv.const.*this.raw_costv.linear)./ ...
          sum(this.raw_costv.linear.^2);
      elseif strcmpi(this.GET_BETA_METHOD, 'OverallAverage')
        
      else
        if strcmpi(this.GET_BETA_METHOD, 'ExponentiaLMovingAverage')
          %           a0 = 0.25;
          %           for i=1:2
          %             n = length(this.raw_beta.(field_names{i}));
          %             alpha = zeros(n,1);
          %             alpha(2:n) = a0 * (1-a0).^(n-2:-1:0);
          %             alpha(1) = (1-a0)^(n-1);
          %             b.(field_names{i}) = dot(alpha, this.raw_beta.(field_names{i}));
          %           end
          a = 0.55;
          for i=1:2
            n = length(this.raw_beta.(field_names{i}));
            a0 = (1-a)/(1+a^n);
            alpha = a0*a.^(n-1:-1:0);
            b.(field_names{i}) = dot(alpha, this.raw_beta.(field_names{i}));
          end
        else
          for i=1:2
            b.(field_names{i}) = mean(this.raw_beta.(field_names{i}));
          end
        end
        % v_nf >= sum_{p}{z_npf}
        n_path = nnz(this.I_dc_path)/this.hs.NumberServiceNodes;
        b.v = b.z/n_path;
      end
    end
		
		function ds = diffstate(this, isfinal)
			ds = diffstate@IDynamicSliceOptimizer(this, isfinal);
			
			%% [TODO]
			
			if isfinal
				this.diff_state = ds;
			end
		end
		
		function identify_change(this, changed_flow_index)
			global DEBUG; %#ok<NUSED>
			if islogical(changed_flow_index)
				changed_flow_index = find(changed_flow_index);
			end
			Nal = this.NumberAugmentedLinks;
			slice = this.hs;
			num_segs = slice.NumberVNFs+1;
			% NOTE: This function should be called after adding new flow, or before removing
			% departing flow. Therefore 'slice.NumberFlows' has consistent meaning here.
			Nf = slice.NumberFlows;
			Nvf = this.NumberFlowSections;  % num_segs*NumberFlows
			Nvnf = slice.NumberVNFs;
			Nsn = slice.NumberServiceNodes;
			this.changed_index.x = logical(sparse(Nal*Nvf,1));
			for fid = changed_flow_index
				this.changed_index.x(Nal*num_segs*(fid-1)+(1:Nal*num_segs)) = true;
			end
			this.changed_index.z = logical(sparse(Nsn*Nvnf*Nf,1));
			for fid = changed_flow_index
				this.changed_index.z(Nsn*Nvnf*(fid-1)+(1:Nsn*Nvnf)) = true;
			end
			this.changed_index.r = logical(sparse(Nf,1));
			this.changed_index.r(changed_flow_index) = true;
			if ~isempty(slice.net_changes)
				%% Situations when there is network updates
				% (1) release resource description after optimization;
				%	(2) After ad-hoc flow arrives, pre-processing (onAddingFlow) before optimization.
				if isfield(slice.net_changes, 'LinkIndex')
					idx_link = slice.net_changes.LinkIndex;
					link_offset = 0;
					for i = 1:Nvf
						this.changed_index.x(link_offset+idx_link) = true;
						link_offset = link_offset + Nal;
					end
				end
				if isfield(slice.net_changes, 'DCIndex')
					dc_idx = slice.net_changes.DCIndex;  % the co-located node's index
					% identifies the augmented links that is changed due to topolgy changes
					Nn = slice.NumberNodes;
					idx_link = ismember(this.augraph.Head, dc_idx+Nn) | ismember(this.augraph.Tail, dc_idx+Nn);
					link_offset = 0;
					for i = 1:Nvf
						this.changed_index.x(link_offset+idx_link) = true;
						link_offset = link_offset + Nal;
					end
					node_offset = 0;
					for i = 1:Nvnf*Nf
						this.changed_index.z(node_offset+dc_idx) = true;
						node_offset = node_offset + Nsn;
					end
					this.changed_index.v = logical(sparse(Nsn*Nvnf,1));
					node_offset = 0;
					for i = 1:Nf
						this.changed_index.v(node_offset+dc_idx) = true;
						node_offset = node_offset + Nsn;
					end
				end
			end
		end
		

	end % end protected method
	
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
				this.x_reconfig_cost = zeros(this.num_vars(1),1);
				this.x_reconfig_cost(~this.I_fake_edgevars) = ...
					repmat(slice.Links.ReconfigCost, slice.NumberFlows*(slice.NumberVNFs+1), 1);
				this.z_reconfig_cost = repmat(slice.ServiceNodes.ReconfigCost, ...
					slice.NumberFlows*slice.NumberVNFs, 1);
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
					beta = this.getbeta();
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
				old_number_flow = height(slice.old_state.flow_table);
				old_num_varx = numel(slice.old_state.variables.x);
				this.x_reconfig_cost = zeros(old_num_varx,1);
				this.x_reconfig_cost(~this.old_state.I_fake_edgevars) = ...
					repmat(slice.Links.ReconfigCost, old_number_flow*(slice.NumberVNFs+1), 1);
				this.z_reconfig_cost = repmat(slice.ServiceNodes.ReconfigCost, ...
					old_number_flow*slice.NumberVNFs, 1);
				%%
				% The removed flow's reconfiguration cost is not considered in the
				% optimization, which is constant.
				if this.ENABLE_DYNAMIC_NORMALIZER
					beta = this.getbeta();
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
		% (super-linear). See also <fcnHessian>, <NormalSliceOptimizer>.<fcnProfitPrimalCC>.
		function [profit, grad] = fcnFastConfigProfit(vars, this, options)
			if isfield(options, 'unit')
				unit = options.unit;
			else
				unit = 1;
			end
			num_vars = length(vars);
			prbm = this.problem;
			% |tx|,|tz| and |tv| are auxilliary variables to transform L1 norm to
			% linear constraints and objective part. The number of auxiliary variables
			% eqauls to the number of the main variables.
			if length(prbm.num_vars) >= 6  % for _fastReconfigure2_, including |v|, and |tv|
				error('error: not implemented!');
			end
			
			num_varxz = sum(prbm.num_vars(1:2));
			num_vars_prim = sum(prbm.num_vars(1:3));
			flow_rate = vars((num_varxz+1):num_vars_prim);
			num_vartx = prbm.num_vars(4);			% the number of x and tx after compress is not the same
			var_tx_index = num_vars_prim+(1:num_vartx);
			var_tz_index = num_vars_prim+num_vartx+(1:prbm.num_vars(5));
			var_tx = vars(var_tx_index);
			var_tz = vars(var_tz_index);
			%% objective value
			slice = this.hs;
			profit = -slice.Weight*sum(fcnUtility(unit*flow_rate));
			%%%
			% calculate reconfiguration cost by |tx| and |tz| (L1 approximation),
			% which is similar to calculating by |x-x0| and |z-z0|.
			profit = profit + dot(var_tx, this.topts.x_reconfig_cost*unit) + ...
				dot(var_tz, this.topts.z_reconfig_cost*unit);
			if length(prbm.num_vars) >= 6 % for _fastConfigure2_
				error('error: not implemented!');
			end
			% If there is only one output argument, return the real profit (positive)
			if options.bFinal
				profit = -profit;
			else
				%% Gradient
				% The partial derivatives are computed by dividing the variables into five parts,
				% i.e., $x,z,r,t_x,t_z$. Since |x,z| do not appear in objective function, the
				% corrsponding derivatives is zero.  
				%
				grad = sparse(num_vars, 1);
				% The partial derivatives of r
				grad((num_varxz+1):num_vars_prim) = -slice.Weight./(1/unit+flow_rate);
				%%%
				% The partial derivatives of t_x is the vector |x_reconfig_cost|;
				% The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
				% since the node reconfiguration cost is only depend on node.
				grad(var_tx_index) = this.topts.x_reconfig_cost*unit;
				grad(var_tz_index) = this.topts.z_reconfig_cost*unit;
				if length(prbm.num_vars) >= 6
					error('error: not implemented!');
				end
			end
		end
		
		%%
		% See also <NormalSliceOptimzier.fcnHessianCC>.
		% The objective function only has non-zero second derivatives on |r|.
		function h = hessReconfigure(vars, lambda, this, options) %#ok<INUSL>
			if nargin>=4 && isfield(options, 'unit')
				unit = options.unit;
			else
				unit = 1;
			end
			num_varxz = sum(this.problem.num_vars(1:2));
			num_vars_prim = sum(this.problem.num_vars(1:3));
			var_r_index = (num_varxz+1):num_vars_prim;
			h = sparse(numel(vars), numel(vars));
			flow_rate = vars(var_r_index);
			h(var_r_index, var_r_index) = spdiag(this.hs.Weight./(1/unit+flow_rate).^2);
		end
		
	end
	
end % end of class