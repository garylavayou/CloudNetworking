classdef NormalSliceOptimizer < SliceOptimizer
	
	properties
		%% Augmented Topology
		augraph DirectedGraph;
		%% Virtual Flow Table
		% We need to segment each flow into multiple segements
		%% TODO: remove this variable, no specific usage.
		flow_section_table;
	end
	
	properties
		%% Incidence Matrices
		I_node_path logical;
		% Excluding the edges/nodes that are not on the candidate paths
		% (Virtual) Edge-Flow Incidence Matrix
		I_flow_edge logical;
		% (Augmented Virtual) Edge-Flow Incidence Matrix
		I_flow_edge_ex logical;
		% (Virtual) Node-Flow Incidence Matrix: the nodes include
		% farwarding nodes and DataCenters
		I_flow_node logical;
		% (Augmented Virtual) Node-Flow Incidence Matrix
		I_flow_node_ex logical;
		%
		I_fake_edgevars;
		
		%% Problem coefficients
		As_flow double;		% coefficient for flow reservation constraints
		Ids double;
		As_proc double;			% processing resource requirements
		As_procz double;
		As_load double;
	end
	
	properties (Dependent)
		NumberAugmentedNodes;
		NumberAugmentedLinks;
		NumberFlowSections;
	end
	
	%% Constructor
	methods
		function this = NormalSliceOptimizer(slice, options)
			if nargin >= 2
				args = {slice, options};
			elseif nargin == 1
				args = {slice};
			else
				args = {};
			end
			this@SliceOptimizer(args{:});
			
			this.problem = struct('Aeq', [], 'A', [], 'beq', [], 'b', []);
			this.flow_section_table = table([], [], [], ...
				'VariableNames', {'FlowIndex', 'Source', 'Destination'});
		end
	end
	
	%% Property Get Methods
	methods
		function n = get.NumberAugmentedNodes(this)
			n = this.augraph.NumberNodes;
		end
		
		function n = get.NumberAugmentedLinks(this)
			n = this.augraph.NumberEdges;
		end
		
		function n = get.NumberFlowSections(this)
			n = this.hs.NumberFlows * (this.hs.NumberVNFs+1);
		end
	end
	
	%% Public Methods
	methods
		[profit, cost, output] = optimalFlowRate(this, options);
		[output, loads, fval] = priceOptimalFlowRate(this, x0, options);
		
		% explictly call this function.
		function [array] = initializeState(this, options)
			initializeState@SliceOptimizer(this);
			%% Build the augmented topology
			slice = this.hs;
			Nn = slice.NumberNodes;
			Nl = slice.NumberLinks;
			Nf = slice.NumberFlows;
			Nvnf = slice.NumberVNFs;
			Nsn = slice.NumberServiceNodes;
			num_segs = Nvnf + 1;
			Nvf = this.NumberFlowSections;
			this.augraph = slice.graph.copy;
			aug_vn = Nn + (1:Nsn)';
			aug_heads = [slice.ServiceNodes.VirtualNode; aug_vn];
			aug_tails = [aug_vn; slice.ServiceNodes.VirtualNode];
			num_aug_links = length(aug_heads);
			props.Weight = eps*ones(num_aug_links,1);
			props.Capacity = inf*ones(num_aug_links,1);
			this.augraph.Update(aug_heads, aug_tails, props);		% TODO: remap
			Nan = this.NumberAugmentedNodes;
			Nal = this.NumberAugmentedLinks;
			
			%% Build virtual flow table
			% Fields include: flow index, source, destination, etc.
			if Nf > 0
				this.updateFlowSectionTable();
			end
			
			%% Initialize edge-flow & node-flow incidence matrix
			this.I_flow_node = sparse(Nf, Nn);
			if ~isempty(this.I_flow_path)
				this.I_flow_edge = double(this.I_flow_path)*double(transpose(this.I_edge_path));
				for i = 1:Nf
					pl = slice.FlowTable{i, 'Paths'};	% PathList
					for j = 1:pl.Width
						this.I_flow_node(i, pl{j}.node_list) = true;
					end
				end
			else
				this.I_flow_edge = true(Nf, Nl);
				this.I_flow_node = true(Nf, Nn);
			end
			this.I_flow_edge_ex = this.I_flow_edge;
			this.I_flow_edge_ex(Nf, Nal) = 0;
			this.I_flow_node_ex = this.I_flow_node;
			this.I_flow_node_ex(Nf, Nan) = 0;
			for k = 1:Nf
				for i = (Nl+1):Nal
					[src,dest] = this.augraph.IndexEdge(i);
					if src<=Nn && this.I_flow_node(k, src)
						this.I_flow_edge_ex(k, i) = 1;
						this.I_flow_node_ex(k, dest) = 1;
					end
					if dest<=Nn && this.I_flow_node(k, dest)
						this.I_flow_edge_ex(k, i) = 1;
						this.I_flow_node_ex(k, src) = 1;
					end
				end
			end

			%% The Flow Reservation Law and Processing Requirement Constraints
			% Serveral point must be noted:
			% (1) The last segment does not go into processing nodes, so there is no flow
			%			reservation constraint on the processing nodes for the last segment.
			%			Hence, total number of constrants for single flow is: Nan*(num_segs-1)+Nn;
			%			For ease of presentation, we use block diag to generate the matrix, while the
			%			last (Nan-Nn) rows is not needed. But we keep it for easy reference of
			%			matrices dimensions.
			%
			% (2) On the other hand, only the first |Nvnf| segments should be processed by the
			%			DC nodes, so the processing requirement constraints only contains |Nvnf| part.
			%			And the constraints is only counted on processing nodes, so the total number
			%			of constraints for single flow is: (Nan-Nn)*Nvnf
			%
			%	(3) The first segment should not come from DC nodes, and the last segment should
			%			not go into DC nodes. This observation make some edge variables to be 0, and
			%			thus can be elimenated. This does not influence the VNF allocation variables,
			%			as it related to the first |Nvnf| segments' incoming traffic and the last
			%			|Nvnf| segments' outgoing traffic on DC nodes. (We eliminate the edge
			%			variables at a later stage) 
			%% TODO
			% VNF type constraint: if a type π of VNF cannot be created on a node i, the
			% corresponding variable z_i(π)=0. Thus, for all flows, the related edge variables
			% x_ji(f,π-1)=0 and x_ij(f,π)=0. 
			%			Find the row in |As| that corresponds to the type π on aug-node n, label the
			%			associated link-variables to be removed, and clear the row. Since this
			%			constraint applies to all flows, perform it before |As| duplicated. Similar
			%			operation is also performed on |Ap|.
			%
			% Flow-VNF assoication constraint: a flow requires that a type of VNF is deployed on
			% specified nodes. As a result, we eliminate those un-related variables z_i(f,π)=0,
			% also indicating x_ji(f,π-1)=0 and x_ij(f,π)=0.
			%			Similar to the first constraint, but this constraint is flow-specific. Hence,
			%			the operation is performed to each flow on the |As_flow_temp|, |As_proc_temp|.
			%
			% Finally, clear the coefficients in |As_proc_z| corresponding to the VNF variables. 
			aug_nid = (Nn+1):Nan;
			aug_eid = (Nl+1):Nal;
			A = this.augraph.GetIncidentMatrix;		% There might exist single direction links in the virtual topology.
			As = block_diag(A, num_segs);
			As(Nvnf*Nan+aug_nid, :) = 0;   % (1)
			Ap = block_diag(A, num_segs);	
			for nidx = aug_nid
				eout = find(A(nidx,:)==-1);
				for j = 1:Nvnf
					As((j-1)*Nan+nidx, (j-1)*Nal+eout) = 0;
					As((j-1)*Nan+nidx,     j*Nal+eout) = -1;		% TODO:  -1 => -α : rate change
					Ap((j-1)*Nan+nidx, (j-1)*Nal+eout) = 0;
				end
			end
			% Extend the coefficent to all flows
			As_flow_temp = block_diag(As, Nf);
			As_proc_temp = block_diag(Ap, Nf);			
			% Specify the demand for the flow reservation constraints
			Ids_temp = spalloc(size(As_flow_temp,1), Nf, 2*Nf);
			nidx_offset = 0;
			for k = 1:Nf
				src = slice.FlowTable{k, 'Source'};
				Ids_temp(nidx_offset+src, k) = -1;			%#ok<SPRIX> % 2 only outgoing traffic at source
				nidx_offset = nidx_offset + (num_segs-1)*Nan;
				dest = slice.FlowTable{k, 'Target'};
				Ids_temp(nidx_offset+dest, k) = 1;			%#ok<SPRIX> % 3
				nidx_offset = nidx_offset + Nan;
			end
			this.Ids = Ids_temp;
			this.As_flow = As_flow_temp;
			this.As_proc = As_proc_temp;
						
			%% Compact
			% 1. Limit each flow's distributing scope based on its candidate path. 
			%		 By default, each flow potentially uses all edges and nodes. As a result, the
			%		 coefficient matrices that are generated from the incidence matrix are with full
			%		 column rank. However, since we pre-clculated candiate paths the flow's scope is
			%		 limited, as indicated by 'I_flow_edge_ex' and 'I_flow_node_ex'. Hence, we can
			%		 utilize the two indicator to eliminate unused variables, and reduce the
			%		 coefficient matrix.
			% 2. Eleminate variables that is predetermined to be 0, which including:
			%		a. incoming edge at source node: source node has no incoming traffic, indicating
			%			 the incoming edge not used;
			%		b. outgoing flow at target node: same as (a);
			%   c. first segment has no outgoing traffic from dc node;
			%		d. last segment has no incoming traffic to dc node;
			%		e. incoming/outgoing flow at the node that cannot create a type
			%			 of VNF. [TODO]
			%
			% The basic variables are given below. The number of variables in each optimization
			% shceme might be different, but they have the common basic variables. The specific
			% optimization scheme will update the information.
			this.num_vars = [...
				Nvf * Nal;...    % number of edge variables: x(f,e)
				Nvnf*Nf*Nsn;...  % number of node variables: z(n,v,f): (v,f)'th section's processing requirement on each service node.
				Nf ...					 % number of flows: r(f)
				];
			I_effect_varx = reshape(repmat(this.I_flow_edge_ex', num_segs, 1), this.num_vars(1), 1);  % (1)
			I_effect_rows = reshape(repmat(this.I_flow_node_ex', num_segs, 1), Nan*Nvf, 1);
			I_unused_varx = logical(sparse(Nal*num_segs,1)); 
			idx = find(sum(As(1:Nn,aug_eid),1)==1); % 2(c) 1st segment has no traffic out from dc
			I_unused_varx(Nl+idx) = true; 			
			idx = find(sum(As((1:Nn)+Nvnf*Nan,aug_eid+Nvnf*Nal),1)==-1);% 2(d) last segment has no traffic into dc
			% Equivalent to:	find(sum(A(Nn+1:Nan,aug_eid),1)==1)																			
			I_unused_varx(Nvnf*Nal+Nl+idx) = true; 
			I_unused_varx = repmat(I_unused_varx, Nf, 1);
			I_effect_varx = I_effect_varx & ~I_unused_varx;
			nidx_offset = 0;
			eidx_offset = 0;
			for k = 1:Nf
				src = slice.FlowTable{k, 'Source'};
				eidx_off = find(A(src, 1:Nal)==1);
				I_effect_varx(eidx_offset+eidx_off) = false;  % 2(a) no incoming flow at source
				nidx_offset = nidx_offset + (num_segs-1)*Nan;
				eidx_offset = eidx_offset + (num_segs-1)*Nal;
				dest = slice.FlowTable{k, 'Target'};
				eidx_off = find(A(dest, 1:Nal)==-1); 
				I_effect_varx(eidx_offset+eidx_off) = false;  % 2(b) no outgoing flow at target
				nidx_offset = nidx_offset + Nan;
				eidx_offset = eidx_offset + Nal;
			end
			%% Fake vars
			this.I_fake_edgevars = logical(sparse(this.num_vars(1),1));
			xoffset = 0;
			for i =1:Nf
				for j = 1:num_segs
					this.I_fake_edgevars(xoffset+aug_eid) = true;
					xoffset = xoffset + Nal;
				end
			end
			b_filter_dc = logical(sparse(Nan*Nvf,1));
			nidx_offset = 0;
			for k = 1:Nf
				for j = 1:Nvnf
					b_filter_dc(nidx_offset+aug_nid) = true; 
					nidx_offset = nidx_offset + Nan;
				end
				nidx_offset = nidx_offset + Nan; % last segment: no processing function required
			end
			this.As_proc = this.As_proc(b_filter_dc,:);
			%% Reduce the matrices
			% We do not remove the cleared coefficient, so that the matrices can be easily used
			% for its known dimensions.
			this.As_flow(~I_effect_rows, :) = 0;
			this.As_flow(:, ~I_effect_varx) = 0;
			this.As_proc(~I_effect_rows(b_filter_dc),:) = 0;
			this.As_proc(:, ~I_effect_varx) = 0;
			this.As_procz = speye(Nsn*Nvnf*Nf);
			this.As_procz(~I_effect_rows(b_filter_dc),:) = 0;
			this.As_load = [repmat(speye(Nl, Nal),1, Nvf), sparse(Nl, Nsn*Nvnf*Nf);...
				sparse(Nsn, Nal*Nvf), repmat(speye(Nsn), 1, Nvnf*Nf)]; % for load-capacity constraints
			%% NOTE: used to pass data in parallel computing
			if nargout >= 1
				array = Dictionary('augraph', this.augraph, ...
					'I_active_variables', this.I_active_variables, ...
					'As_load', this.As_load, ...
					'Ids', this.Ids, ...
					'As_flow', this.As_flow,...
					'As_proc', this.As_proc, ...
					'As_procz', this.As_procz, ...
					'I_flow_node', this.I_flow_node, ...
					'I_flow_edge', this.I_flow_edge, ...
					'I_flow_node_ex', this.I_flow_node_ex, ...
					'I_flow_edge_ex', this.I_flow_edge_ex ...
					);
			end
		end
		
		% function initializeParallel(this, options)

		function setProblem(this, varargin) 
			setProblem@SliceOptimizer(this, varargin{:});
			for i = 1:2:(length(varargin)-1)
				switch varargin{i}
					case 'Array'
						% 						this.augraph = varargin{i+1}.augraph;
						% 						this.As_flow = varargin{i+1}.As_flow;
						% 						this.As_proc = varargin{i+1}.As_proc;
						% 						this.As_procz = varargin{i+1}.As_procz;
						% 						this.As_load = varargin{i+1}.As_load;
						% 						this.Ids = varargin{i+1}.Ids;
						% 						this.I_flow_node = varargin{i+1}.I_flow_node;
						% 						this.I_flow_edge = varargin{i+1}.I_flow_edge;
						% 						this.I_flow_node_ex = varargin{i+1}.I_flow_node_ex;
						% 						this.I_flow_edge_ex = varargin{i+1}.I_flow_edge_ex;
						error('error: not implemented.');
				end
			end
		end
		
		function b = checkFeasible(this, vars, opt_opts)
			b = true;
		end
		
		%%
		% Implement <SliceOptimizer.getFlowRate>
		function r = getFlowRate(this, ~)
			r = zeros(this.hs.NumberFlows, 1);
			for i = 1:this.hs.NumberFlows
				% If a flow does not have available paths, 'Width' is 0.
				for j = 1:this.hs.FlowTable{i, 'Paths'}.Width
					r(i) = r(i) + this.hs.FlowTable{i, 'Paths'}{j}.bandwidth;
				end
			end
		end
		
		function update_options(this, options)
			update_options@SliceOptimizer(this, options);
			fields = fieldnames(options);
			for i = 1:length(fields)
				switch options
					case 'CapacityConstrained'
						this.options.CapacityConstrained = options.CapacityConstrained;
				end
			end
		end
		
		function profit = getProfit(this, options)
			% determine varriables.
			if nargin <= 1
				options = Dictionary();
			else
				options = Dictionary(options);
			end
			if this.hs.isFinal()
				vars = [this.Variables.x; this.Variables.z; this.Variables.r];
			else
				vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.r];
			end
			vars = vars(this.I_active_variables);
			if isfield(options, 'LinkPrice')
				this.prices.Link = options.LinkPrice;
			elseif this.hs.isFinal()
				this.prices.Link = this.hs.Links.Price;
			end
			if isfield(options, 'NodePrice')
				this.prices.Node = options.NodePrice;
			elseif this.hs.isFinal()
				this.prices.Node= this.hs.ServiceNodes.Price;
			end
			options = setdefault(options, this.hs.options, {'PricingPolicy'});
			%%
			% Use class name to avoid dynamic loading of class.
			% Subclasses may override this method to define different ways to calculate profits.
			options.bFinal = true;
			options.bCompact = false;
			profit = NormalSliceOptimizer.fcnProfitPrimalCC(vars, this, options);
			this.setProblem('Price', []);
		end
		
		% 		function saveTempResults(this, result)
		% 			saveTempResults@SliceOptimizer(this, result);
		% 			% FOR DEBUG
		% 			% this.setPathBandwidth(this.temp_vars.x);
		% 		end
		
		function vc = getVNFCapacity(this)
			if ~isfield(this.Variables, 'v') || isempty(this.Variables.v)
				error('error: cannot retrive VNF capacity, the variables not initilized yet.');
			else
				vc = this.Variables.v;
			end
		end
		
		function lc = getVNFLoad(this, z)
			if nargin <= 1
				z = this.Variables.z;
			end
			% sum variable z(n,v,f) on f result in v(n,f)
			slice = this.hs;
			lc = sum(reshape(z, slice.NumberServiceNodes*slice.NumberVNFs, slice.NumberFlows),2);
		end
		
		%%
		% Override <Slice.getLinkLoad>
		function ye = getLinkLoad(this, isfinal, vars)
			% edge_vars:
			if nargin >= 3
				if isnumeric(vars)
					edge_vars = vars;
				else
					edge_vars = vars.x;
				end
			else
				if nargin == 1 || isfinal
					edge_vars = this.Variables.x;
				elseif nargin >= 2
					edge_vars = this.temp_vars.x;
				end
			end
			
			Nve = this.hs.NumberLinks;
			num_vars_edge = this.num_vars(1);
			As = this.As_load(1:Nve, 1:num_vars_edge);
			if length(edge_vars) == num_vars_edge
				ye = full(As * edge_vars);
			else
				b_active_edges = ~this.I_fake_edgevars;
				if length(edge_vars) == nnz(b_active_edges)
					full_edge_vars = sparse(num_vars_edge, 1);
					full_edge_vars(b_active_edges) = edge_vars;
					edge_vars = full_edge_vars;
				end
				ye = full(As * edge_vars(this.I_active_variables(1:num_vars_edge)));
			end
		end
		
		%%
		% Override <Slice.getNodeLoad>
		% |isfinal|: empty is not allowed.
		function vn = getNodeLoad(this, isfinal, vars)
			if nargin >= 3
				if isnumeric(vars)
					node_vars = vars;
				else
					node_vars = vars.z;
				end
			else
				if nargin == 1 || isfinal
					node_vars = this.Variables.z;
				elseif nargin >= 2
					node_vars = this.temp_vars.z;
				end
			end
			
			Nsn = this.hs.NumberServiceNodes;
			Nve = this.hs.NumberLinks;
			num_vars_node = this.num_vars(2);
			num_vars_edge = this.num_vars(1);
			nx = this.num_vars(1);
			nz = this.num_vars(2);
			if length(node_vars) > nz
				node_vars = node_vars(this.I_active_variables(num_vars_edge+(1:num_vars_node)));
			end
			As = this.As_load(Nve+(1:Nsn),nx+(1:nz));
			vn = full(As*node_vars);
		end
		
		function [tf, vars] = postProcessing(this, ~)
			slice = this.hs;
			num_segs = slice.NumberVNFs+1;
			Nan = this.NumberAugmentedNodes;
			Nal = this.NumberAugmentedLinks;
			Nf = slice.NumberFlows;
			nidx_offset = 0;
			eidx_offset = 0;
			tol_zero = this.options.NonzeroTolerance;
			maxx = max(this.temp_vars.x);
      maxz = max(this.temp_vars.z);
      this.temp_vars.x(this.temp_vars.x<tol_zero*maxx) = 0;
      this.temp_vars.z(this.temp_vars.z<tol_zero*maxz) = 0;
      x = spalloc(this.num_vars(1),1, nnz(this.temp_vars.x));
      for fid = 1:Nf
        As = cell(num_segs,1);
        idx_fx = eidx_offset+(1:Nal*num_segs);
        fx_mask = this.I_fake_edgevars(idx_fx);
        fx = this.temp_vars.x(idx_fx);
        gen_paths = struct('idx_edgevar', [], 'nodes', [], 'bandwidth', []);
        p = 0;
        while true
          t_nidx_offset = 0;
          t_eidx_offset = 0;
          %% Copy the incident matrix for the flow
          for j = 1:num_segs
            As{j} = this.augraph.GetIncidentMatrix();
            %% NOTE: As_flow is not appropriate
            % As it should handle the flow conservation law at the Virtual Nodes, the
            % incident coefficient in As_flow has been interchanged.
            % See <initializeState>.
            %
            % As{j} = this.As_flow(nidx_offset+(j-1)*Nan+(1:Nan), ...
            %  eidx_offset+(j-1)*Nal+(1:Nal));		% As{j} are equal for all j
            eidx = fx(t_eidx_offset+(1:Nal))<eps;
            As{j}(:,eidx) = 0;
            t_eidx_offset = t_eidx_offset + Nal;
            t_nidx_offset = t_nidx_offset + Nan;
          end
          src = slice.FlowTable{fid, 'Source'};
          %% build a single path edge-by-edge
          nodes = [src, zeros(1,19)];
          idx_varedge = zeros(1,20);
          t_eidx_offset = 0;
          nidx = 1;
          eidx = 0;
          for j = 1:num_segs
            while true
              eout = find(As{j}(src,:)==-1);
              if ~isempty(eout)
                [~,idx] = max(fx(eout+t_eidx_offset));
                [n_cur, n_next] = this.augraph.IndexEdge(eout(idx));
                assert(n_cur==src);
                nidx = nidx + 1; nodes(nidx) = n_next;
                eidx = eidx + 1; idx_varedge(eidx) = eout(idx)+t_eidx_offset;
                As{j}(src,eout(idx)) = 0;
                src = n_next;
                if src > slice.NumberNodes
                  break;
                end
              else
                break;
              end
            end
            t_eidx_offset = t_eidx_offset + Nal;
          end
          nodes = nodes(1:nidx);
          idx_varedge = idx_varedge(1:eidx);
          if nodes(end) ~= slice.FlowTable{fid, 'Target'}
            break;  % this is not a complete path, and no more paths can be generated.
          end
          % 					idx_aug_nodes = find(path.nodes>Nn);
          % 					idx_aug_edges = reshape([idx_aug_nodes-1,idx_aug_nodes]', ...
          % 						2*length(idx_aug_nodes),1);
          % 					path.idx_simple_edgevar = path_idx_edgevar;
          % 					path.idx_simple_edgevar(idx_aug_edges) = [];
          p = p + 1; gen_paths(p) = struct(...
            'idx_edgevar', idx_varedge, ...
            'nodes', nodes, ...
            'bandwidth', min(fx(idx_varedge)) ...
            );
          % edge_vars on all links (including augmented links) are effective.
          fx(idx_varedge) = fx(idx_varedge) - gen_paths(p).bandwidth;
          fx(fx<tol_zero*maxx) = 0;
          if isempty(find(fx(~fx_mask)>=tol_zero*maxx,1))
            break;
          end
        end
        if isempty(gen_paths(1).nodes)
          warning('No bandwidth allocated to flow %d, create a dummy path.', fid);
          slice.FlowTable{fid, 'Paths'} = PathList(slice.FlowTable{fid, 'Paths'}{1});
          slice.FlowTable{fid, 'Paths'}{1}.bandwidth = 0;
        else
          bd = [gen_paths.bandwidth];
          [~, pidx] = sort(bd, 'descend');
          slice.FlowTable{fid, 'Paths'} = PathList();
          num_paths = min(p,slice.options.NumberPaths);
          for i = 1:num_paths
            p = pidx(i);
            slice.FlowTable{fid, 'Paths'}.Add(Path(gen_paths(p).nodes, gen_paths(p).bandwidth));
            idx = gen_paths(p).idx_edgevar + eidx_offset;
            x(idx) = x(idx) + gen_paths(p).bandwidth; %#ok<SPRIX>
          end
        end
        
        eidx_offset = eidx_offset + Nal*num_segs;
        nidx_offset = nidx_offset + Nan*num_segs;
      end
      
      %% Update Flow-Edge variables, Flow-Node variables and VNF capacity variables
      % After rounding procedure, update the variables.
      link_load = this.getLinkLoad(false);
      vnf_load = this.getVNFLoad(this.temp_vars.z);
      z = this.As_proc * x;
      link_load_r = this.getLinkLoad(false, x);
      vnf_load_r = this.getVNFLoad(z);
      idx_links = (link_load~=0) & (link_load_r~=0);
      idx_vnfs = (vnf_load~=0) & (vnf_load_r~=0);
      af = min([link_load(idx_links)./link_load_r(idx_links); ...
        vnf_load(idx_vnfs)./vnf_load_r(idx_vnfs)]);
      this.Variables.x = af*x;
      this.Variables.z = af*z;
      pid = 0;
      for fid=1:Nf
        path_list = slice.FlowTable{fid, 'Paths'};
        for j = 1:path_list.Width
          pid = pid + 1;
          path_list{j}.local_id = pid;
          path_list{j}.MultiplyBandwidth(af);
        end
      end
      
      this.Variables.r = this.getFlowRate();
      
      if nargout >= 2
        vars = [this.Variables.x; this.Variables.z; this.Variables.r];
      end
      % 			warning('[%s] incomplete!', calledby);
      tf = true;
    end
		
		%% setPathBandwidth
		% Override <SliceOptimizer.setPathBandiwdth>
		function setPathBandwidth(this, x)
		end

		%%
		% Override <Slice.priceOptimalFlowRate>.
		% We introduce the quadratic penalty item of the dual-ADMM method.
		% The variables are the flow-edge variables  and the flow-node
		% variables.
		%
		function [gamma, fval, output] = priceOptimalFlowRate0(this, x0, lambda, q, options)
			global DEBUG INFO; %#ok<NUSED>
			slice = this.hs;
			options = Dictionary(options);
			if strcmpi(this.options.Form, 'compact')
				options.bCompact = true;
			else
				options.bCompact = false;
			end
			Nf = slice.NumberFlows;
			Nl = slice.NumberLinks;
			Nsn = slice.NumberServiceNodes;
			num_var_edge = this.num_vars(1);
			num_var_node = this.num_vars(2);
			% num_var_flow = this.num_vars(3);
			
			%%
			if strcmpi(this.options.OptimizationTool, 'matlab')
				minopts = optimoptions('fmincon');
				minopts.Algorithm = 'interior-point';
				%minopts.Algorithm = 'sqp';
				minopts.SpecifyObjectiveGradient = true;
				minopts.Display = 'off';   %'notify-detailed'; %'notify';
				minopts.OptimalityTolerance = 1e-5;
				minopts.ConstraintTolerance = 1e-3;
				minopts.MaxIterations = 100;
				% minopts.SubproblemAlgorithm = 'cg'; % only take cg steps.
				%minopts.CheckGradients = true;
				%minopts.FiniteDifferenceType = 'central';
				if nargin >= 2 && ~isempty(x0)
					var0 = x0;
				elseif isempty(this.x0)
					% set the initial solution.
					%var0 = rand(this.num_flow_vars,1);
					var0 = zeros(num_vars,1);
				else
					%var0 = this.x0;
					% var0 = rand(this.num_flow_vars,1);
					var0 = zeros(num_vars,1);
				end
				var0 = var0(this.I_active_variables);
				lbs = sparse(length(var0),1);
				if strcmpi(minopts.Algorithm, 'interior-point')
					% The method could be overload by subclasses, so invoking it via class name
					% should be safe.
					minopts.HessianFcn = ...
						@(x,lbd)NormalSliceOptimizer.fcnHessian(x, lbd, this, lambda, q, options);
				end
				warning('off')
				[xs, fval, exitflag, foutput] = ...
					fmincon(@(x)NormalDynamicSliceOptimizer.fcnProfitPrimal(x, this, lambda, q, options), var0, ... % avoid overloading issue
					this.problem.A, this.problem.b, this.problem.Aeq, this.problem.beq, ...
					lbs, [], [], minopts);
				warning('on')
				this.interpretExitflag(exitflag, foutput);
			elseif strcmpi(this.options.OptimizationTool, 'cvx')
				[xs, fval, exitflag, foutput] = cvx_method(lambda, q);
				if exitflag<1
					error('error: CVX failed to solve the problem.');
				end
			else
				error('error: Un-recognized optimization tool.');
			end
			%% CVX
			function [xs, fval, exitflag, output] = cvx_method(lambda, q)
				w = this.hs.Weight;
				p = [this.prices.Link; this.prices.Node];
				plt = options.InterSlicePenalty;
				c = [this.capacities.Link; this.capacities.Node];
				Al = this.As_load(:, this.problem.I_active_xz);
				cvx_problem = []; cvx_optbnd = []; %#ok<NASGU>
				cvx_optval = []; cvx_status = [];
				cvx_slvitr = []; cvx_slvtol = []; cvx_cputime = []; %#ok<NASGU>
				if options.bCompact
					c_nx = this.num_vars(1);
					c_nz = this.num_vars(2);
					c_nr = this.num_vars(3);
				else
					c_nx = this.num_vars(1);
					c_nz = this.num_vars(2);
					c_nr = this.num_vars(3);
				end
				c_nxz = c_nx + c_nz;
				%{
					cvx_begin
					variable xs(n,1);
					%maximize( w*sum(log(1+xs((nxz+1):n)))- p'*Al*xs(1:nxz)-plt/2*norm(max(0,1/plt*(Al*xs(1:nxz)-c+q)))^2 )
					maximize( w*sum(log(1+xs((nxz+1):n)))- p'*Al*xs(1:nxz)-plt/2*sum(pow_pos(lambda+1/plt*(Al*xs(1:nxz)-c+q),2)) )
					subject to
					this.problem.A * xs <= this.problem.b;  %#ok<VUNUS>
					this.problem.Aeq * xs == this.problem.beq;  %#ok<EQEFF>
					xs >= lbs; %#ok<VUNUS>
					cvx_end
				%}
				c_x = []; c_z = []; c_r = [];
				cvx_begin
				variable c_x(c_nx,1);
				variable c_z(c_nz,1);
				variable c_r(c_nr,1);
				minimize( -w*sum(log(1+c_r))+p'*Al*[c_x;c_z]+plt/2*sum(pow_pos(lambda+1/plt*(Al*[c_x;c_z]-c+q),2)) )
				subject to
				this.problem.A(:,1:c_nxz) * [c_x;c_z] <= 0;  %#ok<VUNUS>
				this.problem.Aeq(:,[1:c_nx, c_nxz+(1:c_nr)]) * [c_x;c_r] == 0;  %#ok<EQEFF>
				[c_x;c_z;c_r] >= 0; %#ok<VUNUS>
				cvx_end
				fval = cvx_optval;
				xs = [c_x;c_z;c_r];
				switch cvx_status
					case 'Solved'
						exitflag = 1;
					case 'Inaccurate/Solved'
						exitflag = 2;
					case 'Failed'
						exitflag = -1;
					otherwise
						exitflag = 0;
				end
				output.iterations = cvx_slvitr;
			end
			
			%% output solution
			% assert(this.checkFeasible())
			x = zeros(sum(this.num_vars), 1);
			x(this.I_active_variables) = xs;
			output.temp_vars.x = x(1:num_var_edge);
			output.temp_vars.z = x(num_var_edge+(1:num_var_node));
			output.x0 = x;
			num_varxz = length(xs) - Nf;
			gamma = zeros(slice.Parent.NumberDataCenters+slice.Parent.NumberLinks,1);
			idx_gamma = [slice.Links.PhysicalLink; slice.Parent.NumberLinks+slice.getDCPI];
			gamma(idx_gamma) = ...
				max(0, NormalSliceOptimizer.fcnPenalty(xs(1:num_varxz), this, lambda, q, options.InterSlicePenalty));
			
			% 			output.load.Link = this.As_load(1:Nve,:) * xs(1:num_varxz);
			% 			output.load.Node = this.As_load(Nve+(1:Nd),:) * xs(1:num_varxz);
			output.loads = zeros(slice.Parent.NumberLinks+slice.Parent.NumberDataCenters, 1);
			output.loads(idx_gamma) = [this.As_load(1:Nl,this.problem.I_active_xz) * xs(1:num_varxz);
				this.As_load(Nl+(1:Nsn),this.problem.I_active_xz) * xs(1:num_varxz)];
			output.net_profit = -fval;
			output.iterations = foutput.iterations;
			if strcmpi(this.options.OptimizationTool, 'cvx')
				output.funcCount = foutput.iterations;
			else
				output.funcCount = foutput.funcCount;
			end
		end
		
		function [loads, fval, output] = priceOptimalFlowRateDD(this, x0, lambda, options)
			global DEBUG INFO; %#ok<NUSED>
			options = Dictionary(options);
			slice = this.hs;
			if strcmpi(this.options.Form, 'compact')
				options.bCompact = true;
			else
				options.bCompact = false;
			end
			Nf = slice.NumberFlows;
			Nve = slice.NumberLinks;
			Nd = slice.NumberServiceNodes;
			num_var_edge = this.num_vars(1);
			num_var_node = this.num_vars(2);
			
			%%
			if strcmpi(this.options.OptimizationTool, 'matlab')
				minopts = optimoptions('fmincon');
				minopts.Algorithm = 'interior-point';
				%minopts.Algorithm = 'sqp';
				minopts.SpecifyObjectiveGradient = true;
				minopts.Display = 'off';   %'notify-detailed'; %'notify';
				%minopts.OptimalityTolerance = 1e-5;
				%minopts.ConstraintTolerance = 1e-4;
				minopts.MaxIterations = 100;
				% minopts.SubproblemAlgorithm = 'cg'; % only take cg steps.
				%minopts.CheckGradients = true;
				%minopts.FiniteDifferenceType = 'central';
				if nargin >= 2 && ~isempty(x0)
					var0 = x0;
				elseif isempty(this.x0)
					% set the initial solution.
					%var0 = rand(this.num_flow_vars,1);
					var0 = zeros(num_vars,1);
				else
					%var0 = this.x0;
					% var0 = rand(this.num_flow_vars,1);
					var0 = zeros(num_vars,1);
				end
				var0 = var0(this.I_active_variables);
				lbs = sparse(length(var0),1);
				if strcmpi(minopts.Algorithm, 'interior-point')
					minopts.HessianFcn = ...
						@(x,lbd)NormalSliceOptimizer.fcnHessianDD(x, lbd, this, lambda, options);
				end
				warning('off')
				[xs, fval, exitflag, foutput] = ...
					fmincon(@(x)NormalSliceOptimizer.fcnProfitPrimalDD(x, this, lambda, options), var0, ...
					this.problem.A, this.problem.b, this.problem.Aeq, this.problem.beq, ...
					lbs, [], [], minopts);
				warning('on')
				this.interpretExitflag(exitflag, foutput);
			else
				error('error: Un-recognized optimization tool.');
			end
			
			%% output solution
			% assert(this.checkFeasible())
			x = zeros(num_vars, 1);
			x(this.I_active_variables) = xs;
			output.temp_vars.x = x(1:num_var_edge);
			output.temp_vars.z = x(num_var_edge+(1:num_var_node));
			output.temp_vars.r = x(num_var_edge+num_var_node+(1:Nf));
			output.x0 = x;
			num_varxz = length(xs) - Nf;
			
			loads = zeros(slice.Parent.NumberLinks+slice.Parent.NumberDataCenters, 1);
			idx = [slice.Links.PhysicalLink; slice.Parent.NumberLinks+slice.getDCPI];
			loads(idx) = [this.As_load(1:Nve,this.problem.I_active_xz) * xs(1:num_varxz);
				this.As_load(Nve+(1:Nd),this.problem.I_active_xz) * xs(1:num_varxz)];
			output.net_profit = -fval;
			output.iterations = foutput.iterations;
			output.funcCount = foutput.funcCount;
		end
		
		function [gamma, fval, x_prox, output] = priceOptimalFlowRatePP(this, x0, lambda, q, xp, options)
			global DEBUG INFO; %#ok<NUSED>
			options = Dictionary(options);
			if strcmpi(this.options.Form, 'compact')
				options.bCompact = true;
			else
				options.bCompact = false;
			end
			if options.InterSlicePenalty == 0
				options.InterSlicePenalty = this.hs.Weight;
			end
			slice = this.hs;
			Nf = slice.NumberFlows;
			Nve = slice.NumberLinks;
			Nd = slice.NumberServiceNodes;
			num_var_edge = this.num_vars(1);
			num_var_node = this.num_vars(2);
			
			%%
			xp = xp(this.I_active_variables);
			if strcmpi(this.options.OptimizationTool, 'matlab')
				minopts = optimoptions('fmincon');
				minopts.Algorithm = 'interior-point';
				%minopts.Algorithm = 'sqp';
				minopts.SpecifyObjectiveGradient = true;
				minopts.Display = 'off';   %'notify-detailed'; %'notify';
				minopts.OptimalityTolerance = 1e-5;
				minopts.ConstraintTolerance = 1e-3;
				minopts.MaxIterations = 100;
				% minopts.SubproblemAlgorithm = 'cg'; % only take cg steps.
				%minopts.CheckGradients = true;
				%minopts.FiniteDifferenceType = 'central';
				if nargin >= 2 && ~isempty(x0)
					var0 = x0;
				else
					%var0 = this.x0;
					% var0 = rand(this.num_flow_vars,1);
					%var0 = zeros(this.num_flow_vars,1);
					var0 = xp;
				end
				%var0 = var0(this.I_active_variables);
				lbs = sparse(length(var0),1);
				options.TrueValue = false;
				if strcmpi(minopts.Algorithm, 'interior-point')
					minopts.HessianFcn = ...
						@(x,lbd)NormalSliceOptimizer.fcnHessianPP(x, lbd, this, lambda, q, options);
				end
				warning('off')
				[xs, fval, exitflag, foutput] = ...
					fmincon(@(x)NormalSliceOptimizer.fcnProfitPrimalPP(x, this, lambda, q, xp, options), var0, ...
					this.problem.A, this.problem.b, this.problem.Aeq, this.problem.beq, ...
					lbs, [], [], minopts);
				warning('on')
				this.interpretExitflag(exitflag, foutput);
			elseif strcmpi(this.options.OptimizationTool, 'cvx')
				[xs, fval, exitflag, foutput] = cvx_method(lambda, q, xp);
				if exitflag<1
					error('error: CVX failed to solve the problem.');
				end
			else
				error('error: Un-recognized optimization tool.');
			end
			
			%% CVX
			function [xs, fval, exitflag, output] = cvx_method(lambda, q, xp)
				w = this.hs.Weight;
				p = [this.prices.Link; this.prices.Node];
				r = options.InterSlicePenalty;
				Ns = options.NumberSlices;
				c = [this.capacities.Link; this.capacities.Node];
				Al = this.As_load(:,this.problem.I_active_xz);
				cvx_problem = []; cvx_optbnd = []; %#ok<NASGU>
				cvx_optval = []; cvx_status = [];
				cvx_slvitr = []; cvx_slvtol = []; cvx_cputime = []; %#ok<NASGU>
				if options.bCompact
					c_nx = this.num_vars(1);
					c_nz = this.num_vars(2);
					c_nr = this.num_vars(3);
				else
					c_nx = this.num_vars(1);
					c_nz = this.num_vars(2);
					c_nr = this.num_vars(3);
				end
				c_nxz = c_nx + c_nz;
				c_x = []; c_z = []; c_r = [];
				cvx_begin
				variable c_x(c_nx,1);
				variable c_z(c_nz,1);
				variable c_r(c_nr,1);
				minimize( -w*sum(log(1+c_r))+p'*Al*[c_x;c_z]+ 1/(2*r*Ns^2)*sum(pow_abs([c_x;c_z;c_r]-xp,2)) + r/2*sum(pow_pos(lambda+1/r*(Al*[c_x;c_z]-c+q),2)) )
				subject to
				this.problem.A(:,1:c_nxz) * [c_x;c_z] <= 0;  %#ok<VUNUS>
				this.problem.Aeq(:,[1:c_nx, c_nxz+(1:c_nr)]) * [c_x;c_r] == 0;  %#ok<EQEFF>
				[c_x;c_z;c_r] >= 0; %#ok<VUNUS>
				cvx_end
				fval = cvx_optval;
				xs = [c_x;c_z;c_r];
				switch cvx_status
					case 'Solved'
						exitflag = 1;
					case 'Inaccurate/Solved'
						exitflag = 2;
					case 'Failed'
						exitflag = -1;
					otherwise
						exitflag = 0;
				end
				output.iterations = cvx_slvitr;
			end
			
			%% output solution
			% assert(this.checkFeasible())
			x = zeros(num_vars, 1);
			xs(xs<10^-4*norm(xs)) = 0;
			x(this.I_active_variables) = xs;
			output.temp_vars.x = x(1:num_var_edge);
			output.temp_vars.z = x(num_var_edge+(1:num_var_node));
			output.temp_vars.r = x(num_var_edge+num_var_node+(1:Nf));
			output.x0 = x;
			options.TrueValue = true;
			fval = fval + NormalSliceOptimizer.fcnProfitPrimalPP(xs, this, lambda, q, xp, options);
			num_varxz = length(xs) - Nf;
			x_prox = x;
			
			gamma = zeros(slice.Parent.NumberDataCenters+slice.Parent.NumberLinks,1);
			idx = [slice.Links.PhysicalLink; slice.Parent.NumberLinks+slice.getDCPI];
			gamma(idx) = ...
				max(0, NormalSliceOptimizer.fcnPenalty(xs(1:num_varxz), this, lambda, q, options.InterSlicePenalty));
			output.loads = zeros(slice.Parent.NumberLinks+slice.Parent.NumberDataCenters, 1);
			output.loads(idx) = [this.As_load(1:Nve,this.problem.I_active_xz) * xs(1:num_varxz);
				this.As_load(Nve+(1:Nd),this.problem.I_active_xz) * xs(1:num_varxz)];
			output.net_profit = -fval;
			output.iterations = foutput.iterations;
			if strcmpi(this.options.OptimizationTool, 'cvx')
				output.funcCount = foutput.iterations;
			else
				output.funcCount = foutput.funcCount;
			end
		end
		
	end
	
	methods (Access = protected)
		% 		function newobj = copyElement(this)
		% 			if this.isShallowCopyable
		% 				newobj = copyElement@SliceOptimizer(this);
		% 			else
		% 				newobj = this;
		% 			end
		% 			newobj.augraph = this.augraph.copy();
		% 		end
		
		% this is valid for later procedure to identify active variables.
		function I_active_vars = updateActiveVariables(this)
			prbm = this.problem;
			num_varz = this.num_vars(2);
			num_varxz = sum(this.num_vars(1:2));
			num_vars_prim = sum(this.num_vars(1:3));
			num_lcon_proc = size(this.As_proc,1);
			A = [prbm.Aeq(:,1:num_vars_prim); prbm.A(1:num_lcon_proc, 1:num_vars_prim)];
			this.I_active_variables = sparse(sum(abs(A), 1)~=0);
			prbm.I_fake_vars = [this.I_fake_edgevars; sparse(num_varz,1)];
			prbm.I_active_xz = this.I_active_variables(1:num_varxz);
			if nargout >= 1
				I_active_vars = this.I_active_variables;
			end
		end
		
		%% Optimize
		% This method is separated from <optimalFlowRate> to provide override capability for
		% subclasses.
		function [x, fval] = optimize(this, options)			
			prbm = this.problem;
			[xs, fval, exitflag, output] = ...
				fmincon(@(x)SimpleSliceOptimizer.fcnProfit(x, this, options), ...
				prbm.x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, lb, ub, [], prbm.minopts);
			SimpleSliceOptimizer.interpretExitflag(exitflag, output);
			if options.bCompact
				x = zeros(sum(this.num_vars), 1);
				x(this.I_active_variables) = xs;
			else
				x = xs;
			end
			constr_tol = getstructfields(this.problem.options, 'ConstraintTolerance', ...
				'default-ignore', {this.options.ConstraintTolerance});
			assert(this.checkFeasible(x, constr_tol), 'error: infeasible solution.');
		end
			
		% Call this function after flow is added to slice's FlowTable or before flow is
		% removed from slice's flow table.
		function updateFlowSectionTable(this, fidx)
			slice = this.hs;
			if nargin <= 1
				fidx = 1:slice.NumberFlows;
			end
			num_segs = slice.NumberVNFs + 1;
			Nf = numel(fidx);
			Nvf = Nf * num_segs;
			flow_idxs = repelem(1:Nf, num_segs);
			sources = zeros(Nvf, 1);
			destins = zeros(Nvf, 1);
			sources(1:num_segs:(Nf-1)*num_segs+1) = slice.FlowTable{fidx, 'Source'};
			destins(num_segs:num_segs:Nvf) = slice.FlowTable{fidx, 'Target'};
			this.flow_section_table = [this.flow_section_table;
				table(flow_idxs(:), sources, destins, ...
				'VariableNames', {'FlowIndex', 'Source', 'Destination'})];
		end
	end % end protected methods
	
	methods (Access = protected, Static)
		%% Objective function
		% See also <>
		% NOTE: 'vars' might be compressed.
		function [profit, grad] = fcnProfitPrimal(vars, this, lambda, q, options)
			plt = options.InterSlicePenalty;
			base_price = [this.prices.Link; this.prices.Node];

			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			Nf = pardata.NumberFlows;
			Nv = pardata.NumberVNFs;
			Nvf = this.NumberFlowSections;
			prbm = this.problem;
			
			num_varxz = length(vars) - Nf;
			flow_rate = vars((num_varxz+1):end);
			load = this.As_load(:,prbm.I_active_xz) * vars(1:num_varxz);
			
			profit = -pardata.Weight*sum(fcnUtility(flow_rate)) ...
				+ dot(base_price, load);
			%% TODO
			penalty_item = NormalSliceOptimizer.fcnPenalty(vars(1:num_varxz), this, lambda, q, plt);
			profit = profit	+ plt/2*(norm(max(0,penalty_item)))^2;
			if options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				grad = sparse(sum(this.num_vars), 1);
				grad(~this.I_fake_edgevars) = repmat(this.prices.Link,Nvf,1);
				grad(this.num_vars(1)+(1:this.num_vars(2))) = repmat(this.prices.Node,Nf*Nv,1);
				grad = grad(this.I_active_variables);
				%% update gradient on (x,z)
				idx = penalty_item<=0;
				A = this.As_load(:,prbm.I_active_xz);
				A(idx, :) = 0;
				grad(1:num_varxz) = grad(1:num_varxz) + (A')*max(0, penalty_item);
				%% update gradient on r
				grad((num_varxz+1):end) = -pardata.Weight./(1+flow_rate);
			end
		end
		
		%%
		% Override <Slice.fcnHessian>.
		function h = fcnHessian(vars, ~, this, lambda, q, options)
			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			num_vars = length(vars);
			num_varxz = num_vars - pardata.NumberFlows;
			flow_rate = vars((num_varxz+1):end);
			
			h = sparse(num_vars, num_vars);
			h((num_varxz+1):end, (num_varxz+1):end) = spdiag(pardata.Weight./(1+flow_rate).^2);
			penalty_item = NormalSliceOptimizer.fcnPenalty(vars(1:num_varxz), this, lambda, q, ...
				options.InterSlicePenalty);
			idx = penalty_item<=0;
			A = this.As_load(:,this.problem.I_active_xz);
			A(idx, :) = 0;
			h(1:num_varxz, 1:num_varxz) = h(1:num_varxz,1:num_varxz) + (A')*A;
		end
		
		%% Penalty function
		function p = fcnPenalty(vars, this, lambda, q, factor)
			load = this.As_load(:,this.problem.I_active_xz) * vars;
			capacity = [this.capacities.Link; this.capacities.Node];
			p = lambda+1/factor*(load-capacity+q);
		end
		
		
		%% Functions for Dual Decomposition
		function [profit, grad] = fcnProfitPrimalDD(vars, this, lambda, options)
			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			Nf = pardata.NumberFlows;
			Nvnf = pardata.NumberVNFs;
			Nvf = this.NumberFlowSections;
			prbm = this.problem;
			nx = prbm.num_vars(1);
			nz = prbm.num_vars(2);
			nxz = nx + nz;
			
			flow_rate = vars((nxz+1):end);
			load = this.As_load(:,prbm.I_active_xz) * vars(1:nxz);
			
			profit = -pardata.Weight*sum(fcnUtility(flow_rate)) ...
				+ dot([this.prices.Link; this.prices.Node]+lambda, load);
			profit = profit + (0.001/2)*sum(vars(1:nxz).^2);
			if options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				Nl = length(this.prices.Link);
				grad = sparse(sum(this.num_vars), 1);
				grad(~this.I_fake_edgevars) = repmat(this.prices.Link+lambda(1:Nl),Nvf,1);
				grad(this.num_vars(1)+(1:this.num_vars(2))) = ...
					repmat(this.prices.Node+lambda(Nl+1:end),Nf*Nvnf,1);
				grad = grad(this.I_active_variables);
				%%
				grad(1:nxz) = grad(1:nxz) + 0.001*vars(1:nxz);
				%% update gradient on r
				grad((nxz+1):end) = -pardata.Weight./(1+flow_rate);
			end
		end
		
		%%
		function h = fcnHessianDD(vars, ~, this, ~, options) %#ok<INUSD>
			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			num_vars = length(vars);
			num_varxz = num_vars - slice.NumberFlows;
			flow_rate = vars((num_varxz+1):end);
			
			h = sparse(num_vars, num_vars);
			h(1:num_varxz, 1:num_varxz) = speye(num_varxz)*0.001;
			h((num_varxz+1):end, (num_varxz+1):end) = spdiag(pardata.Weight./(1+flow_rate).^2);
		end
		
		function [profit, grad] = fcnProfitPrimalPP(vars, this, lambda, q, xp, options)
			r = options.InterSlicePenalty;
			base_price = [this.prices.Link; this.prices.Node];
			
			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			Nf = pardata.NumberFlows;
			Nvnf = pardata.NumberVNFs;
			Nvf = this.NumberFlowSections;
			Ns = options.NumberSlices;
			
			num_varxz = length(vars) - Nf;
			flow_rate = vars((num_varxz+1):end);
			load = this.As_load(:,this.problem.I_active_xz) * vars(1:num_varxz);
			
			profit = -slice.Weight*sum(fcnUtility(flow_rate)) + dot(base_price, load);
			if options.TrueValue
				return;
			end
			profit = profit + 1/(2*r*Ns^2)*sum((vars-xp).^2);
			%% TODO
			penalty_item = NormalSliceOptimizer.fcnPenalty(vars(1:num_varxz), this, lambda, q, r);
			profit = profit	+ r/2*sum((max(0,penalty_item)).^2);
			if options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				grad = sparse(sum(this.num_vars), 1);
				grad(~this.I_fake_edgevars) = repmat(this.prices.Link,Nvf,1);
				grad(this.num_vars(1)+(1:this.num_vars(2))) = repmat(this.prices.Node,Nf*Nvnf,1);
				grad = grad(this.I_active_variables);
				%% update gradient on r
				grad((num_varxz+1):end) = -pardata.Weight./(1+flow_rate);
				%% update gradient on (x,z)
				idx = penalty_item<=0;
				A = this.As_load(:,this.problem.I_active_xz);
				A(idx, :) = 0;
				grad(1:num_varxz) = grad(1:num_varxz) + (A')*max(0, penalty_item);
				%% gradient on proximal penalty
				grad = grad + 1/(r*Ns^2)*(vars-xp);
			end
		end
		
		function h = fcnHessianPP(vars, ~, this, lambda, q, options)
			r = options.InterSlicePenalty;
			num_vars = length(vars);
			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			num_varxz = num_vars - pardata.NumberFlows;
			flow_rate = vars((num_varxz+1):end);
			Ns = options.NumberSlices;
			
			h = 1/(r*Ns^2)*speye(num_vars);
			h((num_varxz+1):end, (num_varxz+1):end) = h((num_varxz+1):end, (num_varxz+1):end) ...
				+ spdiag(pardata.Weight./(1+flow_rate).^2);
			penalty_item = NormalSliceOptimizer.fcnPenalty(vars(1:num_varxz), this, lambda, q, ...
				r);
			idx = penalty_item<=0;
			A = this.As_load(:,this.problem.I_active_xz);
			A(idx, :) = 0;
			h(1:num_varxz, 1:num_varxz) = h(1:num_varxz,1:num_varxz) + (A')*A;
		end
		
		function [profit, grad] = fcnProfitPrimalCC(vars, this, options)
			if isfield(options, 'unit')
				unit = options.unit;
			else
				unit = 1;
			end
			% |slice| acts as instance of <Slice> in normal mode or |pardata| in parallel mode.
			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			Nf = pardata.NumberFlows;
			Nvnf = pardata.NumberVNFs;
			Nl = pardata.NumberLinks;
			Nvf = this.NumberFlowSections;
			base_price = [this.prices.Link; this.prices.Node];
			prbm = this.problem;
			
			num_varxz = length(vars) - Nf;
			flow_rate = vars((num_varxz+1):end);
			load = this.As_load(:, prbm.I_active_xz) * vars(1:num_varxz);
			
			profit = -pardata.Weight*sum(fcnUtility(unit*flow_rate));
			switch options.PricingPolicy
				case {'quadratic-price', 'quadratic'}
					link_load = load(1:Nl);
					node_load= load((Nl+1):end);
					[link_payment,link_price_grad] = this.fcnLinkPricing(this.prices.Link, link_load, unit);
					[node_payment,node_price_grad] = this.fcnNodePricing(this.prices.Node, node_load, unit);
					
					profit = profit + link_payment + node_payment;
				case 'linear'
					profit = profit + dot(base_price, load*unit);
				otherwise
					error('%s: invalid pricing policy', calledby);
			end
			if options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				grad = sparse(sum(this.num_vars), 1);		% compress it later
				switch options.PricingPolicy
					case {'quadratic-price', 'quadratic'}
						grad(~this.I_fake_edgevars) = repmat(link_price_grad,Nvf,1);
						grad(this.num_vars(1)+(1:this.num_vars(2))) = repmat(node_price_grad,Nf*Nvnf,1);
					case 'linear'
						grad(~this.I_fake_edgevars) = repmat(this.prices.Link,Nvf,1);
						grad(this.num_vars(1)+(1:this.num_vars(2))) = repmat(this.prices.Node,Nf*Nvnf,1);
				end
				grad = grad(this.I_active_variables);
				%% update gradient on r
				grad((num_varxz+1):end) = -pardata.Weight./(1/unit+flow_rate);
			end
		end
		
		function h = fcnHessianCC(vars, ~, this, options)
			if nargin >= 4 && isfield(options, 'unit')
				unit = options.unit;
			else
				unit = 1;
			end
			num_vars = length(vars);
			num_varxz = sum(this.problem.num_vars(1:2));
			if isempty(this.hs)
				pardata = this.pardata;
			else
				pardata = this.hs;
			end
			Nf = pardata.NumberFlows;
			Nvnf = pardata.NumberVNFs;
			Nl = pardata.NumberLinks;
			flow_rate = vars((num_varxz+1):end);
			h = sparse(num_vars, num_vars);
			h((num_varxz+1):end, (num_varxz+1):end) = spdiag(pardata.Weight./(1/unit+flow_rate).^2);
			if nargin >= 4 && isfield(options, 'PricingPolicy')
				switch options.PricingPolicy
					case {'quadratic-price', 'quadratic'}
						load = this.As_load(:, this.problem.I_active_xz) * vars(1:num_varxz);
						link_load = load(1:Nl);
						node_load = load((Nl+1):end);
						[~,~,de] = this.fcnLinkPricing(this.prices.Link, link_load, unit); % equal to <getLinkCapacity>
						[~,~,dn] = this.fcnNodePricing(this.prices.Node, node_load, unit);
						numvars_edgenodes = sum(this.num_vars(1:2));
						blocks = sparse(numvars_edgenodes, numvars_edgenodes);
						prbm = this.problem;
						blocks(~prbm.I_fake_vars, ~prbm.I_fake_vars) = ...
							blkdiag(block_diag(diag(de), Nf*(Nvnf+1)), block_diag(diag(dn),Nf*Nvnf));
						h(1:num_varxz,1:num_varxz) = blocks(prbm.I_active_xz, prbm.I_active_xz);
				end
			end
		end
		
		function [profit, grad] = fcnProfit(vars, this, options)       	% defined in SliceOptimizer
			warning('call fcnProfitPrimalCC directly!');
			[profit, grad] = NormalSliceOptimizer.fcnProfitPrimalCC(vars, this, options);
		end
		
		%% Net social welfare of a slice
		function [profit, grad] = fcnSocialWelfare(vars, this, options)
			% 			if options.bCompact
			% 				full_vars = sparse(options.num_orig_vars,1);
			% 				full_vars(this.I_active_variables) = vars;
			% 				vars = full_vars;
			% 			end
			prbm = this.problem;
			var_r_index = sum(prbm.num_vars(1:2))+(1:sum(prbm.num_vars(3)));
			flow_rate = vars(var_r_index);
			%% Calculate the cost of network and slices
			% When calculate the cost of the network as a single slice, we can use both
			% _Slice.getSliceCost_ and _CloudNetwork.getNetworkCost_. Otherwise, use
			% _Slice.getSliceCost_ to calculate each slice's cost. *For simplicity, we  only use
			% _getSliceCost_*.
			%
			% Called by _optimalFlowRate_ as objective function;
			%%%
			%     profit = -sum(weight.*fcnUtility(flow_rate)) + ...
			%         S.Parent.totalCost(load);
			if isempty(this.pardata.weight)        % for network as a single slice.    TODO: remove this block
				weight = this.pardata.FlowWeight;
			else                        % for network as multiple slice.
				weight = this.pardata.Weight*ones(this.pardata.NumberFlows, 1);
			end
			num_varxz = sum(prbm.num_vars(1:2));
			loads = this.As_load(:, this.I_active_xz)*vars(1:num_varxz);
			load.Node = loads(Nl+(1:Nsn));
			load.Link = loads(1:Nl);
			profit = -sum(weight.*fcnUtility(flow_rate)) + slice.getResourceCost(load);
			
			if options.bFinal
				profit = -profit;
			else
				link_price = this.pardata.LinkCost;
				node_price = this.pardata.NodeCost;
				grad = sparse(sum(this.num_vars), 1);		% compress it later
				grad(~this.I_fake_edgevars) = repmat(link_price,Nvf,1);
				grad(this.num_vars(1)+(1:this.num_vars(2))) = repmat(node_price,Nf*Nvnf,1);
				grad = grad(this.I_active_variables);
				%% update gradient on r
				grad((num_varxz+1):end) = -weight./(1+flow_rate);
			end
			
		end
		
	end
	
	
end % end of class
%{
				Aeq = sparse(this.NumberAugmentedNodes*(this.NumberVNFs+1)*Nf, this.num_flow_vars);
				beq = sparse(size(Aeq,1),1);
				Aeq(this.I_active_rows_eq, this.I_active_variables) = this.problem.Aeq;
				beq(this.I_active_rows_eq) = this.problem.beq;
				A = sparse(this.NumberDataCenters*this.NumberVNFs*Nf, this.num_flow_vars);
				b = sparse(size(A,1),1);
				A(this.I_active_rows, this.I_active_variables) = this.problem.A;
				b(this.I_active_rows) = this.problem.b;
				options_ = options;
				options_.bCompact = false;
				var0_ = sparse(this.num_flow_vars,1);
				lbs_ = sparse(this.num_flow_vars,1);
				As_load_temp = this.As_load;
				this.As_load = sparse(size(this.As_load,1), num_vars_edge+num_vars_node);
				this.As_load(:,this.I_active_variables(1:num_vars_edge+num_vars_node)) = As_load_temp;
				[xs_, fval_, exitflag_, foutput_] = ...
					fmincon(@(x)NormalSlizeOptimizer.fcnProfitPrimal(x, this, lambda, q, options_), var0_, ...
					A, b, Aeq, beq, lbs_, [], [], minopts);
				this.As_load = As_load_temp;
%}