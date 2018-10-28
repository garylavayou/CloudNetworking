classdef NormalSlice < Slice
	%UNTITLED Summary of this class goes here
	%   Detailed explanation goes here
	
	properties (Access = protected)
		%% Augmented Topology
		augraph DirectedGraph;
		
		%% Virtual Flow Table
		% We need to segment each flow into multiple segements
		%% TODO: remove this variable, no specific usage.
		flow_section_table;
	end
	

	
	properties (Dependent)
		NumberAugmentedNodes;
		NumberAugmentedLinks;
		NumberFlowSections;
	end
	
	properties(Dependent, GetAccess={?FlowEdgeSlice, ?SubstrateNetwork})
		%% Number of variables in the problem
		% The flows are break into segments according to the number of VNFs in
		% the service chain.
		num_flow_vars;	% (original) total number variables
		num_vars_edge;	% number of edge variables
		num_vars_node;	% number of node variables
		num_vars_rate;  % number of flows
	end
	
	methods
		function this = NormalSlice(slice_data)
			if nargin == 0
				slice_data = {};
			end
			this@Slice(slice_data);
			if nargin == 0
				return;
			end
			
			%% Build virtual flow table
			% Fields include: flow index, source, destination, etc.
			if this.NumberFlows > 0
				num_segs = this.NumberVNFs + 1;
				flow_idxs = repelem((1:this.NumberFlows), num_segs);
				num_vflows = length(flow_idxs);
				sources = zeros(num_vflows, 1);
				destins = zeros(num_vflows, 1);
				sources(1:num_segs:(this.NumberFlows-1)*num_segs+1) = this.FlowTable.Source;
				destins(num_segs:num_segs:num_segs*this.NumberFlows) = this.FlowTable.Target;
				this.flow_section_table = table(flow_idxs(:), sources, destins, ...
					'VariableNames', {'FlowIndex', 'Source', 'Destination'});
			end
			%%
			
		end
		
		function n = get.NumberAugmentedNodes(this)
			n = this.augraph.NumberNodes;
		end
		
		function n = get.NumberAugmentedLinks(this)
			n = this.augraph.NumberEdges;
		end
		
		function n = get.NumberFlowSections(this)
			n = this.NumberFlows * (this.NumberVNFs+1);
		end
		
		function n = get.num_flow_vars(this)
			n = this.num_vars_edge + this.num_vars_node + this.num_vars_rate;
		end
		
		function n = get.num_vars_edge(this)
			n = height(this.flow_section_table) * this.NumberAugmentedLinks;
		end
		
		function n = get.num_vars_node(this)
			n = this.NumberVNFs * this.NumberFlows * this.NumberDataCenters;
		end
		
		function n = get.num_vars_rate(this)
			n = this.NumberFlows;
		end
		
		function initializeState(this)
			initializeState@Slice(this);
			%% Build the augmented topology
			this.augraph = this.graph.copy;
			this.ServiceNodes{:, 'AugmentedVirtualNode'} = ...
				this.NumberNodes + (1:this.NumberServiceNodes)';
			this.ServiceNodes{:, 'AugmentedNode'} = ...
				this.Parent.DataCenters.AugmentedNode(this.getDCPI);
			aug_heads = [this.ServiceNodes.VirtualNode; ...
				this.ServiceNodes.AugmentedVirtualNode];
			aug_tails = [this.ServiceNodes.AugmentedVirtualNode;...
				this.ServiceNodes.VirtualNode];
			num_aug_links = length(aug_heads);
			props.Weight = eps*ones(num_aug_links,1);
			props.Capacity = inf*ones(num_aug_links,1);
			this.augraph.Update(aug_heads, aug_tails, props);		% TODO: remap
			
			%% Initialize edge-flow & node-flow incidence matrix
			this.I_flow_node = sparse(this.NumberFlows, this.NumberNodes);
			if ~isempty(this.I_flow_path)
				this.I_flow_edge = double(this.I_flow_path)*double(transpose(this.I_edge_path));
				for i = 1:this.NumberFlows
					pl = this.FlowTable{i, 'Paths'};	% PathList
					for j = 1:pl.Width
						p = pl.paths{j};
						this.I_flow_node(i, p.node_list) = true;
					end
				end
			else
				this.I_flow_edge = true(this.NumberFlows, this.NumberVirtualLinks);
				this.I_flow_node = true(this.NumberFlows, this.NumberVirtualNodes);
			end
			this.I_flow_edge_ex = this.I_flow_edge;
			this.I_flow_edge_ex(this.NumberFlows, this.NumberAugmentedLinks) = 0;
			this.I_flow_node_ex = this.I_flow_node;
			this.I_flow_node_ex(this.NumberFlows, this.NumberAugmentedNodes) = 0;
			for k = 1:this.NumberFlows
				for i = (this.NumberLinks+1):this.NumberAugmentedLinks
					[src,dest] = this.augraph.IndexEdge(i);
					if src<=this.NumberNodes && this.I_flow_node(k, src)
						this.I_flow_edge_ex(k, i) = 1;
						this.I_flow_node_ex(k, dest) = 1;
					end
					if dest<=this.NumberNodes && this.I_flow_node(k, dest)
						this.I_flow_edge_ex(k, i) = 1;
						this.I_flow_node_ex(k, src) = 1;
					end
				end
			end
		end
		
		%%
		% Override <Slice.getFlowRate>
		function r = getFlowRate(this)
			r = this.flow_rate;
		end
		
		%%
		% Override <Slice.setPathBandiwdth>
		% Convert the edge-flow variables into path-flow variables.
		function setPathBandwidth(this, x)
			if nargin == 1
				x = this.Variables.x;
			end
			num_segs = this.NumberVNFs+1;
			Nvn = this.NumberNodes;
			Navn = this.NumberAugmentedNodes;
			Nave = this.NumberAugmentedLinks;
			nidx_offset = 0;
			eidx_offset = 0;
			for i = 1:this.NumberFlows
				As = cell(num_segs,1);
				gen_path = cell(0);
				while true
					fx = x(eidx_offset+(1:Nave*num_segs));
					t_nidx_offset = nidx_offset;
					t_eidx_offset = eidx_offset;
					%% Copy the incident matrix for the flow
					for j = 1:num_segs
						As{j} = this.As_flow(t_nidx_offset+(1:Navn), t_eidx_offset+(1:Nave));		% As{j} are equal for all j
						eidx = fx(t_eidx_offset+(1:Nave))<eps;
						As{j}(:,eidx) = 0;
						t_eidx_offset = t_eidx_offset + Nave;
						t_nidx_offset = t_nidx_offset + Navn;
					end
					src = this.FlowTable{i, 'Source'};
					%% build a single path edge-by-edge
					path.nodes = src;
					path.idx_edgevar = [];		% the index for corresponding edge variables.
					for j = 1:num_segs
						while true
							eout = find(As{j}(src,:)==-1);
							if ~isempty(eout)
								[~,idx] = max(fx(eout));
								[n_cur, n_next] = this.augraph.IndexEdge(eout(idx));
								assert(n_cur==src);
								path.nodes = [path.nodes; n_next];
								path.idx_edgevar = [path.idx_edgevar; eout(idx)];
								As{j}(src,eout(idx)) = 0;
								src = n_next;
							else
								path.idx_edgevar = path.idx_edgevar+Nave*(j-1);
								break;
							end
						end
					end
					assert(length(path.nodes)>1);
					idx_aug_nodes = find(path.nodes>Nvn);
					idx_aug_edges = reshape([idx_aug_nodes-1,idx_aug_nodes]', ...
						2*length(idx_aug_nodes),1);
					path.idx_simple_edgevar = path.idx_edgevar;
					path.idx_simple_edgevar(idx_aug_edges) = [];
					min_bandwidth = min(fx(path.idx_simple_edgevar));
					fx(path.idx_simple_edgevar) = fx(path.idx_simple_edgevar) - min_bandwidth;
					gen_path = [gen_path; {Path(path.nodes, min_bandwidth)}]; %#ok<AGROW>
					if isempty(find(fx>=eps,1))
						break;
					end
				end
				this.FlowTable{i, 'Paths'} = PathList(gen_paths);
				eidx_offset = eidx_offset + Nave*num_segs;
				nidx_offset = nidx_offset + Navn*num_segs;
			end
		end
		
		
		function profit = getProfit(slice, options)
			% determine varriables.
			if nargin >= 2 && isfield(options, 'bFinal') && options.bFinal
				vars = [slice.Variables.x; slice.Variables.z; slice.flow_rate];
			else
				vars = [slice.temp_vars.x; slice.temp_vars.z; slice.flow_rate];
			end
			vars = vars(slice.I_active_variables);
			if nargin >= 2
				if isfield(options, 'LinkPrice')
					slice.prices.Link = options.LinkPrice;
				elseif isfield(options, 'bFinal') && options.bFinal
					slice.prices.Link = slice.VirtualLinks.Price;
				end
				if isfield(options, 'NodePrice')
					slice.prices.Node = options.NodePrice;
				elseif isfield(options, 'bFinal') && options.bFinal
					slice.prices.Node= slice.ServiceNodes.Price;
				end
			end
			if nargin < 2 || ~isfield(options, 'PricingPolicy')
				options.PricingPolicy = 'linear';
				warning('%s: <PricingPolicy> not specifed, set to [%s].', ...
					calledby, options.PricingPolicy);
			end
			%%
			% Here, the invoked method must be <Slice.fcnProfit>.
			% Use class name to avoid dynamic loading of class.
			% Subclasses may override this method to define different ways to calculate profits.
			options.bFinal = true;
			profit = FlowEdgeSlice.fcnProfitPrimalCC(vars, slice, options);
			slice.prices.Node = [];
			slice.prices.Link = [];
		end
	end
	
	methods
		%% Parallel
		function [array, problem, indices] = initializeProblem(this, options)
			if nargin < 2 || ~isfield(options, 'bCompact')
				options.bCompact = false;
			end
			%%
			Nave = this.NumberAugmentedLinks;
			Nve = this.NumberVirtualLinks;
			Navn = this.NumberAugmentedNodes;
			Nvn = this.NumberVirtualNodes;
			Nd = this.NumberDataCenters;
			Nf = this.NumberFlows;
			Nvf = height(this.flow_section_table);
			Nv = this.NumberVNFs;
			num_seg = Nv + 1;
			
			%% Flow reservation constraints and processing resource constraints
			%		[As -Ir] * [F; r]  = 0   ==> As*F  = r
			%		[Al -Iz] * [F; z] <= 0   ==> Al*F <= z
			% ===>
			%								 ┌F┐										┌F┐
			%		[As 0 -Ir] * |z| = 0,  [Al -Iz 0] * |z| <= 0
			%								 └r┘										└r┘
			%
			% Several situation should be considered:
			%	1. intermediate forwarding nodes;
			%	2. source forwarding nodes of the first segment;
			% 3. destiniation forwarding nodes of the last segment;
			%	4. data center nodes that potentially process flows.
			%
			A = this.augraph.GetIncidentMatrix;		% There might exist single direction links in the virtual topology.
			%{
			n = 56;
			eidx = find(A(n,:)==-1);
			[s,t] = this.augraph.IndexEdge(eidx);
			fprintf('outgoing links:\n');
			disp([s,t]);
			eidx = find(A(n,:)==1);
			[s,t] = this.augraph.IndexEdge(eidx);
			fprintf('incoming links:\n');
			disp([s,t]);
			%}
			As = block_diag(A, num_seg);
			As_proc_temp = As;			% for processing requirement constraints
			for nidx = (Nvn+1):Navn
				eout = find(A(nidx,:)==-1);
				for j = 1:(num_seg-1)	
					As((j-1)*Navn+nidx, (j-1)*Nave+eout) = 0;
					As((j-1)*Navn+nidx, j*Nave+eout) = -1;		% TODO:  -1 => -α : rate change
					As_proc_temp(nidx+(j-1)*Navn, eout+(j-1)*Nave) = 0;
				end
				%% The last segment does not go into processing nodes
				% The incomming traffic of the last segment at the processing
				% nodes should be limited to 0. So the corresponding variables
				% could be eliminated (We set the coefficients in A to 0, and later
				% eliminate the variables). 
				As((num_seg-1)*Navn+nidx, :) = 0; 
			end
			idx = find(sum(A(1:Nvn,Nve+1:Nave),1)==1);	% first segment, no outgoing traffic from dc node
			As(:,Nve+idx) = 0;
			idx = find(sum(A(1:Nvn,Nve+1:Nave),1)==-1); % last segment, no incoming traffic to dc node
			As(:,(num_seg-1)*Nave+Nve+idx) = 0;
			As_flow_temp = block_diag(As, Nf);					% 1
			array.As_flow = As_flow_temp;
			As_proc_temp = block_diag(As_proc_temp, Nf);
			array.As_proc = As_proc_temp;
			array.As_procz = speye(Nd*Nv*Nf);
			
			Ids_temp = spalloc(size(As_flow_temp,1), Nf, 2*Nf);
			nidx_offset = 0;
			for k = 1:Nf
				src = this.FlowTable{k, 'Source'};
				Ids_temp(nidx_offset+src, k) = -1;			%#ok<SPRIX> % 2 only outgoing traffic at source
				nidx_offset = nidx_offset + (num_seg-1)*Navn;
				dest = this.FlowTable{k, 'Target'};
				Ids_temp(nidx_offset+dest, k) = 1;			%#ok<SPRIX> % 3
				nidx_offset = nidx_offset + Navn;
			end
			array.Ids = Ids_temp;

			array.As_load = [repmat(speye(Nve, Nave),1, num_seg*Nf), sparse(Nve, Nd*Nv*Nf);...
				sparse(Nd, Nave*num_seg*Nf), repmat(speye(Nd), 1, Nv*Nf)]; % for load-capacity constraints
			
			b_filter_dc = sparse(false(size(As_proc_temp,1),1));
			nidx_offset = 0;
			for k = 1:Nf
				for j = 1:Nv
					b_filter_dc(nidx_offset+((Nvn+1):Navn)) = true; %#ok<SPRIX>
					nidx_offset = nidx_offset + Navn;
				end
				nidx_offset = nidx_offset + Navn; % last segment: no processing function required
			end
			array.As_proc = array.As_proc(b_filter_dc,:);
			%% TODO
			% VNF type constraint: if a type π of VNF cannot be created on a node
			% i, the corresponding variable z_i(π)=0. Thus, for all flows, the
			% related edge variables x_ji(f,π)=0 and x_ij(f,π+1)=0.
			% Flow-VNF assoication constraint: a flow requires that a type of VNF
			% is deployed on specified nodes. As a result, we eliminate those
			% un-related variables z_i(f,π)=0, also indicating x_ji(f,π)=0 and
			% x_ij(f,π+1)=0. 
			
			%% Compact
			% 1. Reduce the coefficient matrix: limit each flow's distributing
			%		 scope according to its path.
			% 2. Eleminate variables that is predetermined to be 0, which including:
			%		a. incoming flow at source node;
			%		b. outgoing flow at target node;
			%		c. incoming/outgoing flow at the node that cannot create a type
			%			 of VNF.
			b_idx_node_temp = true(Navn*Nvf,1);
			b_idx_edge_temp = true(Nave*Nvf,1);
			if options.bCompact
				nidx_offset = 0;
				eidx_offset = 0;
				for k = 1:Nf
					nidx_off = find(~this.I_flow_node_ex(k,:));
					eidx_off = find(~this.I_flow_edge_ex(k,:));
					for j = 1:num_seg
						b_idx_node_temp(nidx_offset+nidx_off) = false;
						b_idx_edge_temp(eidx_offset+eidx_off) = false;
						nidx_offset = nidx_offset + Navn;
						eidx_offset = eidx_offset + Nave;
					end
				end
			end
			nidx_offset = 0;
			eidx_offset = 0;
			for k = 1:Nf
				src = this.FlowTable{k, 'Source'};
				eidx_off = find(array.As_flow(nidx_offset+src, eidx_offset+(1:Nave))==1); % no incoming flow at source
				b_idx_edge_temp(eidx_offset+eidx_off) = false;
				nidx_offset = nidx_offset + (num_seg-1)*Navn;
				eidx_offset = eidx_offset + (num_seg-1)*Nave;
				dest = this.FlowTable{k, 'Target'};
				eidx_off = find(array.As_flow(nidx_offset+dest, eidx_offset+(1:Nave))==-1); % no outgoing flow at target
				b_idx_edge_temp(eidx_offset+eidx_off) = false;
				nidx_offset = nidx_offset + Navn;
				eidx_offset = eidx_offset + Nave;
			end
			array.As_flow(~b_idx_node_temp, :) = 0;
			array.As_flow(:, ~b_idx_edge_temp) = 0;
			array.As_proc(~b_idx_node_temp(b_filter_dc), :) = 0;
			array.As_proc(:, ~b_idx_edge_temp) = 0;
			array.As_procz(~b_idx_node_temp(b_filter_dc),:) = 0;
			
			problem.Aeq = [array.As_flow, sparse(Navn*Nvf, Nd*Nv*Nf), -array.Ids];
			problem.A = [array.As_proc, -array.As_procz, sparse(Nd*Nv*Nf, Nf)]; % TODO: add processing coefficient "alpha"
			problem.beq = sparse(Navn*Nvf, 1);
			problem.b = sparse(Nd*Nv*Nf, 1);
			indices.I_active_variables = sum(abs([problem.Aeq;problem.A]), 1)~=0;
			problem.numvars = [...
				nnz(indices.I_active_variables(1:this.num_vars_edge));...
				nnz(indices.I_active_variables(this.num_vars_edge+(1:this.num_vars_node)));...
				Nf];
			problem.Aeq = problem.Aeq(:, indices.I_active_variables);
			problem.A = problem.A(:, indices.I_active_variables);
			indices.I_active_rows_eq = sum(abs(problem.Aeq), 2)~=0;
			problem.Aeq = problem.Aeq(indices.I_active_rows_eq,:);
			problem.beq = problem.beq(indices.I_active_rows_eq,:);
			indices.I_active_rows = sum(abs(problem.A), 2)~=0;
			problem.A = problem.A(indices.I_active_rows,:);
			problem.b = problem.b(indices.I_active_rows,:);
			num_varxz = Nave*Nvf+Nd*Nv*Nf;
			array.As_load = array.As_load(:, indices.I_active_variables(1:num_varxz));
			
			indices.I_active_edge_vars = true(this.num_vars_edge,1);
			xoffset = 0;
			for i =1:this.NumberFlows
				for j = 1:num_seg
					indices.I_active_edge_vars(xoffset+((Nve+1):Nave)) = false;
					xoffset = xoffset + Nave;
				end
			end
		end

		%%
		% Override <Slice.priceOptimalFlowRate>.
		% We introduce the quadratic penalty item of the dual-ADMM method.
		% The variables are the flow-edge variables and the flow-node
		% variables.
		%
		% Before solving the problem, call <initializeProblem>.
		function [gamma, fval, output] = priceOptimalFlowRate0(this, x0, lambda, q, options)
			global DEBUG INFO; %#ok<NUSED>
			options = structmerge(options, ...
				getstructfields(this.Parent.options, {'Form', 'OptimizationTool'}, ...
				'default-ignore', {'normal', 'matlab'})); % OptimizationTool = {'matlab', 'cvx'}
			if strcmpi(options.Form, 'compact')
				options.bCompact = true;
			else
				options.bCompact = false;
			end
			Nf = this.NumberFlows;
			Nve = this.NumberVirtualLinks;
			Nd = this.NumberDataCenters;
			
			%%
			if strcmpi(options.OptimizationTool, 'matlab')
				fmincon_opt = optimoptions('fmincon');
				fmincon_opt.Algorithm = 'interior-point';
				%fmincon_opt.Algorithm = 'sqp';
				fmincon_opt.SpecifyObjectiveGradient = true;
				fmincon_opt.Display = 'off';   %'notify-detailed'; %'notify';
				fmincon_opt.OptimalityTolerance = 1e-5;
				fmincon_opt.ConstraintTolerance = 1e-3;
				fmincon_opt.MaxIterations = 100;
				% fmincon_opt.SubproblemAlgorithm = 'cg'; % only take cg steps.
				%fmincon_opt.CheckGradients = true;
				%fmincon_opt.FiniteDifferenceType = 'central';
				if nargin >= 2 && ~isempty(x0)
					var0 = x0;
				elseif isempty(this.x0)
					% set the initial solution.
					%var0 = rand(this.num_flow_vars,1);
					var0 = zeros(this.num_flow_vars,1);
				else
					%var0 = this.x0;
					% var0 = rand(this.num_flow_vars,1);
					var0 = zeros(this.num_flow_vars,1);
				end
				var0 = var0(this.I_active_variables);
				lbs = sparse(length(var0),1);
				if strcmpi(fmincon_opt.Algorithm, 'interior-point')
					fmincon_opt.HessianFcn = ...
						@(x,lbd)FlowEdgeSlice.fcnHessian(x, lbd, this, lambda, q, options);
				end
				warning('off')
				[xs, fval, exitflag, foutput] = ...
					fmincon(@(x)FlowEdgeSlice.fcnProfitPrimal(x, this, lambda, q, options), var0, ...
					this.problem.A, this.problem.b, this.problem.Aeq, this.problem.beq, ...
					lbs, [], [], fmincon_opt);
				warning('on')
				this.interpretExitflag(exitflag, foutput);
			elseif strcmpi(options.OptimizationTool, 'cvx')
				[xs, fval, exitflag, foutput] = cvx_method(lambda, q);
				if exitflag<1
					error('error: CVX failed to solve the problem.');
				end
			else
				error('error: Un-recognized optimization tool.');
			end
			%% CVX
			function [xs, fval, exitflag, output] = cvx_method(lambda, q) 
				w = this.weight;
				p = [this.prices.Link; this.prices.Node];
				plt = options.PenaltyFactor; 
				c = [this.capacities.Link; this.capacities.Node]; 
				Al = this.As_load; 
				cvx_problem = []; cvx_optbnd = []; %#ok<NASGU>
				cvx_optval = []; cvx_status = []; 
				cvx_slvitr = []; cvx_slvtol = []; cvx_cputime = []; %#ok<NASGU>
				if options.bCompact
					c_nx = this.problem.numvars(1);
					c_nz = this.problem.numvars(2);
					c_nr = this.problem.numvars(3);
				else
					c_nx = this.num_vars_edge;
					c_nz = this.num_vars_node;
					c_nr = this.num_vars_rate;
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
			x = zeros(this.num_flow_vars, 1);
			x(this.I_active_variables) = xs;
			nx = this.num_vars_edge;
			nz = this.num_vars_node;
			output.temp_vars.x = x(1:nx);
			output.temp_vars.z = x(nx+(1:nz));
			output.flow_rate = x(nx+nz+(1:Nf));
			output.x0 = x;
			num_varxz = length(xs) - Nf;
			gamma = zeros(this.Parent.NumberDataCenters+this.Parent.NumberLinks,1);
			idx_gamma = [this.VirtualLinks.PhysicalLink; this.Parent.NumberLinks+this.getDCPI];
			gamma(idx_gamma) = ...
				max(0, FlowEdgeSlice.fcnPenalty(xs(1:num_varxz), this, lambda, q, options.PenaltyFactor));
			
% 			output.load.Link = this.As_load(1:Nve,:) * xs(1:num_varxz);
% 			output.load.Node = this.As_load(Nve+(1:Nd),:) * xs(1:num_varxz);
			output.loads = zeros(this.Parent.NumberLinks+this.Parent.NumberDataCenters, 1);
			output.loads(idx_gamma) = [this.As_load(1:Nve,:) * xs(1:num_varxz);
				this.As_load(Nve+(1:Nd),:) * xs(1:num_varxz)];
			output.net_profit = -fval;
			output.iterations = foutput.iterations;
			if strcmpi(options.OptimizationTool, 'cvx')
				output.funcCount = foutput.iterations;
			else
				output.funcCount = foutput.funcCount;
			end
		end

		function [loads, fval, output] = priceOptimalFlowRateDD(this, x0, lambda, options)
			global DEBUG INFO; %#ok<NUSED>
			options = structmerge(options, ...
				getstructfields(this.Parent.options, {'Form', 'OptimizationTool'}, ...
				'default-ignore', {'normal', 'matlab'})); % OptimizationTool = {'matlab', 'cvx'}
			if strcmpi(options.Form, 'compact')
				options.bCompact = true;
			else
				options.bCompact = false;
			end
			Nf = this.NumberFlows;
			Nve = this.NumberVirtualLinks;
			Nd = this.NumberDataCenters;
			
			%%
			if strcmpi(options.OptimizationTool, 'matlab')
				fmincon_opt = optimoptions('fmincon');
				fmincon_opt.Algorithm = 'interior-point';
				%fmincon_opt.Algorithm = 'sqp';
				fmincon_opt.SpecifyObjectiveGradient = true;
				fmincon_opt.Display = 'off';   %'notify-detailed'; %'notify';
				%fmincon_opt.OptimalityTolerance = 1e-5;
				%fmincon_opt.ConstraintTolerance = 1e-4;
				fmincon_opt.MaxIterations = 100;
				% fmincon_opt.SubproblemAlgorithm = 'cg'; % only take cg steps.
				%fmincon_opt.CheckGradients = true;
				%fmincon_opt.FiniteDifferenceType = 'central';
				if nargin >= 2 && ~isempty(x0)
					var0 = x0;
				elseif isempty(this.x0)
					% set the initial solution.
					%var0 = rand(this.num_flow_vars,1);
					var0 = zeros(this.num_flow_vars,1);
				else
					%var0 = this.x0;
					% var0 = rand(this.num_flow_vars,1);
					var0 = zeros(this.num_flow_vars,1);
				end
				var0 = var0(this.I_active_variables);
				lbs = sparse(length(var0),1);
				if strcmpi(fmincon_opt.Algorithm, 'interior-point')
					fmincon_opt.HessianFcn = ...
						@(x,lbd)FlowEdgeSlice.fcnHessianDD(x, lbd, this, lambda, options);
				end
				warning('off')
				[xs, fval, exitflag, foutput] = ...
					fmincon(@(x)FlowEdgeSlice.fcnProfitPrimalDD(x, this, lambda, options), var0, ...
					this.problem.A, this.problem.b, this.problem.Aeq, this.problem.beq, ...
					lbs, [], [], fmincon_opt);
				warning('on')
				this.interpretExitflag(exitflag, foutput);
			else
				error('error: Un-recognized optimization tool.');
			end
		
			%% output solution
			% assert(this.checkFeasible())
			x = zeros(this.num_flow_vars, 1);
			x(this.I_active_variables) = xs;
			nx = this.num_vars_edge;
			nz = this.num_vars_node;
			output.temp_vars.x = x(1:nx);
			output.temp_vars.z = x(nx+(1:nz));
			output.flow_rate = x(nx+nz+(1:Nf));
			output.x0 = x;
			num_varxz = length(xs) - Nf;
			
			loads = zeros(this.Parent.NumberLinks+this.Parent.NumberDataCenters, 1);
			idx = [this.VirtualLinks.PhysicalLink; this.Parent.NumberLinks+this.getDCPI];
			loads(idx) = [this.As_load(1:Nve,:) * xs(1:num_varxz);
				this.As_load(Nve+(1:Nd),:) * xs(1:num_varxz)];
			output.net_profit = -fval;
			output.iterations = foutput.iterations;
			output.funcCount = foutput.funcCount;
		end
		
		function [gamma, fval, x_prox, output] = priceOptimalFlowRatePP(this, x0, lambda, q, xp, options)
			global DEBUG INFO; %#ok<NUSED>
			options = structmerge(options, ...
				getstructfields(this.Parent.options, {'Form', 'OptimizationTool'}, ...
				'default-ignore', {'normal', 'matlab'})); % OptimizationTool = {'matlab', 'cvx'}
			if strcmpi(options.Form, 'compact')
				options.bCompact = true;
			else
				options.bCompact = false;
			end
			if options.PenaltyFactor == 0
				options.PenaltyFactor = this.weight;
			end
			Nf = this.NumberFlows;
			
			%%
			xp = xp(this.I_active_variables);
			if strcmpi(options.OptimizationTool, 'matlab')
				fmincon_opt = optimoptions('fmincon');
				fmincon_opt.Algorithm = 'interior-point';
				%fmincon_opt.Algorithm = 'sqp';
				fmincon_opt.SpecifyObjectiveGradient = true;
				fmincon_opt.Display = 'off';   %'notify-detailed'; %'notify';
				fmincon_opt.OptimalityTolerance = 1e-5;
				fmincon_opt.ConstraintTolerance = 1e-3;
				fmincon_opt.MaxIterations = 100;
				% fmincon_opt.SubproblemAlgorithm = 'cg'; % only take cg steps.
				%fmincon_opt.CheckGradients = true;
				%fmincon_opt.FiniteDifferenceType = 'central';
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
				if strcmpi(fmincon_opt.Algorithm, 'interior-point')
					fmincon_opt.HessianFcn = ...
						@(x,lbd)FlowEdgeSlice.fcnHessianPP(x, lbd, this, lambda, q, options);
				end
				warning('off')
				[xs, fval, exitflag, foutput] = ...
					fmincon(@(x)FlowEdgeSlice.fcnProfitPrimalPP(x, this, lambda, q, xp, options), var0, ...
					this.problem.A, this.problem.b, this.problem.Aeq, this.problem.beq, ...
					lbs, [], [], fmincon_opt);
				warning('on')
				this.interpretExitflag(exitflag, foutput);
			elseif strcmpi(options.OptimizationTool, 'cvx')
				[xs, fval, exitflag, foutput] = cvx_method(lambda, q, xp);
				if exitflag<1
					error('error: CVX failed to solve the problem.');
				end
			else
				error('error: Un-recognized optimization tool.');
			end
			
			%% CVX
			function [xs, fval, exitflag, output] = cvx_method(lambda, q, xp)
				w = this.weight;
				p = [this.prices.Link; this.prices.Node];
				r = options.PenaltyFactor;
				Ns = options.NumberSlices;
				c = [this.capacities.Link; this.capacities.Node];
				Al = this.As_load;
				cvx_problem = []; cvx_optbnd = []; %#ok<NASGU>
				cvx_optval = []; cvx_status = [];
				cvx_slvitr = []; cvx_slvtol = []; cvx_cputime = []; %#ok<NASGU>
				if options.bCompact
					c_nx = this.problem.numvars(1);
					c_nz = this.problem.numvars(2);
					c_nr = this.problem.numvars(3);
				else
					c_nx = this.num_vars_edge;
					c_nz = this.num_vars_node;
					c_nr = this.num_vars_rate;
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
			x = zeros(this.num_flow_vars, 1);
			xs(xs<10^-4*norm(xs)) = 0;
			x(this.I_active_variables) = xs;
			
			nx = this.num_vars_edge;
			nz = this.num_vars_node;
			output.temp_vars.x = x(1:nx);
			output.temp_vars.z = x(nx+(1:nz));
			output.flow_rate = x(nx+nz+(1:Nf));
			output.x0 = x;
			options.TrueValue = true;
			fval = fval + FlowEdgeSlice.fcnProfitPrimalPP(xs, this, lambda, q, xp, options);
			num_varxz = length(xs) - Nf;
			x_prox = x;
			
			gamma = zeros(this.Parent.NumberDataCenters+this.Parent.NumberLinks,1);
			idx = [this.VirtualLinks.PhysicalLink; this.Parent.NumberLinks+this.getDCPI];
			gamma(idx) = ...
				max(0, FlowEdgeSlice.fcnPenalty(xs(1:num_varxz), this, lambda, q, options.PenaltyFactor));
			output.loads = zeros(this.Parent.NumberLinks+this.Parent.NumberDataCenters, 1);
			Nve = this.NumberVirtualLinks;
			Nd = this.NumberDataCenters;
			output.loads(idx) = [this.As_load(1:Nve,:) * xs(1:num_varxz);
				this.As_load(Nve+(1:Nd),:) * xs(1:num_varxz)];
			output.net_profit = -fval;
			output.iterations = foutput.iterations;
			if strcmpi(options.OptimizationTool, 'cvx')
				output.funcCount = foutput.iterations;
			else
				output.funcCount = foutput.funcCount;
			end
		end
		
		function [output, loads, fval] = priceOptimalFlowRate(this, x0, options)
			global DEBUG INFO; %#ok<NUSED>
			options = structmerge(...
				getstructfields(options, 'PricingPolicy', 'default', {'quadratic'}),...
				getstructfields(this.Parent.options, {'Form', 'OptimizationTool'}, ...
				'default-ignore', {'normal', 'matlab'})); % OptimizationTool = {'matlab', 'cvx'}
			if strcmpi(options.Form, 'compact')
				options.bCompact = true;
			else
				options.bCompact = false;
			end
			Nf = this.NumberFlows;
			Nve = this.NumberVirtualLinks;
			Nd = this.NumberDataCenters;
			
			%%
			if strcmpi(options.OptimizationTool, 'matlab')
				fmincon_opt = optimoptions('fmincon');
				fmincon_opt.Algorithm = 'interior-point';
				%fmincon_opt.Algorithm = 'sqp';
				fmincon_opt.SpecifyObjectiveGradient = true;
				fmincon_opt.Display = 'off';   %'notify-detailed'; %'notify';
				%fmincon_opt.OptimalityTolerance = 1e-5;
				fmincon_opt.ConstraintTolerance = 1e-3;
				fmincon_opt.MaxIterations = 100;
				% fmincon_opt.SubproblemAlgorithm = 'cg'; % only take cg steps.
				%fmincon_opt.CheckGradients = true;
				%fmincon_opt.FiniteDifferenceType = 'central';
				if nargin >= 2 && ~isempty(x0)
					var0 = x0;
				elseif isempty(this.x0)
					% set the initial solution.
					%var0 = rand(this.num_flow_vars,1);
					var0 = zeros(this.num_flow_vars,1);
				else
					%var0 = this.x0;
					% var0 = rand(this.num_flow_vars,1);
					var0 = zeros(this.num_flow_vars,1);
				end
				var0 = var0(this.I_active_variables);
				lbs = sparse(length(var0),1);
				if strcmpi(fmincon_opt.Algorithm, 'interior-point')
					fmincon_opt.HessianFcn = ...
						@(x,lbd)FlowEdgeSlice.fcnHessianCC(x, lbd, this);
				end
				if isfield(options, 'CapacityConstrained') && options.CapacityConstrained
					A = [this.problem.A; this.As_load, sparse(size(this.As_load,1),Nf)];
					b = [this.problem.b; this.capacities.Link; this.capacities.Node];
				else
					A = this.problem.A;
					b = this.problem.b;
				end
				warning('off')
				[xs, fval, exitflag, foutput] = ...
					fmincon(@(x)FlowEdgeSlice.fcnProfitPrimalCC(x, this, options), var0, ...
					A, b, this.problem.Aeq, this.problem.beq, ...
					lbs, [], [], fmincon_opt);
				warning('on')
				this.interpretExitflag(exitflag, foutput);
			else
				error('error: Un-recognized optimization tool.');
			end
			
			%% output solution
			% assert(this.checkFeasible())
			x = zeros(this.num_flow_vars, 1);
			xs(xs<10^-4*norm(xs)) = 0;
			x(this.I_active_variables) = xs;
			nx = this.num_vars_edge;
			nz = this.num_vars_node;
			output.temp_vars.x = x(1:nx);
			output.temp_vars.z = x(nx+(1:nz));
			output.flow_rate = x(nx+nz+(1:Nf));
			output.x0 = x;
			num_varxz = length(xs) - Nf;
			
			if nargout >= 2
				loads = zeros(this.hs.Parent.NumberLinks+this.hs.Parent.NumberDataCenters, 1);
				idx = [this.VirtualLinks.PhysicalLink; this.Parent.NumberLinks+this.getDCPI];
				loads(idx) = [this.As_load(1:Nve,:) * xs(1:num_varxz);
					this.As_load(Nve+(1:Nd),:) * xs(1:num_varxz)];
			end
			output.net_profit = -fval;
			output.iterations = foutput.iterations;
			output.funcCount = foutput.funcCount;
		end
	end
	
	methods(Access = protected)
		function cp = copyElement(this)
			cp = copyElement@Slice(this);
			cp.augraph = this.augraph.copy();
		end
		
		function [tf, vars] = postProcessing(this)
			tf = true;
			this.Variables.x = this.temp_vars.x;
			this.Variables.z = this.temp_vars.z;
			if nargout >= 2
				vars = [this.Variables.x; this.Variables.z; this.flow_rate];
			end
		end
		
		%%
		% Override <Slice.getLinkLoad>
		% edge_vars: 
		function ye = getLinkLoad(this, edge_vars, bFinal)
			% retrive the link load of the slice, given the path variables.
			if nargin == 1 || isempty(edge_vars)
				edge_vars = this.Variables.x;
			end
			if nargin <= 2
				bFinal = false;
			end
			Nve = this.NumberVirtualLinks;
			if bFinal
				% length(edge_vars) == nnz(this.I_active_edge_vars)
				full_edge_vars = sparse(this.num_vars_edge, 1);
				b_active_edges = this.I_active_edge_vars(1:this.num_vars_edge);
				full_edge_vars(b_active_edges) = edge_vars;
				edge_vars = full_edge_vars;
			end
			nx = this.problem.numvars(1);
			As = this.As_load(1:Nve, 1:nx);
			if length(edge_vars) == nx
				ye = full(As * edge_vars);
			else
				if length(edge_vars) == nnz(this.I_active_edge_vars)
					full_edge_vars = sparse(this.num_vars_edge, 1);
					full_edge_vars(this.I_active_edge_vars) = edge_vars;
					edge_vars = full_edge_vars;
				end
				ye = full(As * edge_vars(this.I_active_variables(1:this.num_vars_edge)));				
			end
		end

		%%
		% Override <Slice.getNodeLoad>
		function vn = getNodeLoad(this, node_vars)
			if nargin == 1 || isempty(node_vars)
				node_vars = this.Variables.z;
			end
			Nd = this.NumberDataCenters;
			Nve = this.NumberVirtualLinks;
			nx = this.problem.numvars(1);
			nz = this.problem.numvars(2);
			if length(node_vars) > nz
				node_vars = node_vars(this.I_active_variables(this.num_vars_edge+(1:this.num_vars_node)));
			end
			As = this.As_load(Nve+(1:Nd),nx+(1:nz));
			vn = full(As*node_vars);
		end
		
	end
	
	methods (Access = protected, Static)
		%% Objective function
		% See also <>
		% NOTE: 'vars' might be compressed.
		function [profit, grad] = fcnProfitPrimal(vars, slice, lambda, q, options)
			plt = options.PenaltyFactor;
			base_price = [slice.prices.Link; slice.prices.Node];
			
			Nf = slice.NumberFlows;
			Nv = slice.NumberVNFs;
			Nvf = slice.NumberFlowSections;
			
			num_varxz = length(vars) - Nf;
			flow_rate = vars((num_varxz+1):end);
			load = slice.As_load * vars(1:num_varxz);
						
			profit = -slice.weight*sum(fcnUtility(flow_rate)) ...
				+ dot(base_price, load);
			%% TODO
			penalty_item = FlowEdgeSlice.fcnPenalty(vars(1:num_varxz), slice, lambda, q, plt);
			profit = profit	+ plt/2*(norm(max(0,penalty_item)))^2;			
			if isfield(options, 'bFinal') && options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				grad = sparse(slice.num_flow_vars, 1);
				grad(slice.I_active_edge_vars) = repmat(slice.prices.Link,Nvf,1);
				grad(slice.num_vars_edge+(1:slice.num_vars_node)) = repmat(slice.prices.Node,Nf*Nv,1);
				grad = grad(slice.I_active_variables);
				%% update gradient on (x,z)
				idx = penalty_item<=0;
				A = slice.As_load;
				A(idx, :) = 0;
				grad(1:num_varxz) = grad(1:num_varxz) + (A')*max(0, penalty_item);
				%% update gradient on r
				grad((num_varxz+1):end) = -slice.weight./(1+flow_rate);
			end
		end
		
		%%
		% Override <Slice.fcnHessian>.
		function h = fcnHessian(vars, ~, slice, lambda, q, options) 
			num_vars = length(vars);
			num_varxz = num_vars - slice.NumberFlows;
			flow_rate = vars((num_varxz+1):end);
			
			h = sparse(num_vars, num_vars);
			h((num_varxz+1):end, (num_varxz+1):end) = spdiag(slice.weight./(1+flow_rate).^2);
			penalty_item = FlowEdgeSlice.fcnPenalty(vars(1:num_varxz), slice, lambda, q, ...
				options.PenaltyFactor);
			idx = penalty_item<=0;
			A = slice.As_load;
			A(idx, :) = 0;
			h(1:num_varxz, 1:num_varxz) = h(1:num_varxz,1:num_varxz) + (A')*A;
		end
		
		%% Penalty function
		function p = fcnPenalty(vars, slice, lambda, q, factor)
			load = slice.As_load * vars;
			capacity = [slice.capacities.Link; slice.capacities.Node];
			p = lambda+1/factor*(load-capacity+q);
		end
		
		
		%% Functions for Dual Decomposition 
		function [profit, grad] = fcnProfitPrimalDD(vars, slice, lambda, options)
			Nf = slice.NumberFlows;
			Nv = slice.NumberVNFs;
			Nvf = slice.NumberFlowSections;
			nx = slice.problem.numvars(1);
			nz = slice.problem.numvars(2);
			nxz = nx + nz;
			
			flow_rate = vars((nxz+1):end);
			load = slice.As_load * vars(1:nxz);
			
			profit = -slice.weight*sum(fcnUtility(flow_rate)) ...
				+ dot([slice.prices.Link; slice.prices.Node]+lambda, load);
			profit = profit + (0.001/2)*sum(vars(1:nxz).^2);
			if isfield(options, 'bFinal') && options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				nl = length(slice.prices.Link);
				grad = sparse(slice.num_flow_vars, 1);
				grad(slice.I_active_edge_vars) = ...
					repmat(slice.prices.Link+lambda(1:nl),Nvf,1);
				grad(slice.num_vars_edge+(1:slice.num_vars_node)) = ...
					repmat(slice.prices.Node+lambda(nl+1:end),Nf*Nv,1);
				grad = grad(slice.I_active_variables);
				%%
				grad(1:nxz) = grad(1:nxz) + 0.001*vars(1:nxz);
				%% update gradient on r
				grad((nxz+1):end) = -slice.weight./(1+flow_rate);
			end
		end
		
		%%
		function h = fcnHessianDD(vars, ~, slice, ~, options) %#ok<INUSD>
			num_vars = length(vars);
			num_varxz = num_vars - slice.NumberFlows;
			flow_rate = vars((num_varxz+1):end);
			
			h = sparse(num_vars, num_vars);
			h(1:num_varxz, 1:num_varxz) = speye(num_varxz)*0.001;
			h((num_varxz+1):end, (num_varxz+1):end) = spdiag(slice.weight./(1+flow_rate).^2);
		end
		
		function [profit, grad] = fcnProfitPrimalPP(vars, slice, lambda, q, xp, options)
			r = options.PenaltyFactor;
			base_price = [slice.prices.Link; slice.prices.Node];
			
			Nf = slice.NumberFlows;
			Nv = slice.NumberVNFs;
			Nvf = slice.NumberFlowSections;
			Ns = options.NumberSlices;
			
			num_varxz = length(vars) - Nf;
			flow_rate = vars((num_varxz+1):end);
			load = slice.As_load * vars(1:num_varxz);
						
			profit = -slice.weight*sum(fcnUtility(flow_rate)) + dot(base_price, load);
			if isfield(options, 'TrueValue') && options.TrueValue
				return;
			end
			profit = profit + 1/(2*r*Ns^2)*sum((vars-xp).^2);
			%% TODO
			penalty_item = FlowEdgeSlice.fcnPenalty(vars(1:num_varxz), slice, lambda, q, r);
			profit = profit	+ r/2*sum((max(0,penalty_item)).^2);
			if isfield(options, 'bFinal') && options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				grad = sparse(slice.num_flow_vars, 1);
				grad(slice.I_active_edge_vars) = repmat(slice.prices.Link,Nvf,1);
				grad(slice.num_vars_edge+(1:slice.num_vars_node)) = repmat(slice.prices.Node,Nf*Nv,1);
				grad = grad(slice.I_active_variables);
				%% update gradient on r
				grad((num_varxz+1):end) = -slice.weight./(1+flow_rate);
				%% update gradient on (x,z)
				idx = penalty_item<=0;
				A = slice.As_load;
				A(idx, :) = 0;
				grad(1:num_varxz) = grad(1:num_varxz) + (A')*max(0, penalty_item);
				%% gradient on proximal penalty
				grad = grad + 1/(r*Ns^2)*(vars-xp);
			end
		end
		
		function h = fcnHessianPP(vars, ~, slice, lambda, q, options) 
			r = options.PenaltyFactor;
			num_vars = length(vars);
			num_varxz = num_vars - slice.NumberFlows;
			flow_rate = vars((num_varxz+1):end);
			Ns = options.NumberSlices;
			
			h = 1/(r*Ns^2)*speye(num_vars);
			h((num_varxz+1):end, (num_varxz+1):end) = ...
				h((num_varxz+1):end, (num_varxz+1):end) + spdiag(slice.weight./(1+flow_rate).^2);
			penalty_item = FlowEdgeSlice.fcnPenalty(vars(1:num_varxz), slice, lambda, q, ...
				options.PenaltyFactor);
			idx = penalty_item<=0;
			A = slice.As_load;
			A(idx, :) = 0;
			h(1:num_varxz, 1:num_varxz) = h(1:num_varxz,1:num_varxz) + (A')*A;
		end
		
		function [profit, grad] = fcnProfitPrimalCC(vars, slice, options)
			base_price = [slice.prices.Link; slice.prices.Node];
			
			Nf = slice.NumberFlows;
			Nv = slice.NumberVNFs;
			Nvf = slice.NumberFlowSections;
			
			num_varxz = length(vars) - Nf;
			flow_rate = vars((num_varxz+1):end);
			load = slice.As_load * vars(1:num_varxz);
						
			profit = -slice.weight*sum(fcnUtility(flow_rate));
			switch options.PricingPolicy
				case {'quadratic-price', 'quadratic'}
					link_load = load(1:slice.NumberVirtualLinks);
					node_load= load((slice.NumberVirtualLinks+1):end);
					[link_payment,link_price_grad] = slice.fcnLinkPricing(slice.prices.Link, link_load);
					[node_payment,node_price_grad] = slice.fcnNodePricing(slice.prices.Node, node_load);
					profit = profit + link_payment + node_payment;
				case 'linear'
					profit = profit + dot(base_price, load);
				otherwise
					error('%s: invalid pricing policy', calledby);
			end
			if isfield(options, 'bFinal') && options.bFinal
				profit = -profit;
			end
			
			%% Gradient
			if nargout == 2
				grad = sparse(slice.num_flow_vars, 1);
				switch options.PricingPolicy
					case {'quadratic-price', 'quadratic'}
						grad(slice.I_active_edge_vars) = repmat(link_price_grad,Nvf,1);
						grad(slice.num_vars_edge+(1:slice.num_vars_node)) = repmat(node_price_grad,Nf*Nv,1);						
					case 'linear'
						grad(slice.I_active_edge_vars) = repmat(slice.prices.Link,Nvf,1);
						grad(slice.num_vars_edge+(1:slice.num_vars_node)) = repmat(slice.prices.Node,Nf*Nv,1);
				end
				grad = grad(slice.I_active_variables);
				%% update gradient on r
				grad((num_varxz+1):end) = -slice.weight./(1+flow_rate);
			end
		end
		
		function h = fcnHessianCC(vars, ~, slice, options)
			num_vars = length(vars);
			num_varxz = num_vars - slice.NumberFlows;
			flow_rate = vars((num_varxz+1):end);
			h = sparse(num_vars, num_vars);
			h((num_varxz+1):end, (num_varxz+1):end) = spdiag(slice.weight./(1+flow_rate).^2);
			if nargin >= 4 && isfield(options, 'PricingPolicy')
					switch options.PricingPolicy
						case {'quadratic-price', 'quadratic'}
							load = slice.As_load * vars(1:num_varxz);
							link_load = load(1:slice.NumberVirtualLinks);
							node_load= load((slice.NumberVirtualLinks+1):end);
							[~,~,de] = this.fcnLinkPricing(this.prices.Link, link_load); % equal to <getLinkCapacity>
							[~,~,dn] = this.fcnNodePricing(this.prices.Node, node_load);
							blocks = blkdiag(block_diag(de, slice.NumberFlows*(slice.NumberVNFs+1)),...
								block_diag(dn,slice.NumberFlows*slice.NumberVNFs));
							blocks(slice.I_active_variables,:) = [];
							blocks(:,slice.I_active_variables) = [];
							h(1:num_varxz,1:num_varxz) = blocks;
					end
			end
		end
		
	end
	
end

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
				this.As_load = sparse(size(this.As_load,1), this.num_vars_edge+this.num_vars_node);
				this.As_load(:,this.I_active_variables(1:this.num_vars_edge+this.num_vars_node)) = As_load_temp;
				[xs_, fval_, exitflag_, foutput_] = ...
					fmincon(@(x)FlowEdgeSlice.fcnProfitPrimal(x, this, lambda, q, options_), var0_, ...
					A, b, Aeq, beq, lbs_, [], [], fmincon_opt);
				this.As_load = As_load_temp;
%}