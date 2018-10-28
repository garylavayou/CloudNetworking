classdef NormalSliceOptimizer < SliceOptimizer

  properties (Access={?NormalSlice})
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
		% mask of the original flow-edge variables;
		I_active_edge_vars logical;
		% 
		I_active_rows logical;
		I_active_rows_eq logical;
		%% Problem coefficients
		As_flow double;		% coefficient for flow reservation constraints
		Ids double;
		As_proc double;			% processing resource requirements
		As_procz double;
		As_load double;
		problem = struct('Aeq', [], 'A', [], 'beq', [], 'b', []);
  end

  methods
    initializeProblem(this);
    
    function setProblem(this, capacities, array, problem, indices)
      if nargin >= 2 && ~isempty(capacities)
        setProblem@SliceOptimizer(capacities);
      end
			if nargin >= 3 && ~isempty(array)
				this.As_flow = array.As_flow;
				this.As_proc = array.As_proc;
				this.As_procz = array.As_procz;
				this.As_load = array.As_load;
				this.Ids = array.Ids;
			end
			if nargin >= 4 && ~isempty(problem)
				this.problem = problem;
			end
			if nargin >= 5 && ~isempty(indices)
				this.I_active_variables = indices.I_active_variables;
				this.I_active_edge_vars = indices.I_active_edge_vars;
				this.I_active_rows = indices.I_active_rows;
				this.I_active_rows_eq = indices.I_active_rows_eq;
			end
		end
  end
  

  methods

  end
end