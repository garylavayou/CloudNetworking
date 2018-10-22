classdef SliceOptimizer
	properties
    Variables;  % Variables after post processing.
		prices;
  end
  
  properties (Access = protected)
    hs;         % reference to the host slice.
		x0;         % start point
		temp_vars = struct; % Temporary variables that directly resulted from optimization.
		flow_rate;          % Temporary results of flow rate.
		I_active_variables logical;  % indicator of active variables
  end
  
  properties(Abstract)
    % Final results of VNF instance capacity on each node, this is configured by
		% inter/intra-slicing, the values are stored in |Variables.v|; This property is a
		% wrapper.
		VNFCapacity;
  end
  
  properties (Dependent)
    NumberVariables int32;	% Number of optimization variables.
    NumberLinearConstraints int32; % Number of linear constraints
  end
  
  methods 
    function this = SliceOptimizer(slice)
      this.hs = slice;
    end
    
		function n = get.NumberVariables(this)
			n = this.get_number_variables();
    end
    
    function n = get.NumberLinearConstraints(this)
			n = this.get_number_linear_constraints();
		end    
  end
  
  methods (Access = protected, Abstract)
    %% Override the method to provide property overriding
    n = get_number_variables(this);
    n = get_number_linear_constraints(this);
  end
  
  methods
    Clear();      % clear variables.
  end
end