classdef SimpleSliceOptimizer < SliceOptimizer
  properties
		As_res;   % Coefficient matrix of the processing constraints.
		Hrep;     % Coefficient used to compute node resource consumption.
		Hdiag;		% Coefficient used to compute VNF intance capacity.
    num_varz; % number of variables in vector |z_npf|.
    
  end

  properties (Dependent)
		VNFCapacity;
    NumberVariables int32;	% Number of optimization variables.
    NumberLinearConstraints int32; % Number of linear constraints
  end
  
  methods
    function this = SimpleSliceOptimizer(slice)
      this@SliceOptimizer(slice);

      this.getAs_res;
      this.getHrep;
      this.getHdiag;
    end
  end
  
  methods
    function c = get.VNFCapacity(this)
      if ~isfield(this.Variables, 'v') || isempty(this.Variables.v)
        warning('VNF capacity not set, set to VNF load.');
        c = this.getVNFCapacity;
        this.Variables.v = c;
      else
        c = this.Variables.v;
      end
    end
    
		function n = get.num_varz(this)
			n = this.NumberVNFs*this.NumberServiceNodes*this.NumberPaths;
		end
    
  end
  
  methods (Access = protected)
		function n = get_number_variables(this)
			n = (this.hs.NumberVNFs*this.hs.NumberServiceNodes+1)*this.hs.NumberPaths;
    end
    
    function n = get_number_linear_constraints(this)
			n = size(this.As_res,1);
    end

    
  end
end