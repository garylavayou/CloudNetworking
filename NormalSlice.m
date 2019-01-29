classdef NormalSlice < Slice
	
	%% Constructor
	methods
		function this = NormalSlice(slice_data)
			if nargin == 0
				args = {};
			else
				args = {slice_data};
			end
			this@Slice(args{:});
			if nargin == 0
				return;
			end

		end
		
	end
	
	%% Public Methods
	methods
		function op = getOptimizer(this, options)
			if nargin >= 2 && isfield(options, 'Optimizer')
				Optimizer = options.Optimizer;
			else
				Optimizer = 'NormalSliceOptimizer';
			end
			if nargin == 1
				this.op = instantiateclass(Optimizer, this);
			else
				this.op = instantiateclass(Optimizer, this, options);
			end
			op = this.op;
			op.initializeState();
		end
		
	end
	
end