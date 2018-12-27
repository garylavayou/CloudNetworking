%%  Network Slice
% Support network resource allocation scheme.
%%
classdef SimpleSlice < Slice
  properties (Access = protected)
    % 		link_load;          % temporary results of link load.
    % 		node_load;          % temporary results of node load.
  end
  properties(Dependent, GetAccess={?CloudNetwork})
    % local_path_id;
	end
  
  methods
    function this = SimpleSlice(slice_data)
      if nargin == 0
        args = cell(0);
			else
				args = {slice_data};
      end
      this@Slice(args{:});
    end
    
	end
  
	%% Public Methods
  methods
		function op = getOptimizer(this, options)
			if nargin == 1
				this.op = SimpleSliceOptimizer(this);
			else
				this.op = SimpleSliceOptimizer(this, options);
			end
			op = this.op;
			op.initializeState();
		end
    
	end
  
  methods (Access = private)
    %%
    %         function sc = getResourceCost(this, load, model)
    %             if nargin <= 1 || isempty(load)
    %                 load.Node = this.Nodes.Load;
    %                 load.Link = this.Links.Load;
    %             end
    %             if nargin <= 2
    %                 warning('model is set as Approximate.');
    %                 model = 'Approximate';
    %             end
    %
    %             pn = this.Parent;
    %             link_uc = pn.getLinkCost(this.Links.PhysicalLink); % the virtual links's unit cost
    %             node_uc = pn.getNodeCost(this.getDCPI); % the virtual nodes's unit cost
    %             epsilon = pn.unitStaticNodeCost;
    %
    %             if strcmp(model, 'Approximate')
    %                 sc = dot(link_uc, load.Link) + dot(node_uc, load.Node) ...
    %                     + pn.phis_n*sum(load.Node)+pn.phis_l*sum(load.Link)+...
    %                     pn.static_factor*epsilon/pn.NumberSlices;
    %             elseif strcmp(model, 'Accurate')
    %                 sc = dot(link_uc, load.Link) + dot(node_uc, node_load);
    %             else
    %                 error('error: invalid model %s', model);
    %             end
    %         end
    %%%
    % When compute the static cost, the Capacity of all physical nodes and links is
    % included.
    % |epsilon/pn.NumberSlices| is a constant. To keep consistence with other methods,
    % this part should not be ignored. See also <totalCost> and <getStaticCost>.
    % The calculation is not absolutely precise, since it cannot be decide that the
    % static cost should be arributed to which slices.
  end
  
end
