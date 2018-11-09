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
    profit = getProfit(slice, options);
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
		end
    
    function vc = getVNFCapacity(this, z)
      %       znpf = reshape(full(this.Variables.z), this.NumberServiceNodes, ...
      %       this.NumberPaths, this.NumberVNFs);
      %       znpf = znpf.* full(this.I_dc_path);  % compatible arithmetic operation
      %       this.VNFCapacity = reshape(sum(znpf,2), this.NumberServiceNodes*this.NumberVNFs,1);
      if nargin <= 1
        z = this.op.Variables.z;
      end
      vc = full(this.op.Hdiag * z);
    end
    
    function sc = getCost(this, load, model)
      sc = this.getResourceCost(this, load, model);
    end
    
    function r = getRevenue(this)
      if isempty(this.op.flow_rate)
        r = this.weight*sum(fcnUtility(this.FlowTable.Rate));
      else
        r = this.weight*sum(fcnUtility(this.op.flow_rate));
      end
    end
    
    function tf = isDynamicFlow(~)
      tf = false;
    end
    
    function [omega, sigma, alpha] = utilizationRatio(this)
      n_idx = this.ServiceNodes.Capacity>eps;
      e_idx = this.Links.Capacity>eps;
      c_node = sum(this.ServiceNodes.Capacity(n_idx));
      c_link = sum(this.Links.Capacity(e_idx));
      alpha = [c_node c_link]./(c_node+c_link);
      theta_v = sum(this.ServiceNodes.Load(n_idx))/c_node;
      theta_l = sum(this.Links.Load(e_idx))/c_link;
      omega = dot(alpha, [theta_v, theta_l]);
      
      if nargout == 2
        sigma = std([this.Links.Load(e_idx)./this.Links.Capacity(e_idx);...
          this.ServiceNodes.Load(n_idx)./this.ServiceNodes.Capacity(n_idx)]);
      end
    end
  end
  
  %% Protected Methods
  methods (Access=protected)
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
