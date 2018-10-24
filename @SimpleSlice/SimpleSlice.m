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
        return;
      end
      this@Slice(slice_data);
      this.op = SimpleSliceOptimizer(this);
    end
    
  end
  
  methods
    %% Capacity and Load
    % Public interface for network to inquire the resource occupation of the slices.
    % In class <Slice>, <getLinkCapacity> and <getNodeCapacity> are equal to the
    % protected methods <getLinkLoad> and <getNodeLoad> respectively. But in
    % subclasses of <Slice>, the slice load might be less than its capacity, so that
    % the two group of methods return different results.
    function c = getLinkCapacity(this, isfinal)
      if nargin == 1 || isfinal
        c = this.getLinkLoad;
      else
        c = this.getLinkLoad(this.op.temp_vars.x);
      end
    end
    
    function c = getNodeCapacity(this, isfinal)
      if nargin == 1 || isfinal
        c = this.getNodeLoad;
      else
        c = this.getNodeLoad(this.op.temp_vars.z);
      end
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
    function setPathBandwidth(this, x)
      if nargin == 1
        x = this.Variables.x;
      end
      p = 1;
      for i = 1:this.NumberFlows
        pathlist = this.FlowTable{i,'Paths'}.paths;
        for l = 1:length(pathlist)
          pathlist{l}.bandwidth = x(p);
          p = p + 1;
        end
      end
    end
    
    %% Get Load
    % <getLinkLoad> and <getNodeLoad> are only used inside the <Slice> class. The
    % substrate network cares how much resources are ocuppied by the slice, that is
    % the capacity of resources. See also <getLinkCapacity> and <getNodeCpacity>.
    function ye = getLinkLoad(this, path_vars)
      % retrive the link load of the slice, given the path variables.
      if nargin == 1
        ye = this.op.I_edge_path * this.op.Variables.x;
      else
        ye = this.op.I_edge_path * path_vars;
      end
      ye = full(ye);
    end
    
    %%%
    % retrive the node load of the slice, given the node variables. The node variables
    % represent the resource allocation of data center nodes.
    %
    % |v_n|: data center's resource consumption.
    function v_n = getNodeLoad(this, node_vars)
      if nargin == 1
        node_vars = this.op.Variables.z;
      end
      %%
      % |node_vars| is index by |(node,path,function)|.
      v_n = full(this.op.Hrep*node_vars);
      % assert(isempty(find(abs(this.Hrep*node_vars-v_n)>10^-6,1)), 'error: unequal node load');
      
      %% Alternative way to compute the node load.
      % _Hrep_ replace the following procedure.
      %             node_load = sum(f, node_vars(:,:,f).*I_dc_path).
      %             Nsn = this.op.NumberServiceNodes;
      %             Np = this.op.NumberPaths;
      %             v_n = zeros(Nsn,1);
      %             np = Nsn * Np;
      %             z_index = 1:np;
      %             for i = 1:this.NumberVNFs
      %                 node_vars_fi = reshape(node_vars(z_index), Nsn, NP);
      %                 v_n = v_n + sum(this.op.I_dc_path.*node_vars_fi,2);
      %                 z_index = z_index + np;
      %             end
    end
  end
  
  %% Static Methods
  methods (Static)
    %% fcnProfit
    % Evalute the objective function and gradient.
    %
    %      [profit, grad] = fcnProfit(vars, slice, options)
    %
    % |grad|: the gradient value of the objective function.
    % The upper bound number of non-zero elements in the gradient vector: the gradient on path
    % variable is nonzeros, so there is |P| components; whether the gradient on node variable
    % is zeros is decided by the node-path incidence matrix, i.e. |nnz(I_dc_path)*F|.
    [profit, grad] = fcnProfit(vars, slice, options);
    %% fcnProfitCompact
    % Evaluate the objective function and gradient
    %
    %     [profit, grad] = fcnProfitCompact(act_vars, slice, options)
    %
    % only active independent variables are passed into the objective function.
    % Considering the constraint's coefficient matrix, if the corresponding column of the
    % coefficient matrix for a variable is all zero, then this variable is inactive and can be
    % directly set as 0. So we can remove it from the optimization problem.
    %
    % NOTE: we can also remove the all-zero rows of the coefficient matrix, which do not
    % influence the number of variables. See also <optimalFlowRate>.
    function [profit, grad] = fcnProfitCompact(act_vars, slice, options)
      vars = zeros(options.num_orig_vars,1);
      vars(slice.I_active_variables) = act_vars;
      
      % we extend the active variables by adding zeros to the inactive ones.
      if nargout <= 1
        profit = SimpleSlice.fcnProfit(vars, slice, options);
      else
        [profit, grad] = SimpleSlice.fcnProfit(vars, slice, options);
        % eliminate the inactive variable's derivatives.
        grad = grad(slice.I_active_variables);
      end
    end
    
    [profit, grad] = fcnSocialWelfare(x_vars, S, options);
    function [profit, grad] = fcnSocialWelfareCompact(act_vars, slice)
      vars = zeros(options.num_orig_vars,1);
      vars(slice.I_active_variables) = act_vars;
      
      if nargout <= 1
        profit = SimpleSlice.fcnSocialWelfare(vars, slice);
      else
        [profit, grad] = SimpleSlice.fcnSocialWelfare(vars, slice);
        grad = grad(slice.I_active_variables);
      end
    end
    
    %% fcnHessian
    % Hessian matrix of the Largrangian.
    %
    %      hs = fcnHessian(var_x, ~, slice, options)
    %
    % Since the problem only contains linear constraint, the hessian matrix of the
    % Largrangian is equal to the second derivatives of the objective function, and the
    % Largrangian multipliers $\lambda$ takes no effect.
    % The Hessian matrix contains only $P^2$ nonzeros elements on the diagonal,
    % which is the second derviatives on path variables.
    hs = fcnHessian(var_x, ~, slice, options);
    
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
    % this part should not be ignored. See also <getNetworkCost> and <getStaticCost>.
    % The calculation is not absolutely precise, since it cannot be decide that the
    % static cost should be arributed to which slices.
  end
  
end
