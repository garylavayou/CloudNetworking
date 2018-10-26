classdef SliceEx < SimpleSlice
    properties
        constant_profit = 0;	% consant profit added to avoid slices with negative profit.
    end
    
    methods
        function this = SliceEx(slice_data)
            if nargin == 0
                slice_data = [];
            end
            this@SimpleSlice(slice_data);
            if nargin == 0
                return;
            end
            if isfield(slice_data, 'ConstantProfit')
                this.constant_profit = slice_data.ConstantProfit;
            end
        end
        
        %         function pc = getPathCost(this, lambda_e, lambda_n)
        %             link_uc = this.Parent.readLink('UnitCost', this.VirtualLinks.PhysicalLink); % the virtual links's unit cost
        %             node_uc = this.Parent.readNode('UnitCost', this.VirtualNodes.PhysicalNode); % the virtual nodes's unit cost
        %             pc = cell(this.NumberFlows, 1);
        %             p = 1;
        %             for i = 1:this.NumberFlows
        %                 num_paths = this.FlowTable{i,'Paths'}.Width;
        %                 pc{i} = zeros(num_paths, 1);
        %                 for j = 1:num_paths
        %                     eid = this.I_edge_path(:,p)~=0;
        %                     pc{i}(j) = pc{i}(j) + sum(link_uc(eid)+lambda_e(eid)) + this.Parent.phis_l*nnz(eid);
        %                     nid = this.I_dc_path(:,p)~=0;
        %                     pc{i}(j) = pc{i}(j) + min(node_uc(nid)+lambda_n(nid)) + this.Parent.phis_n;
        %                     p = p + 1;
        %                 end
        %             end
        %         end
        function pc = getPathCost(this, lambda_e, lambda_n) 
            link_uc = this.LinkCost;
            node_uc = this.NodeCost;
            pc = cell(this.NumberFlows, 1);
            p = 1;
            for i = 1:this.NumberFlows
                num_paths = this.FlowTable{i,'Paths'}.Width;
                pc{i} = zeros(num_paths, 1);
                for j = 1:num_paths
                    eid = this.I_edge_path(:,p)~=0;
                    pc{i}(j) = pc{i}(j) + sum(link_uc(eid)+lambda_e(eid)) + this.Parent.phis_l*nnz(eid);
                    nid = this.I_dc_path(:,p)~=0;
                    pc{i}(j) = pc{i}(j) + min(node_uc(nid)+lambda_n(nid)) + this.Parent.phis_n;
                    p = p + 1;
                end
            end
        end
        
        %%%
        % When compute the static cost, the Capacity of all physical nodes and links is included.
        % |epsilon/pn.NumberSlices| is a constant. To keep consistence with other methods,
        % this part should not be ignored. See also <totalCost> and <getStaticCost>.
        % The calculation is not absolutely precise, since it cannot be decide that the
        % static cost should be arributed to which slices.
        %
        % When calculate network cost as a single slice, this method equals to
        % _getNetworkCost_ .
        function sc = getCost(this, node_load, link_load, model)
            if nargin <= 1 || isempty(node_load)
                node_load = this.VirtualDataCenters.Load;
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.VirtualLinks.Load;
            end
            if nargin <= 3
                warning('model is set as Approximate.');
                model = 'Approximate';
            end
            
            pn = this.Parent;
            link_uc = this.link_unit_cost;
            node_uc = this.node_unit_cost;
            epsilon = pn.unitStaticNodeCost;
            
            if strcmp(model, 'Approximate')
                %% A temporary slice should be assigned the parent network
                % so that |this.Parent| is valid.
                %
                % $$static\_cost =
                %   \frac{\epsilon\delta(N-1)\sum_{n\in\mathcal{N}^{(s)}}{v_n}}{\sum_{n\in\mathcal{N}}{V_n}}
                % + \frac{\epsilon(1-\delta)(N-1)\sum_{e\in\mathcal{L}^{(s)}}{y_e}}{\sum_{e\in\mathcal{L}}{C_e}}
                % + \frac{\epsilon}{|\mathcal{S}|}$$
                sc = dot(link_uc, link_load) + dot(node_uc, node_load) ...
                    + pn.phis_n*sum(node_load)+pn.phis_l*sum(link_load)+...
                    pn.static_factor*epsilon/pn.NumberSlices;
                %%%
                % equal to:
                %
                %   sc = dot(link_uc+pn_phis_l, link_load) + dot(node_uc+pn.phis_n, node_load) + ...
                %       pn.static_factor*epsilon/pn.NumberSlices;
            elseif strcmp(model, 'Accurate')
                %% Accurate Model
                % A slice cannot decide how to devide the static cost between slices.
                % one method is devision by usage, using the physic network's load data
                % of each slice. So this function only calculate the dynamic part of the
                % cost, and the slices should further calculate the static cost outside
                % this method.
                sc = dot(link_uc, link_load) + dot(node_uc, node_load);
            else
                error('error: invalid model %s', model);
            end
        end
        
        %         function [nc, nid] = getPathNodeCost(this, pid, lambda_n)
        %            node_uc = this.Parent.readNode('UnitCost', this.VirtualNodes.PhysicalNode);
        %            nid = find(this.I_dc_path(:,pid)~=0);
        %            nc = node_uc(nid) + lambda_n(nid) + this.Parent.phis_n;
        %         end
        %% Get the node cost on a path.
        function [nc, nid] = getPathNodeCost(this, pid, lambda_n)
            node_uc = this.node_unit_cost;
            nid = find(this.I_dc_path(:,pid)~=0);
            nc = node_uc(nid) + lambda_n(nid) + this.Parent.phis_n;
        end

        function r = getRevenue(this)
            r = getRevenue@SimpleSlice(this) + this.constant_profit;
        end
        
        [fval, node_load, link_load] = subproblemNetSocialWelfare( this, lambda );
        [fval, node_load, link_load] = subproblemNetSocialWelfare2( this, lambda );

    end
    
    methods(Static)
        % * *fcnProfitCompact* |static| : Evalute the objective function and gradient.
        %
        %      [profit, grad]= fcnProfitCompact(act_var_x, S)
        %
        % Only active independent variables are passed into the objective function.
        %
        [profit, grad]= fcnProfitCompact(act_vars, slice);
        [fval,  grad] = subproblemObjective(vars, lambda, slice);
        %%%
        % * *fcnHessianCompact* |static| : Compact form of Hessian matrix of the Largrangian.
        %
        %      hess = fcnHessianCompact(act_var_x, ~, S)
        hess = fcnHessianCompact(act_var_x, ~, S);
    end
    
    methods(Access=private)
        function rc = getResourceCost(this, node_load, link_load, model)
            if nargin <= 1 || isempty(node_load)
                node_load = this.VirtualDataCenters.Load;
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.VirtualLinks.Load;
            end
            if nargin <= 3
                warning('model is set as Approximate.');
                model = 'Approximate';
            end
            rc = getResourceCost@SimpleSlice(this, node_load, link_load, model);
            epsilon = pn.unitStaticNodeCost;
            switch model
                case 'Approximate'
                    %%
                    % $$static\_cost =
                    %   \frac{\epsilon\delta(N-1)\sum_{n\in\mathcal{N}^{(s)}}{v_n}}{\sum_{n\in\mathcal{N}}{V_n}}
                    % + \frac{\epsilon(1-\delta)(N-1)\sum_{e\in\mathcal{L}^{(s)}}{y_e}}{\sum_{e\in\mathcal{L}}{C_e}}
                    % + \frac{\epsilon}{|\mathcal{S}|}$$
                    rc = rc + pn.static_factor*epsilon/pn.NumberSlices;
                    %%%
                    % equal to:
                    %
                    %   sc = dot(link_uc+pn_phis_l, link_load) + dot(node_uc+pn.phis_n,
                    %   node_load) + ... 
                    %       pn.static_factor*epsilon/pn.NumberSlices; 
            end
        end
    end
    
    methods(Access=protected)
        function [x, fval, exitflag] = optimize(this, params, options)
            if options.SlicingMethod.IsPricing
							[x, fval, exitflag] = optimize@SimpleSlice(this, params, options);
						else
							[x, fval, exitflag] = ...
								fmincon(@(x)SimpleSlice.fcnSocialWelfare(x,this,'Approximate'), ...
								params.x0, params.As, params.bs, params.Aeq, params.beq, ...
								params.lb, params.ub, [], fmincon_opt);
						end
        end
    end
end
