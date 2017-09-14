classdef SliceEx < Slice
    properties
        NumberDataCenters;
        num_vars;
    end
    
    methods
        function n = get.NumberDataCenters(this)
            n = nnz(this.VirtualDataCenters);
        end
        function n = get.num_vars(this)
            n = (this.NumberVNFs*this.NumberDataCenters+1)*this.NumberPaths;
        end
        %%%
        % retrive the node load of the slice, given the node variables. The node variables
        % represent the resource allocation of data center nodes.
        %
        % |v_n|: data center's resource consumption.
        function v_n = getNodeLoad(this, node_vars)
            if nargin == 1
                node_vars = this.Variables.z;
            end
            %%
            % |node_vars| is index by |(node,path,function)|.
            % node_load = sum(f, node_vars(:,:,f).*I_node_path).
            NC = this.NumberDataCenters;
            NP = this.NumberPaths;
            v_n = zeros(NC,1);
            np = NC * NP;
            z_index = 1:np;
            for i = 1:this.NumberVNFs
                node_vars_fi = reshape(node_vars(z_index), NC, NP);
                v_n = v_n + sum(this.I_node_path.*node_vars_fi,2);
                z_index = z_index + np;
            end
            
            %% Alternative way to compute the node load.
            %             col_index = (1:NC:((NP-1)*NC+1))';
            %             col_index = repmat(col_index, 1, NV);
            %             for c = 2:NV
            %                 col_index(:,c) = col_index(:,c-1) + NC*NP;
            %             end
            %             col_index = col_index(:);
            %             As = zeros(NC, NP*NV);
            %             for row_index = 1:NC
            %                 As(row_index, col_index) = repmat(this.I_node_path(row_index,:),1, NV);
            %                 col_index = col_index + 1;
            %             end
            %             vn = As*node_vars;
        end
        
        function pc = getPathCost(this, lambda_e, lambda_n)
            link_uc = this.link_unit_cost;
            node_uc = this.node_unit_cost;
            pc = cell(this.NumberFlows, 1);
            p = 1;
            for i = 1:this.NumberFlows
                num_paths = this.FlowTable{i,'Paths'}.Width;
                pc{i} = zeros(num_paths, 1);
                for j = 1:num_paths
                    eid = this.I_edge_path(:,p)~=0;
                    pc{i}(j) = pc{i}(j) + sum(link_uc(eid)+lambda_e(eid)) + this.Parent.phis_l*nnz(eid);
                    nid = this.I_node_path(:,p)~=0;
                    pc{i}(j) = pc{i}(j) + min(node_uc(nid)+lambda_n(nid)) + this.Parent.phis_n;
                    p = p + 1;
                end
            end
        end
        
        %%%
        % getSliceCost  When compute the static cost, the Capacity of all physical nodes
        % and links is included.
        % |epsilon/pn.NumberSlices| is a constant. To keep consistence with other methods,
        % this part should not be ignored. See also getNetworkCost and getStaticCost.
        % The calculation is not absolutely precise, since it cannot be decide that the
        % static cost should be arributed to which slices.
        %
        % When calculate network cost as a single slice, this method equals to
        % _getNetworkCost_ .
        function sc = getSliceCost(this, node_load, link_load, model)
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
        %% Get the node cost on a path.
        function [nc, nid] = getPathNodeCost(this, pid, lambda_n)
            node_uc = this.node_unit_cost;
            nid = find(this.I_node_path(:,pid)~=0);
            nc = node_uc(nid) + lambda_n(nid) + this.Parent.phis_n;
        end


    end
end