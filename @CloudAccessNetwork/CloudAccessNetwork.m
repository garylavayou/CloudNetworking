%% CloudAccessNetwork 
% Incorporate access network into cloud network.
% NOTE: CloudAccessNetwork can be derived from <CloudNetwork> and <AccessNetwork> (TODO),
% which are both derived from <PhysicalNetwork>.
classdef CloudAccessNetwork < CloudNetwork
    properties (GetAccess = public, SetAccess = protected)
        base_stations;      % node index of base stations
    end 
    properties(Dependent)
        NumberBaseStations;
    end
    
    methods
        function this = CloudAccessNetwork(node_opt, link_opt, VNF_opt, options)
            %% Initialize as SD-RAN
            % node is indexed by the order: router, GBS NBS.
            net = TwoLayerNetwork('SD-RAN');
            dc.core = [2; 3; 6; 7; 9; 12]; 
            dc.gateway = net.index.GBS;
            %% Set link capacity
            if isfield(link_opt, 'CapacityFactor')
                warning('CapacityFactor (%.2f) changes the link capacity configuration.',...
                    link_opt.CapacityFactor);
            end
            link_opt.Capacity = zeros(net.Size);
            link_opt.link_cost = zeros(net.Size);
            link_capacity_map = [1000, 1000, 300, 100];
            link_cost_map = [1, 2, 4, 6];
            link_opt.cost = LinkCostOption.NetworkSpecified;
            for eid = 1:length(net.link)
                link_type = net.link(eid, 3);
                if link_type <= 3
                    link_capacity = link_capacity_map(1+link_type);
                    link_cost = link_cost_map(1+link_type);
                else
                    link_capacity = 0;
                    link_cost = inf;
                end
                link_opt.Capacity(net.link(eid,1),net.link(eid,2)) = link_capacity;
                link_opt.Capacity(net.link(eid,2),net.link(eid,1)) = link_capacity;
                link_opt.link_cost(net.link(eid,1),net.link(eid,2)) = link_cost;
                link_opt.link_cost(net.link(eid,2),net.link(eid,1)) = link_cost;
            end
            [~, ~, link_opt.link_cost] = find(link_opt.link_cost);
%             if isfield(link_opt, 'CostUnit')
%                 link_opt.link_cost = link_opt.link_cost * link_opt.CostUnit;
%             end
            node_opt.model = NetworkModel.SD_RAN;
            node_opt.location = net.node(:,1:2);
            node_opt.node_capacity = zeros(net.Size,1);
            node_opt.node_capacity(dc.core) = 5000;
            node_opt.node_capacity(dc.gateway) = 500;
            if isfield(node_opt, 'capacity_factor')
                warning('CapacityFactor (%.2f) changes the node capacity configuration.',...
                    node_opt.capacity_factor);
            end
            node_opt.capacity = NodeCapacityOption.NetworkSpecified;
            node_opt.cost = NodeCostOption.NetworkSpecified;
            node_opt.node_cost = inf * ones(net.Size,1);
            node_opt.node_cost(dc.core) = 1;
            node_opt.node_cost(dc.gateway) = 3;
            
%             if isfield(node_opt, 'CostUnit')
%                 node_opt.node_cost = node_opt.node_cost * node_opt.CostUnit;
%             end
            %             [this.index.Router; this.index.GBS]
            this@CloudNetwork(node_opt, link_opt, VNF_opt, options);
            this.base_stations = [net.index.NBS'; dc.gateway']; % the node should be naturally sorted
            this.Topology.Nodes.BaseStation = zeros(this.NumberNodes,1);
            this.Topology.Nodes{this.base_stations, 'BaseStation'} = ...
                (1:length(this.base_stations))';        
%             this.BaseStations = table([net.index.NBS'; dc.gateway'], ...
%                 'VariableName', {'NodeIndex'});
        end
    end
    
    methods 
        function n = get.NumberBaseStations (this)
            n = length(this.base_stations);
        end
        function value = getNodeField(this, name, node_id)
            if strcmp(name, 'BaseStation')
                value = this.Topology.Nodes{node_id, 'BaseStation'};
            else
                if nargin == 2
                    value = getNodeField@CloudNetwork(this, name);
                else
                    value = getNodeField@CloudNetwork(this, name, node_id);
                end
            end
        end
    end
    
    methods(Access = protected)
        [flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt )
    end
end

