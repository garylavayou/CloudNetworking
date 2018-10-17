classdef DynamicAccessSlice < DynamicSlice
    
    properties (Access = {?DynamicAccessSlice,?CloudNetwork})
        VirtualBaseStations;
    end
    properties(Dependent)
        NumberBaseStations;
    end
    
    methods
        function this = DynamicAccessSlice(slice_data)
            this@DynamicSlice(slice_data);
            %%%
            % *Virtual Base Stations*
            % Select the base station nodes from all the virtual nodes of this slice.
            bs_virtual_node_index = ...
                find(this.Parent.readNode('BaseStation', this.VirtualNodes.PhysicalNode)>0);
            this.VirtualBaseStations = ...
                table(bs_virtual_node_index, 'VariableNames', {'VirtualNode'});
            this.VirtualNodes.BaseStation = zeros(this.NumberVirtualNodes,1);
            this.VirtualNodes{bs_virtual_node_index, 'BaseStation'} = ...
                (1:this.NumberBaseStations)';
        end
        function n = get.NumberBaseStations(this)
            n = height(this.VirtualBaseStations);
        end
        %%%
        % Get physical node index of base stations
        function bs_node_id = getBSNI(this, bs_index)
            if nargin == 1
                bs_index = 1:this.NumberBaseStations;
            end
            vn_id = this.VirtualBaseStations.VirtualNode(bs_index);
            bs_node_id = this.VirtualNodes.PhysicalNode(vn_id);
        end
        %%%
        % Get the physical index of base station.
        function bs_phy_index = getBSPI(this, bs_index)
            if nargin == 1
                bs_node_id = this.getBSNI;
            else
                bs_node_id = this.getBSNI(bs_index);
            end
            bs_phy_index = this.Parent.readNode('BaseStation', bs_node_id);
        end
    end
    
end

