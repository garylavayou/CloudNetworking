%% Virtual Network
% define resource allocation and information in a virtual network.
classdef VirtualNetwork < matlab.mixin.Copyable 
    % Specify the properties that can only be modified by Physcial Network directly
    properties (SetAccess = {?VirtualNetwork,?PhysicalNetwork,?SliceFlowEventDispatcher})
        Identifier;     %
        Parent;
        Type;           % Type from slice template
        Topology;       % topology information of the slice <DirectedGraph>
        VirtualLinks;   % information of virtual links in the slice
        VirtualNodes;   % information of virtual nodes
        VirtualDataCenters; % data centers in the slice 
        FlowTable;      % flow information in the slice
        VNFList;        % List of virtual network functions in the slice
        Options;
        
        I_edge_path;    % Edge-Path Incidence Matrix
        I_node_path;    % Node-Path Incidence Matrix
        I_flow_path;    % Flow-Path Incidence Matrix
        
        
        % I_path_function;        % used by <singleSliceOptimization>.
        %         constant_profit = 0;	% consant profit added to avoid slices with negative profit.
        % => move to SliceEx
    end
    
    properties (GetAccess = {?PhysicalNetwork,?Slice})
        PhysicalLinkMap;
        PhyscialNodeMap;    %
        path_owner;         % the associated flow of path;
        %         I_active_variable;  % indicator of active variables
    end
   
    properties (Dependent)
        NumberVirtualNodes;    % Number of nodes in the slice
        NumberDataCenters;     % Number of virtual data centers.
        NumberVirtualLinks;    % Number of edges in the slice
        NumberFlows;    % Number of flows in the slice
        NumberPaths;    % Number of paths in the slice
        NumberVNFs;     % Number of virtual network functions in the slice
    end
        
    methods
        function this = VirtualNetwork(net_data)
            if nargin == 0
                return;
            elseif isa(net_data, 'VirtualNetwork')
                this = net_data.copy;
                return;
            end
            this.Parent = net_data.parent;
            this.Topology = DirectedGraph(net_data.adjacent);
            if isfield(net_data, 'Type')
                this.Type = net_data.Type;
            end
            %%%
            % Virtual Links
            this.VirtualLinks = array2table(net_data.link_map_S2P, ...
                'VariableNames', {'PhysicalLink'});
            % Link capacity
            if isfield(net_data, 'link_capacity')
                this.VirtualLinks.Capacity = net_data.link_capacity;
            end
            % Link load
            this.VirtualLinks.Load = zeros(this.NumberVirtualLinks, 1);
            this.PhysicalLinkMap = array2table(net_data.link_map_P2S,...
                'VariableNames', {'VirtualLink'});
            %%%
            % Virtual Nodes
            this.VirtualNodes = array2table(net_data.node_map_S2P,...
                'VariableNames', {'PhysicalNode'});
            %%%
            % Virtual Data Center Node
            % Select the data center nodes from all the virtual nodes of this slice.
            dc_virtual_node_index = ...
                find(this.Parent.getNodeField('Capacity', this.VirtualNodes.PhysicalNode) > 0);
            this.VirtualDataCenters = ...
                array2table(dc_virtual_node_index, 'VariableNames', {'VirtualNode'});
            this.VirtualNodes.DataCenter = zeros(this.NumberVirtualNodes,1); 
            this.VirtualNodes{dc_virtual_node_index,'DataCenter'} = ...
                (1:this.NumberDataCenters)';
            %%%
            % Data center node capacity
            if isfield(net_data, 'node_capacity')
                this.VirtualDataCenters.Capacity = net_data.node_capacity;
            else
                this.VirtualDataCenters.Capacity = zeros(this.NumberDataCenters,1);
            end
            %%%
            % Data center node load
            this.VirtualDataCenters.Load = zeros(this.NumberDataCenters, 1);
            this.PhyscialNodeMap = array2table(net_data.node_map_P2S,...
                'VariableNames', {'VirtualNode'});
            %%%
            % Flow Table
            % convert node index to virtual node index
            for k = 1:height(net_data.flow_table)
                path_list = net_data.flow_table{k,{'Paths'}};
                for p = 1:path_list.Width
                    path = path_list.paths{p};
                    path.node_list = this.PhyscialNodeMap{path.node_list, 'VirtualNode'};
                end
            end
            this.FlowTable = net_data.flow_table;
            
            this.VNFList = net_data.VNFList;
            this.initializeState;

            if isfield(net_data,'Identifier')
                this.Identifier = net_data.Identifier;
            end
            
%             if isfield(net_data, 'ConstantProfit')  => move to SliceEx
%                 this.constant_profit = net_data.ConstantProfit;
%             end
            
            if isfield(net_data, 'Pattern')
                this.Options.FlowPattern = net_data.Pattern;
            else
                warning('FlowPattern option is not provided.');
            end
            if isfield(net_data, 'DelayConstraint')
                this.Options.DelayConstraint = net_data.DelayConstraint;
            else
                warning('DelayConstraint option is not provided.');
            end
            if isfield(net_data, 'NumberPaths')
                this.Options.NumberPaths = net_data.NumberPaths;
            else
                warning('NumberPaths option is not provided.');
            end    
        end
        
        function n = get.NumberVirtualNodes(this)
            n = this.Topology.NumberNodes;
        end
        
        function n = get.NumberDataCenters(this)
            n = height(this.VirtualDataCenters);
        end

        function m = get.NumberVirtualLinks(this)
            m = this.Topology.NumberEdges;
        end
        
        function f = get.NumberFlows(this)
            f = height(this.FlowTable);
        end
        
        function p = get.NumberPaths(this)
            p = 0;
            for i=1:this.NumberFlows
                p = p + this.FlowTable.Paths(i).Width;
            end
        end
        
        function n = get.NumberVNFs(this)
            n = length(this.VNFList);
        end
        %         function p = get.Parent(this)
        %             p = this.Parent;
        %         end
        %         function n = get.num_vars(this)
        %             n = (this.NumberVNFs*this.NumberVirtualNodes+1)*this.NumberPaths;
        %         end
        
        %%%
        % Get physical node index of data center
        function dc_node_id = getDCNI(this, dc_index)
            if nargin == 1
                dc_index = 1:this.NumberDataCenters;
            end
            vn_id = this.VirtualDataCenters.VirtualNode(dc_index);
            dc_node_id = this.VirtualNodes.PhysicalNode(vn_id);
        end
        
        %%%
        % Get the physical index of data center.
        function dc_phy_index =getDCPI(this,dc_index)
            if nargin == 1
                dc_node_id = this.getDCNI;
            else
                dc_node_id = this.getDCNI(dc_index);
            end
            dc_phy_index = this.Parent.getNodeField('DataCenter', dc_node_id);
        end
        
        function pid = getLocalPathId(this, path)
            %%%
            % {need override}
            cprintf('comment', '%s%s%s', ...
                'getLocalPathId should be overrided,', ...
                'if subclasses dynamically manage flows. ', ...
                'Instead, using path.local_id, which should be dynamically maintained');
            pid = path.id - this.FlowTable.Paths(1).paths{1}.id + 1;
        end

%         function pl = getPathLength(this)
%             pl = zeros(this.NumberPaths,1);
%             pid = 1;
%             for i = 1:this.NumberFlows
%                 pathlist = this.FlowTable{i,'Paths'}.paths;
%                 for l = 1:length(pathlist)
%                     pl(pid) = pathlist{l}.Length;
%                 end
%             end
%         end
                
        function tf = isDynamicFlow(~)
            tf = false;
        end
        
        function initializeState(this)
            NC = this.NumberDataCenters;
            NP = this.NumberPaths;
            NL = this.NumberVirtualLinks;
            NF = this.NumberFlows;
            
            this.I_node_path = sparse(NC, NP);
            this.I_edge_path = sparse(NL, NP);
            this.I_flow_path = sparse(NF, NP);
            this.path_owner = zeros(NP,1);
            % this.local_path_id = zeros(this.NumberPaths, 1);
            pid = 0;
            for fid=1:NF
                path_list = this.FlowTable{fid,{'Paths'}};
                for j = 1:path_list.Width
                    pid = pid + 1;
                    this.I_flow_path(fid,pid) = 1;
                    this.path_owner(pid) = fid;
                    path = path_list.paths{j};
                    path.local_id = pid;    % record the local path in the slice.
                    for k = 1:(path.Length-1)
                        e = path.Link(k);
                        eid = this.Topology.IndexEdge(e(1),e(2));
                        this.I_edge_path(eid, pid) = 1;
                        dc_index = this.VirtualNodes{e(1),'DataCenter'};
                        if dc_index~=0
                            this.I_node_path(dc_index, pid) = 1;
                        end
                    end
                    dc_index = this.VirtualNodes{e(2),'DataCenter'}; % last node
                    if dc_index~=0
                        this.I_node_path(dc_index, pid) = 1;
                    end
                end
            end
        end
    end
        
    methods (Access = protected)
        %% Deep Copy
        function this = copyElement(sl)
            % Make a shallow copy of all properties
            this = copyElement@matlab.mixin.Copyable(sl);
            % Make a deep copy of the DeepCp object
            this.Topology = sl.Topology.copy;
            for f = 1:height(sl.FlowTable)
                % path_list is handle object, is should be copyed to the new table.
                this.FlowTable{f,'Paths'} = sl.FlowTable{f,'Paths'}.copy;
            end
        end
    end
    
    %     methods (Static)
    %         slice_template = loadSliceTemplate(index);
    %     end
end

%% Properties
% * *Topology*: Normally, the network slice will not run shrtest path algorithm, so the
% absolute value of the adjacent matrix of Topology does not matter. On the other hand,
% the link and node capacity of the slice is also not determined until the substrate
% network allocate the resource to the slice.
% * *VirtualLinks* : fields include _PhysicalLink_, _Price_, _Load_.
% * *VitrualNodeMap* : fields include _PhysicalNode_, _Price_, _Load_.
%
% * *NumberVirtualNodes* |get|
%
% * *NumberVirtualLinks* |get|
%
% * *NumberFlows* |get|
%
% * *NumberPaths* |get|

%% Methods
% * *getLocalPathId* : find the path's local identifier.
%
%      pid = getLocalPathId(slice, path)
%