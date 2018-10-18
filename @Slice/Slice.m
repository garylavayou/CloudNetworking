classdef Slice < VirtualNetwork	
	properties
		Type;           % Type from slice template
		
		FlowTable;      % flow information in the slice
		VNFList;        % List of virtual network functions in the slice
	end
	
	properties
		I_edge_path logical;    % Edge-Path Incidence Matrix
		I_dc_path logical;    % Node-Path Incidence Matrix
		I_flow_path logical;    % Flow-Path Incidence Matrix
		% the associated flow of path, to provide a fast inquiry method for associated
		% flow than using |I_flow_path|.
		path_owner;
	end
	
	properties (Dependent)
		NumberFlows;					% Number of flows in the slice
		NumberVNFs;						% Number of virtual network functions in the slice
	end
	
	methods (Abstract)
		
	end
	
	methods
		function this = Slice(slice_data)
			if nargin == 0
				return;
			end
			this@VirtualNetwork(slice_data);

			if isfield(slice_data, 'Type')
				this.Type = slice_data.Type;
			end
			
			% Flow Table
			% convert node index to virtual node index
			phy_nodes = this.PhysicalNodeMap;
			this.FlowTable = slice_data.FlowTable;
			for k = 1:height(this.FlowTable)
				path_list = this.FlowTable{k,{'Paths'}};
				for p = 1:path_list.Width
					path = path_list.paths{p};
					path.node_list = phy_nodes(path.node_list);
				end
			end
			
			this.VNFList = vnet_data.VNFList;
			
			% Flow options
			this.options = getstructfields(vnet_data, ...
				{'FlowPattern', 'DelayConstraint', 'NumberPaths'});
		end
		
		function delete(this)
			for f = 1:this.NumberFlows
				delete(this.FlowTable{f,'Paths'});
			end
		end
		
	end
	
	methods (Access = protected)
		function newobj = copyElement(this)
			newobj = copyElement@VirtualNetwork(this);
			% Deep Copy
			% [FlowTable.Paths]
			for f = 1:this.NumberFlows
				% 'Paths' is a handle object, is should be copyed to the new table.
				newobj.FlowTable{f,'Paths'} = this.FlowTable{f,'Paths'}.copy;
			end
		end
		
	end
	
	methods 
		function f = get.NumberFlows(this)
			f = height(this.FlowTable);
		end
		
		function n = get.NumberVNFs(this)
			n = length(this.VNFList);
		end
	end
	
	methods
		function initializeState(this)
			NC = this.NumberServiceNodes;
			NP = this.NumberPaths;
			NL = this.NumberLinks;
			NF = this.NumberFlows;
			
			this.I_dc_path = sparse(NC, NP);
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
						dc_index = this.Nodes{e(1),'ServiceNode'};
						if dc_index~=0
							this.I_dc_path(dc_index, pid) = 1;
						end
					end
					dc_index = this.Nodes{e(2),'ServiceNode'}; % last node
					if dc_index~=0
						this.I_dc_path(dc_index, pid) = 1;
					end
				end
			end
		end
		[payment, grad, pseudo_hess] = fcnLinkPricing(this, link_price, link_load);
		[payment, grad, pseudo_hess] = fcnNodePricing(this, node_price, node_load);
	end
	
	methods (Static)
		slice_template = loadSliceTemplate(index);
	end
	
end

