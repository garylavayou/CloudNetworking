classdef Slice < VirtualNetwork	
	properties
		Type;           % Type from slice template
		
		FlowTable;      % flow information in the slice
		VNFList;        % List of virtual network functions in the slice
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
				args = cell(0);
			else
				args = {slice_data};
			end
			this@VirtualNetwork(args{:});

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
	
	methods (Static)
		slice_template = loadSliceTemplate(index);
	end
	
end

