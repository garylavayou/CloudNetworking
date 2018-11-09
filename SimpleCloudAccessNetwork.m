%% CloudAccessNetwork
% Incorporate access network into cloud network.
% NOTE: CloudAccessNetwork can be derived from <CloudNetwork> and <AccessNetwork> (TODO),
% which are both derived from <PhysicalNetwork>.
classdef SimpleCloudAccessNetwork < SimpleCloudNetwork & IAccessNetwork
	methods
		function this = SimpleCloudAccessNetwork(node_opt, link_opt, VNF_opt, options)
			%% Initialize as SD-RAN
			[node_opt, link_opt] = IAccessNetwork.loadNetworkData(node_opt, link_opt);
			this@SimpleCloudNetwork(node_opt, link_opt, VNF_opt, options);
			this@IAccessNetwork(node_opt);
		end
	end
	
	methods
		%% readNode
		% override <SimpleCloudNetwork>.readNode and <IAccessNetwork>.readNode.
		function value = readNode(this, name, node_id)
			if nargin <= 2
				args = {name};
			else
				args = {name, node_id};
			end
			value = readNode@IAccessNetwork(this, args{:});
			if isempty(value)
				value = readNode@SimpleCloudNetwork(this, args{:});
			end
		end
	end
	
	methods(Access = protected)		
		function slice_opt = preAddingSlice(this, slice_opt)
			slice_opt = preAddingSlice@SimpleCloudNetwork(this, slice_opt);
			slice_opt = preAddingSlice@IAccessNetwork(this, slice_opt);
		end
		
		function end_points = generateEndPoints(this, slice_opt)
			end_points = generateEndPoints@IAccessNetwork(this, slice_opt);
			if isempty(end_points)
				end_points = generateEndPoints@SimpleCloudNetwork(this, slice_opt);
			end
		end
	end
	
	methods (Static)
		%% Override
		function [node_opt, link_opt] = loadNetworkData(node_opt, link_opt)
			[node_opt, link_opt] = loadNetworkData@IAccessNetwork(node_opt, link_opt);
		end
	end
	
end
