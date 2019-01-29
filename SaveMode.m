classdef SaveMode < uint64	
	enumeration
		None(0),			% No save operation
		Save(1),			% Save data to file
		Append(2),		% Append data to file
		Plot(4)				% Plot graph
	end
	
	methods
	end
end

