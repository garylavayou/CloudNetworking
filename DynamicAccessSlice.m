classdef SimpleDynamicAccessSlice < SimpleDynamicSlice & IAccessSlice
    methods
        function this = SimpleDynamicAccessSlice(slice_data)
					this@SimpleDynamicSlice(slice_data);
					this@IAccessSlice(slice_data);
        end
    end
    
end

