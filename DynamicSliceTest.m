classdef DynamicSliceTest < DynamicSlice
%% 
% Add a warm-up phase for experiment, in the warm-up phase we only use the 'dimconfig' method to
% reconfigure slices. 
    properties (Access = protected)
        saved_options;
    end
    
    methods
        function this = DynamicSliceTest(slice_data)
            this@DynamicSlice(slice_data);
            this.saved_options = this.options;
            this.options = structmerge(this.options, ...
                getstructfields(slice_data, 'NumberEventWarmUp','default-ignore', 50));
            this.options.ReconfigMethod = 'dimconfig';
        end
    end
    
    
    methods(Access=protected)
        function [exitflag,fidx] = executeMethod(this, action)
            global event_num DEBUG;  %#ok<NUSED>
            [exitflag,fidx] = executeMethod@DynamicSlice(this, action);
            if event_num == this.getOption('NumberEventWarmUp')
                this.options.ReconfigMethod = this.saved_options.ReconfigMethod;
            end
        end
    end
    
end

