classdef StaticProperties < handle
    
    properties
        props;
    end
    
    methods
        function set(this, field, value)
            this.props.(field) = value;
        end
        
        function v = get(this, field)
            v = this.props.(field);
        end
    end
end

