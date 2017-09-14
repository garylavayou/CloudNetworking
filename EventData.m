classdef EventData < event.EventData
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        userdata;
    end
    
    methods
        function this = EventData(userdata)
            if nargin >= 1 && ~isempty(userdata);
                this.userdata = userdata;
            end
        end
    end 
end