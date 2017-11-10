%% Global State
% Backup and restore the global variables.
classdef GlobalState < matlab.mixin.Copyable
    
    properties (Access = private)
        builder_id = 0;
        entity_id = 0;
        slice_id = 0;
        event_id = 0;
        rand_state = rng;
    end
    
    methods
        function Save(this)
            global builder_id entity_id slice_id eid;
            if ~isempty(builder_id)
                this.builder_id = builder_id;
            end
            if ~isempty(entity_id)
                this.entity_id = entity_id;
            end
            if ~isempty(slice_id)
                this.slice_id = slice_id;
            end
            if ~isempty(eid)
                this.event_id = eid;
            end
            this.rand_state = rng;
        end
        
        function Restore(this)
            global builder_id entity_id slice_id eid;
            builder_id = this.builder_id;
            entity_id = this.entity_id;
            slice_id = this.slice_id;
            eid = this.event_id;
            rng(this.rand_state);
        end
        
        function Display(this)
            for i=1:length(this)
            fprintf('Builder ID: %d, Entity ID: %d, Event ID: %d, Slice ID: %d, Random Seed: %d.\n',...
                this(i).builder_id, this(i).entity_id, this(i).event_id, this(i).slice_id, this(i).rand_state.Seed);
            end
        end
    end
    
    methods (Static)
        function Initialize()
            global builder_id entity_id slice_id eid;
            builder_id = 0;
            entity_id = 0;
            slice_id = 0;
            eid = 0;
        end
    end
end

