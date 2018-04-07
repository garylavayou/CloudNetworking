%% Global State
% Backup and restore the global variables.
classdef GlobalState < matlab.mixin.Copyable
    
    properties (Access = private)
        builder_id = 0;
        entity_id = 0;
        slice_id = 0;
        event_id = 0;
        rand_state;
    end
    
    methods
        function this = GlobalState()
            this.rand_state = rng;
        end
        function Save(this)
            this.builder_id = EntityBuilder.getGlobalBuilderId();
            this.entity_id = Entity.getGlobalEntityId();
            this.slice_id = SliceEntity.getGlobalSliceId();
            this.event_id = Event.getGlobalEventId();
            this.rand_state = rng;
        end
        
        function Restore(this)
            EntityBuilder.setGlobalBuilderId(this.builder_id);
            Entity.setGlobalEntityId(this.entity_id);
            SliceEntity.setGlobalSliceId(this.slice_id);
            Event.setGlobalEventId(this.event_id);
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
            EntityBuilder.setGlobalBuilderId(0);
            Entity.setGlobalEntityId(0);
            SliceEntity.setGlobalSliceId(0);
            Event.setGlobalEventId(0);
        end
    end
end

