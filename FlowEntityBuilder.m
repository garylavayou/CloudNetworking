%% Flow Entity Builder
% In addition to functionality of <EntityBuilder>, *FlowEntityBuilder* can generate the
% local and global identifiers for flow entities.
classdef FlowEntityBuilder < EntityBuilder
    
    properties (Access = private)
        local_flow_identifier = SerialNumber(1, [], true);
    end
    properties (SetAccess = private)
        Parent;
    end
    properties (Dependent)
        % Type is derived from <EnityBuilder>
        Type;
    end
    
    methods
        %%
        % See also <EntityBuilder>.
        % The |arrival_rate| and |service_interval| of flow in a slice is parameters of
        % the slice. 
        % |Parent|: parent slice of the builder aand flow entities produced by this builder.
        function this = FlowEntityBuilder(varargin)
            if isstruct(varargin{1})
                flow_opt = varargin{1};
                arrive_rate = flow_opt.ArrivalRate;
                service_interval = flow_opt.ServiceInterval;
                if nargin >= 2
                    parent = varargin{2};
                else
                    parent = [];
                end
%                 flow_opt = rmfield(flow_opt, {'ArrivalRate', 'ServiceInterval'});
            elseif nargin >= 2
                arrive_rate = varargin{1};
                service_interval = varargin{2};
                if nargin >= 3
                    parent = varargin{3};
                else
                    parent = [];
                end
                flow_opt = [];
            else
                error('error: arguments not enough.');
            end
            this@EntityBuilder(arrive_rate, service_interval, flow_opt);
            this.Parent = parent;
        end
        function t = get.Type(~)
            t = EntityType.Flow;
        end
    end
    
    methods(Access=protected)
        function entity = buildentity(this, time_arrive, time_serve, flow_id)
            % build the entity
            entity = FlowEntity(time_arrive, time_serve, this);
            entity.LocalIdentifier = this.local_flow_identifier.next;
            if nargin >= 4
                % might be unset when construct the entity.
                entity.GlobalIdentifier = flow_id;
            end
        end
    end
end

