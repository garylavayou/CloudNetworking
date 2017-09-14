%% EventDrivenNetwork 
% This class extend <PhysicalNetwork> to enable handling of network events.
%
% If we do not derive this class from <PhysicalNetwork>, then this class can also be
% treated as an interface, if we used it with subclasses of <PhysicalNetwork>.
classdef DynamicNetwork < PhysicalNetwork & EventSender & EventReceiver

    methods
        % EventDrivenNetwork(node_opt, link_opt, VNF_opt, net_opt)
        function this = DynamicNetwork(varargin)
            this@PhysicalNetwork(varargin{:});
        end
        
        %%%
        % Since CloudNetwork and EventDrivenNetwork has different definition of AddSlice,
        % we need override the two superclass methods.
        function sl = AddSlice(this, slice_opt, varargin)
            if ~cellstrfind(this.bypass, 'DynamicNetwork')
                AddSlice@PhysicalNetwork(this, slice_opt, varargin{:});
            end
            %% Perform admiting control and resource allocation
            % Compute the resource allocation for the slice, and decide if this slice can
            % be admitted.
            % NOTE: currently we admit all slice request, so this function also does not
            % do anything.
            %
            % Subclasses should implement <onAddingSlice>.
            sl = this.onAddingSlice(this.slices{end});
            %% TODO: provide options when constructing the network.
        end        
        function sl = RemoveSlice(this, arg1)
            sl = RemoveSlice@PhysicalNetwork(this, arg1);
            sl = this.onRemovingSlice(sl);
        end
                
    end
    
    methods
        function eventhandler(this, source, eventData)
            global DEBUG; %#ok<NUSED>
            % target = eventData.targets;
            % where target should be the <DynamicNetwork> object.
            %% TODO: adding dynamic slice dimensioning for flow arrival/departure scale.
            ev = eventData.event;
            switch eventData.EventName
                case 'SliceArrive'
                    et = ev.Entity;
                    sl = this.AddSlice(et.Options);
                    % 2. Allocate flow id: should it be allocated by the Network or
                    % allocated by the Event Dispatcher.
                    % Inform the <EventDispatcher> to add existing flows and the flow
                    % entity builder.
                    if isempty(sl)
                        notify(this, 'AddSliceFailed');
                    else
                        sl.Identifier = et.SliceIdentifier;
                        if sl.isDynamicFlow
                            this.AddListener(sl, {'FlowArrive', 'FlowDepart'}, @sl.eventhandler);
                            sl.AddListener(this, {'AddFlowSucceed', 'AddFlowFailed', ...
                                'RemoveFlowSucceed', 'RemoveFlowFailed'}, @this.eventhandler);
                        end
                        data = DispatchEventData(ev, sl);
                        % The |entity| and |slice| information are needed to create flow
                        % entitities.
                        notify(this, 'AddSliceSucceed', data);
                    end
                    % Who should listen to slice and who should slice send events to?
                    % sl.AddListeners()
                case 'SliceDepart'
                    % Remove the slice if 'mandatorydepart' is enable; otherwise, if
                    % 'naturaldepart' is enabled, we set a flag. and until the last one
                    % flow of this slice departs, we remove the slice.
                    sl = ev.userdata;
                    if strcmpi(source.flow_depart_option, 'naturaldepart')
                        sl.b_ondepart = true;
                    else
                        this.RemoveSlice(sl);
                    end
                case 'FlowArrive'
                    %%%
                    % Notify slice that a flow arrives.
                    % |SliceIdentifier| of <SliceEntity> is equal to the |Identifier| of
                    % <Slice>.
                    % NOTE: handle may also be used in place of |Identifier|.
                    % TODO: pass the slice identifier directly through event data.
                    slice_id = ev.Entity.Parent.SliceIdentifier;
                    sl = this.slices{this.findSlice(slice_id, 'Identifier')};
                    ft = this.createflow(sl);
                    data = FlowEventData(ev, sl, ft);
                    %                     data.targets = this.FindSlice(et.Parent.SliceIdentifier);
                    notify(this, 'FlowArrive', data);
                case 'FlowDepart'
                    % notify slice
                    % TODO: remove flow entity child form slice entity
                    slice_id = ev.Entity.Parent.SliceIdentifier;
                    sl = this.slices{this.findSlice(slice_id, 'Identifier')};
                    flow_id = ev.Entity.GlobalIdentifier;
                    data = FlowEventData(ev, sl, flow_id);
                    notify(this, 'FlowDepart', data);
                    if strcmpi(source.flow_depart_option, 'naturaldepart')
                        if sl.b_ondepart && sl.NumberFlows == 0
                            % Finally remove the slice from the network
                            this.RemoveSlice(sl);
                            % then remove the slice entity form the entity list.
                            data = FlowEventData(ev, sl, []);
                            notify(this, 'RemoveSliceSucceed', data);
                        end
                    end
                case 'AddFlowSucceed'
                    % allocate flow id
                    % notify <EventDispatcher> to assign identifier to flow entry
                    fidx = eventData.flow;
                    sl = eventData.slice;
                    identifier = this.flow_identifier_generator.next(length(fidx));
                    sl.FlowTable{fidx, 'Identifier' } = identifier;
                    eventData = DispatchEventData(eventData.event, identifier);
                    notify(this, 'AddFlowSucceed', eventData);
                case 'AddFlowFailed'
                    % notify <EventDispatcher> to remove the invalid flow entry.
                case 'RemoveFlowSucceed'
                    data = EventData(eventData.entity);
                    notify(this, 'RemoveFlowSucceed', data);
                case 'RemoveFlowFailed'
                otherwise
                    error('error: cannot handle event %s.', eventData.EventName);
            end
        end
    end
        
    methods (Abstract, Access = protected)
        %% Perform admitting control and resource allocation
        % Compute the resource allocation for the slice, and decide if this slice can
        % be admitted.
        % NOTE: currently we admit all slice request, so this function also does not
        % do anything.
        %         function tf = OnAddlingSlice(this)
        %             tf = true;
        %         end
        tf = onAddingSlice(this, sl);
        tf = onRemovingSlice(this, sl);
    end
    methods (Access =protected)
        function sl = createslice(this, slice_opt)
            % examine flow arrival parameters.
            % usage of <Slice>: if a slice without |ArrivalRate| or |ServiceInterval| or
            % their values are invalid, the slice (<Slice> or <DynamicSlice>) is treated
            % as no dynamics of flow, and it will not handle flow events. Therefore
            % initilize it as class <Slice> is OK.
            if ~isfield(slice_opt, 'ArrivalRate') || ~isfield(slice_opt, 'ServiceInterval')
                this.slices{end+1} = Slice(slice_opt);
            elseif isempty(slice_opt.ArrivalRate) || isempty(slice_opt.ServiceInterval)
                this.slices{end+1} = Slice(slice_opt);
                warning('slice created with type Slice.');
            else
                this.slices{end+1} = DynamicSlice(slice_opt);
            end
            sl = this.slices{end};
        end
        %%%
        % *Create new flows*
        % Creating new flows in the slice could guarantee no extra node or link would be
        % needed. If we enable new flows from new locations, we should create the flow
        % in the network.
        % |ft|: return flow table entries.
        function ft = createflow(this, slice, numflow)
            % map virtual network to physical network
            A = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
            C = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
            for i = 1:slice.NumberVirtualLinks
                h = slice.Topology.Head(i);
                t = slice.Topology.Tail(i);
                ph = slice.VirtualNodes{h, 'PhysicalNode'};
                pt = slice.VirtualNodes{t, 'PhysicalNode'};
                A(ph, pt) = slice.Topology.Adjacent(h,t); %#ok<SPRIX>
                C(ph, pt) = slice.Topology.Capacity(h,t); %#ok<SPRIX>
            end
            graph = DirectedGraph(A, C);
            slice_opt.Pattern = slice.Options.FlowPattern;
            slice_opt.DelayConstraint = slice.Options.DelayConstraint;
            slice_opt = update_slice_option(slice, slice_opt);
            if nargin <= 2
                slice_opt.NumberFlows = 1;
            else
                slice_opt.NumberFlows = numflow;
            end
            slice_opt.NumberPaths = slice.Options.NumberPaths;
            slice_opt.method = 'dynamic-slicing';
            b_vailid_flow = false;
            while ~b_vailid_flow
                try
                    b_vailid_flow = true;
                    ft = this.generateFlowTable(graph, slice_opt);
                catch ME
                    disp(ME)
                    b_vailid_flow = false;
                end
            end
            %%%
            % Update slice information as creating slice.
            %% TODO
            % When new nodes/edges should be added.
            ft.Properties.VariableNames = ...
                {'Source', 'Target', 'Rate', 'Delay', 'Paths'};
            for k = 1:height(ft)
                path_list = ft{k,{'Paths'}};
                for p = 1:path_list.Width
                    path = path_list.paths{p};
                    path.node_list = slice.PhyscialNodeMap{path.node_list,'VirtualNode'};
                    path.id = this.path_identifier_generator.next;
                end
            end
        end
    end
    events
        FlowArrive;
        FlowDepart;
        AddSliceSucceed;
        AddSliceFailed;
        RemoveSliceSucceed;
        RemoveSliceFailed;          % NOT used.
        AddFlowSucceed;
        AddFlowFailed;
        RemoveFlowSucceed;
        RemoveFlowFailed;          % NOT used.
    end
        
    methods(Static)
        %%%
        % subclass can override this method.
        function slice_opt = update_slice_option(slice, slice_opt)
            switch slice.Options.FlowPattern
                case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
                    slice_opt.NodeSet = slice.VirtualNodes.PhysicalNode;
                otherwise
                    error('error: cannot handle the flow pattern <%s>.', ...
                        slice.Options.FlowPattern.char);
            end
        end
    end
end

