%% Physical Network
% network resource abstraction and allocation description; No resource allocation method
% defined, which can be realized by subclasses, see also <CloudNetwork>.
classdef PhysicalNetwork < INetwork
	properties
		DataCenters = table;         % the forwarding node mapping of data center.
		% Inclduing fields: _Capacity_, _StaticCost_, _Load_, _Price_.
		VNFTable = table;    % Meta data of virtual network function
		%%%
		% * *VNFTable*: including following fields.
		%
		% _ProcessEfficiency_: The coefficient for converting service rate to processing resource
		%  requirement, i.e. $ProcessLoad = ServiceRate \times ProcessEfficiency$;
		slices;              % a list of <Slice> objects
		AggregateLinkUsage;
		AggregateNodeUsage;
		slice_template;
	end
	
	%% Properties
	properties (Dependent)
		NumberSlices;        % Number of network slices
		NumberDataCenters;   % Number of data centers
		NumberFlows;         % Number of flows of all slices
		NumberPaths;         % Number of paths in all slices
		NumberVNFs;          % Number of VNF types
		
		LinkOptions;         % Options for link properties
		NodeOptions;         % Options for node properties
	end
	
	properties (Access = protected)
		flow_identifier_generator;
		path_identifier_generator;
  end
	
	methods (Abstract, Access = protected)
		%%%
		% Called by <AddSlice>.
		sl = createslice(this, slice_opt, varargin);
  end  
  
  %% Constructor and Destructor
  %   PhysicalNetwork(node_opt, link_opt, VNF_opt, options)
  % TODO: alternativelt input arguments by struct.
	methods
    function this = PhysicalNetwork(varargin)
      if nargin >= 2
        node_opt = varargin{1};
        link_opt = varargin{2};
        netdata = PhysicalNetwork.loadNetworkData(node_opt, link_opt);
        args = {netdata};
      else
        args = cell(0);
      end
      this@INetwork(args{:});
			if isempty(varargin)
				return;
			end
      
      VNF_opt = varargin{3};
      this.Nodes = netdata.Nodes;
      this.Links = netdata.Edges;
      dc_node_index = find(this.Nodes.Capacity>0);
      this.DataCenters = table(dc_node_index, 'VariableName', {'Node'});
      c = 1;
      while c<=width(this.Nodes)
        name = this.Nodes.Properties.VariableNames{c};
        switch name
          case {'Name', 'Location'}
            c = c + 1;
          otherwise
            this.DataCenters{:,{name}} = this.Nodes{dc_node_index, {name}};
            this.Nodes(:,{name}) = [];
        end
      end
      this.Nodes.DataCenter = zeros(this.NumberNodes,1);
      this.Nodes{dc_node_index, 'DataCenter'} = (1:height(this.DataCenters))';
      this.Links{:,{'Load','Price'}} = 0;
      this.DataCenters{:,{'Load', 'Price'}} = 0;
      
      % Initialize VNF Specification
      this.initializeVNF(VNF_opt);
      
			this.slices = Slice.empty(0,1);
      
			this.flow_identifier_generator = SerialNumber(1, [], true);
      this.path_identifier_generator = SerialNumber(1, [], true);
      
      this.AggregateLinkUsage = zeros(this.NumberLinks,1);
      this.AggregateNodeUsage = zeros(this.NumberNodes,1);
      
			this.options = setdefault(this.options, ...
				struct(...
				'SlicingMethod', SlicingMethod.AdjustPricing,...
				'PricingPolicy', 'linear'));
			if nargin >= 4
				new_opts = varargin{4};
				this.options = structupdate(this.options, new_opts);
			else
				new_opts = struct;
			end
      
      if this.options.SlicingMethod.IsFactorPricing
        % specified for _pricingFactorAjustment_ and .....
        % the pricing factor can be specified by <CloudNetwork.SetOptions>.
        this.options = structmerge(this.options, ...
          getstructfields(new_opts, 'PricingFactor'));
      end
      if this.options.SlicingMethod == SlicingMethod.StaticPricing
        this.options = structmerge(this.options, ...
          getstructfields(new_opts, 'AdmitPolicy'));
      end
    end
		
    function delete(this)
      for i = 1:length(this.slices)
        delete(this.slices(i));
      end
      delete@INetwork(this);
    end
  end

  %% Property Access Methods
  methods
    function opt = get.LinkOptions(this)
      opt = this.Links.Properties.UserData;
    end
    
    function opt = get.NodeOptions(this)
      opt = this.Nodes.Properties.UserData;
    end
    
    function n = get.NumberSlices(this)
      n = length(this.slices);
    end
    
    function n = get.NumberDataCenters(this)
      n = height(this.DataCenters);
    end
    
    function n = get.NumberVNFs(this)
      n = height(this.VNFTable);
    end
    
    function n = get.NumberFlows(this)
      n = 0;
      for i = 1:this.NumberSlices
        n = n + this.slices(i).NumberFlows;
      end
    end
    
    function n = get.NumberPaths(this)
      n = 0;
      for i = 1:this.NumberSlices
        n = n + this.slices(i).NumberPaths;
      end
    end
    
  end
  
  %% Public Methods
  methods
		[output, runtime] = optimizeResourcePrice(this, sub_slices, options);
		
    %% CountSlices
    % statistics of slices.
    function type_count = CountSlices(this)
      type_list = [this.slice_template.Type];
      type_count = zeros(1,length(type_list));
      type = zeros(1,this.NumberSlices);
      for i = 1:this.NumberSlices
        type(i) = this.slices(i).Type;
      end
      [~, tid] = ismember(type, type_list);
      for t = tid
        assert(t~=0, 'error: unknown slice type.');
        type_count(t) = type_count(t) + 1;
      end
    end
    
    %% LinkId
    % link index of links in the edge table
    %
    %      [head, tail] = LinkId(this, idx);
    %      idx = LinkId(this, s, t)
    %
    % |s|: source nodes of links;
    % |t|: tail nodes of links;
    % |idx|: the numeric index of links.
    %
    % See Also _DirectedGraph.IndexEdge_.
    function [argout_1, argout_2] = LinkId(this, argin_1, argin_2)
      switch nargin
        case 1
          argout_1 = this.Links.EndNodes(:,1);
          argout_2 = this.Links.EndNodes(:,2);
        case 2
          argout_1 = this.Links.EndNodes(argin_1,1);
          argout_2 = this.Links.EndNodes(argin_1,2);
        otherwise
          argout_1 = this.graph.IndexEdge(argin_1,argin_2);
      end
    end
    
    %% writeLink
    % set the value for a field in the Edge Table.
    %
    %      writeLink(this, name, value)
    %
    % |name|: a character array represents the field name.
    % |value|: a column vector stores values to be set for the target field.
    function writeLink(this, name, value)
      this.Links{:,{name}} = value;
    end
    
    %% readLink
    % get the value from a field in the Edge Table. See also
    % writeLink.
    %
    %      value = readLink(this, name, link_id)
    %
    %    link field corresponds to column indexed links.
    % |link_id| can be integers or logical numbers.
    function value = readLink(this, name, link_id)
      if nargin < 2 || isempty(name)
        error('input arguments are not enough (name, link_id).');
      end
      if nargin <= 2
        link_id = 1:this.NumberLinks;
      end
      switch name
        case 'ResidualCapacity'
          value = this.Links{link_id,{'Capacity'}} - this.Links{link_id,{'Load'}};
        otherwise
          value = readLink@INetwork(this, name, link_id);
      end
    end
    
    function writeNode(this, name, value, node_index)
      if nargin < 3 || isempty(name) || isempty(value)
        error('input arguments are not enough (name, value).');
      end
      if nargin == 3
        node_index = 1:this.NumberNodes;
      end
      if strcmp(name, 'DataCenter')
        warning('does not support change data center configurations.');
      elseif contains(name, this.Nodes.Properties.VariableNames)
        this.Nodes{node_index,{name}} = value;
      end
    end
    
    function writeDataCenter(this, name, value, dc_index)
      if nargin < 3 || isempty(name) || isempty(value)
        error('input arguments are not enough (name, value).');
      end
      if nargin == 3
        dc_index = 1:this.NumberDataCenters;
      end
      if strcmp(name, 'ResidualCapacity')
        error('error: ResidualCapacity cannot be set.');
      else
        this.DataCenters{dc_index,name} = value;
      end
    end
    
    % * *node_id*: physical node index.
    function value = readNode(this, name, node_id)
      if nargin < 2 || isempty(name)
        error('input arguments are not enough (name, node_id).');
      end
      if nargin < 3
        node_id = 1:this.NumberNodes;
      end
      if contains(name, this.Nodes.Properties.VariableNames)
        value = readNode@INetwork(this, name, node_id);
      elseif contains(name, this.DataCenters.Properties.VariableNames)
        value = zeros(this.NumberNodes, 1);
        value(this.DataCenters.Node) = this.readDataCenter(name);
        value = value(node_id);
        cprintf('SystemCommands', '[%s] read information from data center nodes.\n', calledby);
      else
        error('error:[%s] ''%s'' does not exist.', calledby, name);
      end
    end
    
    % * *dc_id*: data center index.
    function value = readDataCenter(this, name, dc_id)
      if nargin < 2 || isempty(name)
        error('input arguments are not enough (name, node_id).');
      end
      if nargin < 3
        dc_id = 1:this.NumberDataCenters;
      end
      if strcmp('ResidualCapacity', name)
        value = this.DataCenters{dc_id,{'Capacity'}} - ...
          this.DataCenters{dc_id,{'Load'}};
      else
        value = this.DataCenters{dc_id, {name}};
      end
    end

    %% AddSlice
    % Add Slice to Substrate Network
    %
    %      AddSlice(phy_network, slice_opt)
    %
    % |slice_opt|:  option for the added slice;
    sl = AddSlice(this, slice_opt, varargin);

    %% RemoveSlice:
    % Remove the slice with identifier |id|.
    % if no slice with identifier |id|, this method do not perform any operation.
    % |b_update| should be 1, if using static slicing method.
    function sl = RemoveSlice(this, arg1)
      if isnumeric(arg1)
        id = arg1;
        sid = this.findSlice(id, 'Identifier');
      elseif isa(arg1, 'Slice')
        sl = arg1;
        sid = this.findSlice(sl);
      end
      %             if nargin <= 2
      %                 b_update = false;
      %             end
      if ~isempty(sid)
        sl = this.slices(sid);
        link_id = sl.Links.PhysicalLink;
        this.AggregateLinkUsage(link_id) = this.AggregateLinkUsage(link_id) - 1;
        node_id = sl.Nodes.PhysicalNode;
        this.AggregateNodeUsage(node_id) = this.AggregateNodeUsage(node_id) - 1;
        %                 if b_update
        %                     %%%
        %                     % Update the load of the substrate network.
        %                     % for calculate the residual capacity.
        %                     node_load = zeros(this.NumberNodes, 1);
        %                     node_load(sl.Nodes.PhysicalNode) = sl.Nodes.Load;
        %                     this.writeDataCenter('Load', node_load - node_load);
        %                     link_load = zeros(this.NumberLinks, 1);
        %                     link_load(sl.Links.PhysicalLink) = sl.Links.Load;
        %                     this.writeLink('Load', link_load - link_load);
        %                 end
        this.slices(sid) = [];
        %
        % Since the flow/path id might be used by other associate entities, the
        % cost to reallocate flow id is large. And reallocation is not necessary,
        % since the space of identifier is large enough.
        %    this.AllocateFlowId(sid);
        %    this.AllocatePathId(sid);
      else
        sl = [];
      end
    end
    
    %% findSlice
    %   sid = findSlice(this, slice)
    %   sid = findSlice(this, key, field)
    % (1) Find the given slice's Index.
    % (2) Find the index of slices with Type |key|, return a row vector of index.
    function sid = findSlice(this, key, field)
      sid = [];
      if isa(key, 'Slice')
        for s = 1:this.NumberSlices
          if key == this.slices(s)
            sid = s;
            return;
          end
        end
      else
        if nargin <= 2
          field = 'Type';
        end
        switch field
          case 'Type'
            sid = false(1,this.NumberSlices);
            for s = 1:this.NumberSlices
              if key == this.slices(s).Type
                sid(s) = true;
              end
            end
            sid = find(sid);
          case 'Identifier'
            for s = 1:this.NumberSlices
              if key == this.slices(s).Identifier
                sid = s;
                break;
              end
            end
          otherwise
            warning('undefined key type [%s].', key);
        end
      end
    end
    
    %% Modify Options by User
    function setOptions(this, opt_name, opt_value)
      if ischar(opt_name)
        % string type can be indexed by {} operator, returning char array.
        opt_name = string(opt_name);
      elseif iscell(opt_name) || isstring(opt_name)
      else
        error('[%s]error: %s', calledby(0), ...
          '''opt_name'' must be specified as character array, string array',...
          'or cell array with characters');
      end
      
      if isnumeric(opt_value)
        opt_value = num2cell(opt_value);
      end
      
      this.update_options(opt_name, opt_value);
    end
    
    %% plot
    % Visualize Substrate Network and Network Slices
    %
    %      plot(phy_network)
    plot(this, b_undirect);
    
    function V = totalNodeCapacity(this)
      V = sum(this.readDataCenter('Capacity'));
    end
    
    function C = totalLinkCapacity(this)
      C = sum(this.readLink('Capacity'));
    end

    function [r_mean, r_max, r_min, r_std] = nodeUtilization(this)
      node_load = this.readDataCenter('Load');
      node_capacity = this.readDataCenter('Capacity');
      % 			node_index = node_load > 1;
      %%%
      % Another method: ratio = sum(node_load)/sum(node_capacity);
      ratio = node_load ./ node_capacity;
      r_mean = mean(ratio);
      %%%
      % The range of node utilization may large, since the load of nodes depends on
      % the flow's location, the node's cost, and our objective is not to balancing
      % the node load.
      if nargout >= 2
        r_max = max(ratio);
      end
      if nargout >= 3
        r_min = min(ratio);
      end
      if nargout >= 4
        r_std = std(ratio);
      end
    end
    
    function [r_mean, r_max, r_min, r_std] = linkUtilization(this)
      link_load = this.readLink('Load');
      link_capacity = this.readLink('Capacity');
      % link_index = link_load > 1;
      ratio = link_load ./ link_capacity;
      r_mean = mean(ratio);
      if nargout >= 2
        r_max = max(ratio);
      end
      if nargout >= 3
        r_min = min(ratio);
      end
      if nargout >= 4
        r_std = std(ratio);
      end
    end
    
    function theta = utilizationRatio(this, load)
      if nargin == 1
        load.Node = this.readDataCenter('Load');
        load.Link = this.readLink('Load');
      end
      theta_v = sum(load.Node)/this.totalNodeCapacity;
      theta_l = sum(load.Link)/this.totalLinkCapacity;
      theta = 0.5*(theta_v + theta_l);
    end
    
    %%%
    % * *Network Operation Cost*:
    % There are two methods to calculate network cost,
    %
    % # Calculate with the approximate model, where the static node cost is computed
    % by the approximate formula.
    % # Calculate with the accurate model, where the static node cost is computed by
    % the solution of VNF deployment.
    %
    % When the network only include a single slice, this method equals to
    % _getSliceCost_ .
    %         function c = totalCost(this, nload, model)
    function c = totalCost(this, load)
      if nargin <=1 || isempty(load)
        load.Node = this.readDataCenter('Load');
        load.Link = this.readLink('Load');
      end
      
      c = this.totalNodeCost(load.Node) + this.totalLinkCost(load.Link);
    end
    
    %%% compute link cost. Subclass may override this to provide cost.
    function link_uc = getLinkCost(this, link_id)
      if nargin == 1
        link_uc = this.readLink('UnitCost');
      else
        link_uc = this.readLink('UnitCost', link_id);
      end
    end
    
    %%% compute node cost. Subclass may override this to provide cost.
    % * *dc_id*: data center index (not the node index of the substrate physical node).
    function node_uc = getNodeCost(this, dc_id)
      if nargin == 1
        node_uc = this.readDataCenter('UnitCost');
      else
        node_uc = this.readDataCenter('UnitCost', dc_id);
      end
		end
        
    %% statistics of the output
    % type_index is a scalar.
    % Profit, Rate
    % TODO: formulate the output as a table
    function [p,r] = statSlice(this, type_index, profit)
      s_index = this.findSlice(type_index);
      if isempty(s_index)
        p = [0, 0, 0, 0];
        if nargout >= 2
          r = [0, 0, 0, 0];
        end
      else
        %%%
        % Only the statistics of the admitted slices are counted.
        p = [mean(profit(s_index)), max(profit(s_index)), min(profit(s_index)), ...
          std(profit(s_index))];
        if nargout >= 2
          rate = zeros(this.NumberFlows,1);
          total_num_flow = 0;
          for s = s_index     % s_index is a row vector
            num_flow = this.slices(s).NumberFlows;
            rate(total_num_flow + (1:num_flow)) = this.slices(s).FlowTable.Rate;
            total_num_flow = total_num_flow + num_flow;
          end
          rate = rate(1:total_num_flow);
          r = [mean(rate), max(rate), min(rate), std(rate)];
        end
      end
    end

	end
  
	%% Friend Methods
	methods (Access ={?PhysicalNetwork, ?NetworkOptimizer})
		    %% getNetworkLoad
    % Network load equals to the sums of occupied capacity from all slices. See also
    % <Slice.getLinkCapacity>, <Slice.getNodeCapacity> and
    % <DynamicSlice.getLinkCapacity>, <DynamicSlice.getNodeCapacity>.
    %
    % |slices|: If the argument is provided, calculate load from the
    %						set of |slices|, otherwise, calculate from all slices.
    % |options|
    %			_Stage_
    %					'final': directly copy from 'Capacity' field of each slice.
    %					'temp': use the temporary variables of each slice, to
    %							calculate temporary capacity of each slice. Makesure
    %							the temporary variables is up-to-date, when calling
    %							this method.
    function load = getNetworkLoad(this, slices, options)
			defaultopts = struct('Stage', 'final');
			if nargin <= 2
				options = defaultopts;
			else
				options = structupdate(defaultopts, options);
			end
			
			if nargin <= 1 || isempty(slices)
				slices = this.slices;
			end
			
			load.Node = zeros(this.NumberDataCenters, 1);
			load.Link = zeros(this.NumberLinks, 1);
			for i = 1:length(slices)
				sl = slices(i);
				link_id = sl.Links.PhysicalLink;
				dc_id = sl.getDCPI;
				if sl.isFinal() || strcmpi(options.Stage, 'final')
					load.Node(dc_id) = load.Node(dc_id) + sl.getNodeCapacity;
					load.Link(link_id) = load.Link(link_id) + sl.getLinkCapacity;
				elseif strcmpi(options.Stage, 'temp')
					load.Node(dc_id) = load.Node(dc_id) + sl.getNodeCapacity(false);
					load.Link(link_id) = load.Link(link_id)+ sl.getLinkCapacity(false);
				else
					error('error:[%s] Invalid value for options ''Stage''=''%s''.', options.Stage);
				end
			end
		end

    %% Price Adjustment Methods
    [prices, runtime] = pricingFactorAdjustment(this, new_opts);
	end
	
  %% Protected Methods
	methods (Access = protected)
		[sp_profit, b_violate, violates] = SolveSCPCC(this, slices, prices, options);
		
		argout = calculateOutput(this, argin, new_opts);

		function newobj = copyElement(this)
			%% Make a shallow copy of all properties
			newobj = copyElement@INetwork(this);
			%% Deep Copy
			% *slice.Parent*: update this link to the copyed network object.
			for i = 1:this.NumberSlices
				newobj.slices(i) = this.slices(i).copy;
				newobj.slices(i).Parent = newobj;
			end
		end
    
		% calledby <setOptions>.
    function update_options(this, opt_name, opt_value)
      for i = 1:length(opt_name)
        if contains(opt_name{i}, {'SlicingMethod', 'PricingFactor'})
          this.options.(opt_name{i})= opt_value{i};
        end
      end
    end
    
    function [cap, load] = readCapacityLoad(this)
      cap.node = this.readDataCenter('Capacity');
      cap.link = this.readLink('Capacity');
      load.node = this.readDataCenter('Load');
      load.link = this.readLink('Load');
		end
    
    %%% compute the total link cost.
    function c = totalLinkCost(this, link_load)
      if nargin == 1
        c = dot(this.readLink('Load'), this.getLinkCost);
      else
        c = dot(link_load, this.getLinkCost);
      end
    end
    
    %%% compute the total node cost.
    function c = totalNodeCost(this, node_load)
      if nargin == 1
        c = dot(this.DataCenters.Load, this.getNodeCost);
      else
        c = dot(node_load, this.getNodeCost);
      end
    end
    
		%%%
		% Called by _Constructor_;
		function initializeVNF(this, VNF_opt)
			if isfield(VNF_opt, 'ProcessEfficiency')
				assert(VNF_opt.Number==length(VNF_opt.ProcessEfficiency), ...
					'error: process efficiency data is missing.');
				this.VNFTable.ProcessEfficiency = VNF_opt.ProcessEfficiency;
			else
				if isfield(VNF_opt, 'RandomSeed')
					rng(VNF_opt.RandomSeed(1));
				else
					rng(floor(now));
					warning('random seed is not sepecifed for VNF, set as %d', floor(now));
				end
				this.VNFTable.ProcessEfficiency = 0.5 + rand([VNF_opt.Number, 1]);
			end
    end
    
    %%%
		% Called by <AddSlice>.
    function slice_opt = preAddingSlice(this, slice_opt)
      %% Random Seed
      if ~isfield(slice_opt, 'RandomSeed') || isempty(slice_opt.RandomSeed)
        slice_opt.RandomSeed = floor(now);
        warning('random number seed is not specified (set as %d).', slice_opt.RandomSeed);
        % this value should be the same on the same day.
      end
      rng(slice_opt.RandomSeed);
      %% Weight
			getstructfields(slice_opt, 'Weight', 'error');
      
      %% VNF List
      if ~isfield(slice_opt, 'VNFList') || isempty(slice_opt.VNFList)
        if ~isfield(slice_opt, 'NumberVNFs') || isempty(slice_opt.NumberVNFs) ||...
            slice_opt.NumberVNFs > this.NumberVNFs
          error('error: invalid VNF list options.'); % slice_opt.VNFList = unique_randi(this.NumberVNFs, randi([1 4]));
        else
          slice_opt.VNFList = unique_randi(this.NumberVNFs, slice_opt.NumberVNFs);
        end
      else
        if max(slice_opt.VNFList) > this.NumberVNFs
          error('error: invalid VNF list options.');
        end
        % TODO check VNF list
        if isfield(slice_opt, 'NumberVNFs') && ~isempty(slice_opt.NumberVNFs)
          slice_opt.VNFList = slice_opt.VNFList(1:slice_opt.NumberVNFs);
          %     else
          %         slice_opt.VNFList = slice_opt.VNFList;
        end
      end
      
			defaultopts = struct(...
				'FlowPattern', FlowPattern.RandomMultiFlow, ...
				'DuplicateFlow',false, ...
				'DelayConstraint', inf, ...
				'NumberPaths', 1,...
				'SlicingMethod', this.options.SlicingMethod, ...
				'PricingPolicy', this.options.PricingPolicy);
      slice_opt = setdefault(slice_opt, defaultopts, 'warning');
      
      %% Path Constraints
      slice_opt.DelayModel = this.LinkOptions.DelayModel;
      if this.NumberDataCenters < this.NumberNodes
        % if only part of the forwarding nodes is VNF-capable, we should make sure that the
        % path at least transit one VNF-capable node.
        % no matter when, the DataCenters is the middle nodes. However, if the MiddleNodes
        % option is not provided, the route calculation will be performed in a different way.
        slice_opt.MiddleNodes = this.DataCenters.Node;
      end
      
      %% Slicing Method
      % override the slice's specification.
      % By default, the options including 'SlicingMethod' and 'AdmitPolicy' is inherited
      % from the network. But slice can use its own options in the configuration file.
      %% Pricing policy
			% each slice can specify their own pricing, but the network determines whether
			% to adopt this policy or use the network specified pricing policy.
			% (currently, we assume that network's setting override the slice setting.)
			if this.options.SlicingMethod.IsStatic
				slice_opt.AdmitPolicy = this.options.AdmitPolicy;
			end
		
			%% Demand Information
			if slice_opt.FlowPattern ~= FlowPattern.RandomSingleFlow && ...
					~isfield(slice_opt, 'NumberFlows')
				error('[%s] error: <%s> must be specified with <NumberFlows>.', calledby(0), FlowPattern.RandomSingleFlow.char);
			end
			switch slice_opt.FlowPattern
				case FlowPattern.RandomSingleFlow
					slice_opt.NumberFlows = 1;
				case FlowPattern.RandomMultiFlow
					if ~slice_opt.DuplicateFlow
						slice_opt.NumberFlows = min(slice_opt.NumberFlows, this.NumberNodes*(this.NumberNodes-1));
					end
				case FlowPattern.RandomInterDataCenter
					if ~slice_opt.DuplicateFlow
						slice_opt.NumberFlows = min(slice_opt.NumberFlows, this.NumberDataCenters*(this.NumberDataCenters-1));
					end
				otherwise
			end
			switch slice_opt.FlowPattern
				case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
					if ~isfield(slice_opt, 'NodeSet')
						slice_opt.NodeSet = 1:this.NumberNodes;
					end
				case FlowPattern.RandomInterDataCenter
					if ~isfield(slice_opt, 'NodeSet')
						slice_opt.NodeSet = this.DataCenters.Node;
					end
			end
		end
		
		%%%
		% Called by _AddSlice_, return value is used by _generateFlowTable_;
		function graph = residualgraph(this, slice_opt)
      if slice_opt.SlicingMethod.IsStatic
        % If a link's residual capacity is zero, then this link should be removed
        % from the grpah.
        % If a node's residual capacity is zero, then this node and the adjacent
        % links should be removed from the graph.
        link_capacity = this.readLink('ResidualCapacity');
        node_capacity = this.readDataCenter('ResidualCapacity');
        link_capacity(link_capacity<1) = 0;
        node_capacity(node_capacity<1) = 0;
        A = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
        C = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
        for i = 1:this.NumberLinks
          if link_capacity(i) <= 0
            continue;
          end
          h = this.graph.Head(i);
          dc_h = this.Nodes.DataCenter(h);
          if dc_h && node_capacity(dc_h) <= 0
            continue;
          end
          t = this.graph.Tail(i);
          dc_t = this.Nodes.DataCenter(t);
          if dc_t && node_capacity(dc_t) <= 0
            continue;
          end
          A(h, t) = this.graph.Adjacent(h, t); %#ok<SPRIX>
          C(h, t) = link_capacity(i); %#ok<SPRIX>
        end
        graph = DirectedGraph(A, C);
      else
        graph = this.graph;
      end
    end
    
    %%%
		% Called by <AddSlice>
		[flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt );

    %%%
    % Called by <generateFlowTable>
    function end_points = generateEndPoints(this, slice_opt) %#ok<INUSL>
      end_points = slice_opt.NodeSet(unique_randi(numel(slice_opt.NodeSet), 2, 'stable'));
    end
	
    %%%
    % Calculate the profit ratio of slices and network.
    % No matter whether all slices are reconfigured, the profit ratio of all slices
    % will be checked to ensure the profit of ratio of network higher than the
    % threshold.
    %
    % |options|: 'PricingPolicy','Epsilon'.
    %
    % NOTE: _checkProfitRatio_ is a stop condition, which is not directly related to
    % the optimization problem.
    function [b, profit_gap] = checkProfitRatio(this, prices, options)
      global DEBUG;
      
      slice_profit_ratio = zeros(this.NumberSlices,1);
			for s = 1:this.NumberSlices
				sl = this.slices(s);
				revenue = sl.getRevenue();        % get utility
				% Prices announced to each slice.
				if ~sl.isFinal()
					sl.Optimizer.setProblem('Price', prices);
				end
				slice_profit_ratio(s) = sl.op.getProfit(options)/revenue;
			end
      clear revenue;
      
      [sp_profit, sp_revenue] = this.getSliceProviderProfit([], prices, options);
      network_profit_ratio = sp_profit/sp_revenue;
      % a = 1;        % {0.5|0.75|1}
      switch this.op.options.Threshold
        case 'min'
          profit_threshold = min(slice_profit_ratio);
        case 'average'
          profit_threshold = mean(slice_profit_ratio);
        case 'max'
          profit_threshold = max(slice_profit_ratio);
        otherwise
          error('error: invalid option (Threshold = %s)', this.op.options.Threshold);
      end
      if nargin == 4  && isfield(options, 'Epsilon')
        if abs(network_profit_ratio - profit_threshold) < options.Epsilon
          b = true;
        else
          b = false;
        end
      else
        if network_profit_ratio >= profit_threshold
          b = true;
        else
          b = false;
        end
      end
      if nargout == 2
        profit_gap = network_profit_ratio - profit_threshold;
      end
      if ~isempty(DEBUG) && DEBUG
        disp('Profit ratio {slices|network}:');
        disp([slice_profit_ratio; network_profit_ratio]);
      end
    end
    
    %% ISSUE: VNFlist is not conmmonly shared.
    function runtime = priceIteration(this, prices, options)
      if nargout == 1
        slice_runtime = 0;
        runtime.Serial = 0;
			end
			warning('TODO: parallel, see <SolveSCPCC>.');
      for s = 1:this.NumberSlices
        sl = this.slices(s);
        sl.Optimizer.setProblem('Price', prices);
        %%%
        % optimize each slice with price and resource constraints.
        if nargout == 1
          tic;
        end
        sl.Optimizer.optimalFlowRate(options);
        if nargout == 1
          t = toc;
          slice_runtime = max(slice_runtime, t);
          runtime.Serial = runtime.Serial + t;
        end
        sl.Optimizer.setProblem('Price', [])
      end
      if nargout == 1
        runtime.Parallel = slice_runtime;
      end
		end
    
		%% getSliceProviderProfit
    % |slices|: if |slices| is provided, only calsulate the revenue and
    %						cost of the specified |slices|.
    % |prices|: if |prices| are not provided, the stored price are used.
    % |options|: |PricingPolicy| must be specified.
    %
    % Reconfiguration cost does not influence the profit of Slice Provider,
    % see also <optimizeResourcePrice>.
    function [profit, revenue] = getSliceProviderProfit(this, slices, prices, options)
      defaultopts = struct(...
        'PricingPolicy', this.options.PricingPolicy, ... % {linear|quadratic}
				'Stage', 'temp');  
      if nargin <= 3
        options = defaultopts;
      else
        options = structupdate(defaultopts, options, fieldnames(defaultopts), 'silent');
      end
      
      if nargin <= 1 || isempty(slices)
        slices = this.slices;
      end
      if nargin <= 2 || isempty(prices)
        prices.Node = this.readDataCenter('Price');
        prices.Link = this.readLink('Price');
      end
      load = this.getNetworkLoad(slices, options);
      revenue = 0;
      switch options.PricingPolicy
        case {'quadratic-price', 'quadratic'}
          for s = 1:length(slices)
            sl = slices(s);
            link_id = sl.Links.PhysicalLink;
            dc_id = sl.getDCPI;
            % To get the revenue of slice provider, we need to how much
            % resource the slices occupy.
            revenue = revenue + ...
              sl.Optimizer.fcnLinkPricing(prices.Link(link_id), sl.getLinkCapacity(false)) + ...
              sl.Optimizer.fcnNodePricing(prices.Node(dc_id), sl.getNodeCapacity(false));
          end
        case 'linear'
          revenue = dot(load.Node, prices.Node) + dot(load.Link, prices.Link);
        otherwise
          error('%s: invalid pricing policy', calledby);
      end
      profit = revenue - this.totalCost(load);
		end
		
    %%%
    % * *Finalize substrate network*
    %
    % # Record the resource allocation variables, flow rate, virtual node/link load of
    %   each slice.
    % # Virtual Nodes/Links' capacity is derived from node/link load;
    % # Calculate and announce the resource prices to each slice.
    % # Record/update the substrate network's node/link load, price.
    %
    % Usually, this function should be provided with 3 arguments, except that it is
    % called by
    % <file:///E:/workspace/MATLAB/Projects/Documents/CloudNetworking/singleSliceOptimization.html singleSliceOptimization>.
    % NOTE: the price here might be only prcing parameters (for varing pricing
    % policy). To calculate the payment, using _fcnLinkPricing_ and _fcnNodePricing_
    % function.
    function finalize(this, prices, slices)
      if nargin <= 2
        slices = this.slices;
      end
      num_slices = length(slices);
      for i = 1:num_slices
        slices(i).finalize(prices);
      end
      load = this.getNetworkLoad;
      this.writeLink('Load', load.Link);
      this.writeDataCenter('Load', load.Node);
      if nargin >= 3
        % NOTE: prices in the substrate network is updated, while the
        % links/nodes that are not involved in the update procedure, do not
        % change their prices.
        % See also <DynamicCloudNetwork>.<optimizeResourcePrice>.
        pre_link_idx = prices.Link==0;
        prices.Link(pre_link_idx) = this.readLink('Price', pre_link_idx);
        this.writeLink('Price', prices.Link);
        pre_node_idx = prices.Node==0;
        prices.Node(pre_node_idx) = this.readDataCenter('Price', pre_node_idx);
        this.writeDataCenter('Price', prices.Node);
      end
    end
    
  end
	
  %% Static Methods
	methods (Static)
		%%%
		% * *loadNetworkData* |static| : generate graph data.
		%
		%       graph_data = LoadNetworkData(link_opt, node_opt)
		%
		%% TODO: override it for subclass.
		graph_data = loadNetworkData(node_opt, link_opt);
		%%%
		% * *LinkDelay* |static| : convert bandwidth to link delay.
		%
		%       dt = PhysicalNetwork.LinkDelay(delay_opt, bandwidth)
		%
		% |delay_opt|: enumeration type of LinkDelayOption.
		%
		% |bandwidth|: bandwidth with unit of |Mbps|.
		%
		% |dt|: delay with unit of |ms|.
		%
		function dt = LinkDelay(delay_opt, bandwidth)
			switch delay_opt 
				case LinkDelayOption.BandwidthPropotion
					dt = 0.001*bandwidth;
				case delay_opt == LinkDelayOption.BandwidthInverse
					dt = 100./bandwidth;
				otherwise
					error('error: the delay option (%s) cannot be handled.', delay_opt.char);
			end
		end
		
	end
	methods (Static, Access=protected)
		function flag = assert_path_list(end_points, path_list, slice_opt)
			global DEBUG;
			if ~exist('DEBUG', 'var')
				DEBUG = false;
			end
			
			if path_list.Width == 0
				if slice_opt.SlicingMethod.IsStatic
					% two choice: reject the slice or reject the flow.
					if isfield(slice_opt, 'AdmitPolicy') && ...
							strcmp(slice_opt.AdmitPolicy, 'reject-slice')
						message = 'Reject the slice request.';
						flag = 2;
					else
						% slice_opt.AdmitPolicy = 'reject-slice', the actual number
						% of generated flow may less than |number_flow|
						message = sprintf('Reject the flow (%d,%d) in the slice request.',...
							end_points(1), end_points(2));
						flag = 1;
					end
					if DEBUG
						warning(['[', calledby, ']', message]);
					else
						cprintf('SystemCommands', 'Warning: [%s] %s\n', calledby, message);
					end
				else
					flag = -1;
				end
			else
				flag = 0;
			end
		end
  end		
	
	properties (Constant)
		EdgeColor = [ % Predefined set of color data, used to specify the edge color when visualize network.
			0     0     1;    %blue
			0     0.5   0;
			0     0.7   0.7;
			0.078 0.169 0.549;
			0     0.447 0.741;
			0.494 0.184 0.557;
			0.467 0.675 0.188;
			0.302 0.745 0.933;
			1     0.6   0.78;
			];
		NodeColor = [ % Predefined set of color data, used to specify the node color when visualize network.
			1     0     0;  %red
			1     0     1;  %magenta
			0.75  0.75  0;
			0.871 0.5   0;
			0.851 0.325 0.098;
			0.25  0.25  0.25];
	end
	
	methods  % Deprecated
		%%%
		% * *AllocateFlowId* : Allocate flow identifier.
		%
		%      AllocateFlowId(phy_network, start_slice)
		%
		% |start_slice|: the slice index or the handle of the slice. If this argument is
		% not provided, we allocate flow id from the first slice.
		%         function id = AllocateFlowId(this, start_slice)
		%             % Allocate flow identifier
		%             if nargin <= 1
		%                 slice_id = 1;
		%             else
		%                 if isa(start_slice, 'Slice')
		%                     slice_id = this.findSlice(start_slice);
		%                 else
		%                     slice_id = start_slice;
		%                 end
		%             end
		%             if isempty(slice_id) || slice_id == 1
		%                 flow_id = uint64(0);
		%             else
		%                 flow_id = this.slices(slice_id-1).FlowTable.Identifier(end);
		%             end
		%             for s = slice_id:this.NumberSlices
		%                 slice = this.slices(s);
		%                 % Since identifier space is large enough (64-bit), no need to worry that
		%                 % the identifier will duplicate.
		%                 slice.FlowTable.Identifier = flow_id + (1:height(slice.FlowTable))';
		%                 flow_id = flow_id + height(slice.FlowTable);
		%             end
		%         end
		
		%%%
		% * *AllocatePathId* : Allocate path identifier
		%
		%      AllocatePathId(phy_network)
		function AllocatePathId(this, start_slice)
			if nargin <= 1
				slice_id = 1;
			else
				if isa(start_slice, 'Slice')
					slice_id = this.findSlice(start_slice);
				else
					slice_id = start_slice;
				end
			end
			if isempty(slice_id) || slice_id == 1
				path_id = uint64(0);
			else
				path_list = this.slices(slice_id-1).FlowTable{end, 'Paths'};
				path_id = path_list{path_list.Width}.id;		% end cannot be used for {}
			end
			for s = slice_id:this.NumberSlices
				for j = 1:this.slices(s).NumberFlows
					path_list = this.slices(s).FlowTable{j, 'Paths'};
					for k = 1:path_list.Width
						path_id = path_id + 1;
						path_list{k}.id = path_id;
					end
				end
			end
		end
	end
end
