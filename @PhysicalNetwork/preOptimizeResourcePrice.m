function [slices, capacities, unitcosts, fields, options] = ...
	preOptimizeResourcePrice(this, slices, options)
%% Initialization
Ne = this.NumberLinks;
Ndc = this.NumberDataCenters;
if nargin <= 1 || isempty(slices)
  % all slices are involved in slice dimensioning
  slices = this.slices;
	capacities = [this.readLink('Capacity');	this.readDataCenter('Capacity')];
else
  % residual capacity + reallocatable capacity.
	capacities = [this.readLink('ResidualCapacity');	...
		this.readDataCenter('ResidualCapacity')];
	capacities(capacities<0) = 0;
	for i = 1:length(slices)
    sl = slices(i);
    capacities(sl.Links.PhysicalLink) = ...
			capacities(sl.Links.PhysicalLink) + sl.Links.Capacity;
		dc_id = sl.getDCPI();
    capacities(Ne+dc_id) = capacities(Ne+dc_id) + sl.ServiceNodes.Capacity;
  end
end

if nargin <= 2 || ~isfield(options, 'PricingMethod')
  cost_init_method = 'RandomizeCost';
else
  cost_init_method = options.PricingMethod;
  options = rmfield(options, 'PricingMethod');
end
switch cost_init_method
  case 'RandomizeCost'
    st = rng;
    rng(20180909);
    unitcosts = [this.getLinkCost(); this.getNodeCost()].*((rand(Ne+Ndc,1)-0.5)/100+1);
    rng(st);
  case 'UniformCost'
    st = rng;
		rng(20180909);
    unitcosts = [ones(Ne,1)*mean(this.getLinkCost());
			ones(Ndc,1)*mean(this.getNodeCost())].*((rand(Ne+Ndc,1)-0.5)/100+1);
    rng(st);
  case 'OriginCost'
    unitcosts = [this.getLinkCost(); this.getNodeCost()];
	case 'NormalizeCost'
		linkcosts = mean(this.getLinkCost());
		nodecosts = mean(this.getNodeCost());
		unitcosts = [linkcosts*Ne./(sum(1./capacities(1:Ne)).*capacities(1:Ne));
			nodecosts*Ndc./(sum(1./capacities(Ne+(1:Ndc))).*capacities(Ne+(1:Ndc)))];
  otherwise
    error('error: %s.',cost_init_method);
end
defaultopts = Dictionary('ResidualCapacity', this.convertParameter(capacities), ...
	'Stage', 'temp', ...
	'bFinal', false, ...
	'CapacityConstrained', false,...
	'bCountTime', false);
if nargin <= 2
	options = defaultopts;
else
	options = structmerge(defaultopts, options);
end

if isfield(options, 'Fields')
	fields = options.Fields;
	options = rmfield(options, 'Fields');
else
	fields = 'all';
end
options.bInitialize = true;
end