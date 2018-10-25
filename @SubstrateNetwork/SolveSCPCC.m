function [sp_profit, b_violate, violates] = SolveSCPCC(this, slices, prices, options)
Ns = length(slices);
Ne = this.NumberLinks;
Ndc = this.NumberDataCenters;
num_process = Ns;
output_k = cell(num_process,1);
loads = zeros(Ne+Ndc, num_process);
for si = 1:Ns
  sl = slices{si};
  dc_id = sl.getDCPI;
  link_id = sl.Links.PhysicalLink;
  sl.prices.Link = prices.Link(link_id);
  sl.prices.Node = prices.Node(dc_id);  % |prices.Node| only contain the price of data center nodes.
end
nout = nargout;
M = getParallelInfo();
if nargin >= 4
  if options.CapacityConstrained
    for si = 1:Ns
      caps = struct('Link', options.Capacities(1:Ne,si), ...
        'Node', options.Capacities(Ne+(1:Ndc),si));
      slices{si}.op.setProblem(caps);
    end
  end
  for si = 1:Ns
    slices{si}.op.update_options(options);
  end
end
parfor (sj = 1:Ns,M)
  %for sj = 1:num_process
  sl = slices{sj};
  if nout >= 2
    [output_k{sj}, loads(:,sj)] = sl.op.priceOptimalFlowRate([], options);
  else
    output_k{sj} = sl.op.priceOptimalFlowRate([], options);
  end
end
for si = 1:Ns
  sl = slices{si};
  sl.op.saveTempResults(output_k{si});
end

%% Output
sp_profit = this.getSliceProviderProfit(prices, ...
  getstructfields(options, {'Slices','PricingPolicy'}, 'error'));
if nargout >= 2
  if ~options.CapacityConstrained
    violates = sum(loads,2) > options.Capacities;
    b_violate = ~isempty(find(violates,1));
  else
    violates = [];
    b_violate = false;
  end
end
end
