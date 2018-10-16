function [sp_profit, b_violate, violates] = SolveSCPCC(this, slices, prices, capacities, options)
    Ns = length(slices);
    Ne = this.NumberLinks;
    Ndc = this.NumberDataCenters;
    num_process = Ns;
    output_k = cell(num_process,1);
    loads = zeros(Ne+Ndc, num_process);
    for si = 1:Ns
        sl = slices{si};
        dc_id = sl.getDCPI;
        link_id = sl.VirtualLinks.PhysicalLink;
        sl.prices.Link = prices.Link(link_id);
        sl.prices.Node = prices.Node(dc_id);
    end
    nout = nargout;
		if options.CapacityConstrained
			for si = 1:Ns
				slices{si}.setProblem([],[],[],...
					capacities.Link(1:Ne,si), capacities.Node(Ne+(1:Ndc),si));
			end
		end
		M = getParallelInfo();
    parfor (sj = 1:Ns,M)
    %for sj = 1:num_process
        sl = slices{sj};
        if nout >= 2
            [output_k{sj}, loads(:,sj)] = sl.priceOptimalFlowRateCC([], options);
        else
            output_k{sj} = sl.priceOptimalFlowRateCC([], options);
        end
    end
    for si = 1:Ns
        sl = slices{si};
        sl.saveResults(output_k{si});
        sl.prices.Link = [];
        sl.prices.Node = [];
        sl.capacities = [];
    end
    
    %% Output
    sp_profit = this.getSliceProviderProfit(prices, ...
			struct({'Slices','PricingPolicy'}, {'linear', slices}));
    if nargout >= 2
			if ~options.CapacityConstrained
        violates = sum(loads,2) > options.capacities;
        b_violate = ~isempty(find(violates,1)); 
			else
				violates = [];
				b_violate = false;
			end
    end
end
