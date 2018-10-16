classdef SlicingMethod < uint32
	
	enumeration
		% The value is 32-bit unsigned integer.
		StaticPartition(hex2dec('0000'));
		StaticPricing(hex2dec('0310'));
		SingleNormal(hex2dec('0313'));
		SingleFunction(hex2dec('10313'));
		AdjustPricing(hex2dec('0111'));
		DualPricing(hex2dec('0211'));
		%FactorPricing(hex2dec('0211'));		% factor cost to pricing
%%
		DynamicPartition(hex2dec('0001'));
		PartitionPricing(hex2dec('0021'));
	end
	
	
	methods
		%% Static v.s Dynamic (1-bit)
		% If resource is re-allocate.
		function tf = IsStatic(this)
			tf = bitand(this, hex2dec('0001')) == 0;
		end
		
		%% Single (1-bit)
		% all flows from each slices are considered together.
		function tf = IsSingle(this)
			tf = bitand(this, hex2dec('0002')) ~= 0;
		end		
		
		%% Resource Allocation (4-bit)
		% Partition, Pricing, ... 
		function tf = IsPartition(this)
			tf = bitshift(bitand(this, hex2dec('00F0')),-4) == 0; 
		end
		function tf = IsPricing(this)
			tf = bitshift(bitand(this, hex2dec('00F0')),-4) == 1; 
		end
		
		%% Pricing (4-bit)
		% Using pricing mechanism to allocate resources.
		function tf = IsFixedPricing(this)
			b = bitshift(bitand(this, hex2dec('0F00')), -8);
			tf = this.IsPricing && (b == 0 || b==3); 
		end
		function tf = IsAdjustPricing(this)
			tf = this.IsPricing && (bitshift(bitand(this, hex2dec('0F00')), -8) == 1); 
		end
		function tf = IsDualPricing(this)
			tf = this.IsPricing && (bitshift(bitand(this, hex2dec('0F00')), -8) == 2); 
		end
		function tf = IsFactorPricing(this)
			tf = this.IsPricing && (bitshift(bitand(this, hex2dec('0F00')), -8) == 3); 
		end

		%% Single Function (TEST)
		function tf = IsSingleFunction(this)
			tf = this.IsSingle && (bitshift(bitand(this, hex2dec('10000')), -16) ~= 0); 
		end
		
	end
end

