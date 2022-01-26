	function [vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT] = ...
			getTraceOffset(vecTraceT,vecTraceAct,vecStimUseOnTime,dblSamplingInterval,dblUseMaxDur)
	%getTraceOffset Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT,vecMeanTrace] = ...
	%	getTraceOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblSamplingInterval,dblUseMaxDur)
	%
	%This is a subfunction for getZeta().
	
	
	%get data
	[vecRefT,matTracePerTrial] = getTraceInTrial(vecTraceT,vecTraceAct,vecStimUseOnTime,dblSamplingInterval,dblUseMaxDur);
	vecMeanTrace = nanmean(matTracePerTrial,1)';
	vecThisFrac = cumsum(vecMeanTrace) / sum(vecMeanTrace);
	
	%get linear fractions
	vecThisFracLinear = linspace(mean(vecMeanTrace),sum(vecMeanTrace),numel(vecMeanTrace))' / sum(vecMeanTrace);
	
	%assign data
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

