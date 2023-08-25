	function [vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT] = ...
			getTraceOffsetOne(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur)
	%getTraceOffsetOne Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT] = ...
	%	getTraceOffsetOne(vecTimestamps,vecData,vecEventStartT,vecRefT,dblUseMaxDur)
	%
	%This is a subfunction for getZeta().
	
	%% prepare
	vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur);
	
	%get data
	%build interpolated data
	[vecRefT,matTracePerTrial] = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur,vecRefT);
	indRemPoints = vecRefT<0 | vecRefT>dblUseMaxDur;
	vecRefT(indRemPoints) = [];
	matTracePerTrial(:,indRemPoints)=[];
	vecMeanTrace = nanmean(matTracePerTrial,1)';
	vecThisFrac = cumsum(vecMeanTrace) / sum(vecMeanTrace);
	
	%get linear fractions
	vecThisFracLinear = linspace(mean(vecMeanTrace),sum(vecMeanTrace),numel(vecMeanTrace))' / sum(vecMeanTrace);
	
	%assign data
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

