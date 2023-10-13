	function [vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT] = ...
			getTraceOffsetOne(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur,dblSuperResFactor)
	%getTraceOffsetOne Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT] = ...
	%	getTraceOffsetOne(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur,dblSuperResFactor)
	%
	%This is a subfunction for getZeta().
	
	%% prepare
	vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur,dblSuperResFactor);
	
	%get data
	%build interpolated data
	matTracePerTrial = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT);
	indRemPoints = vecRefT<0 | vecRefT>dblUseMaxDur;
	vecRefT(indRemPoints) = [];
	matTracePerTrial(:,indRemPoints)=[];
	vecMeanTrace = nanmean(matTracePerTrial,1)';
	dblMin = min(vecMeanTrace);
	dblMax = max(vecMeanTrace);
	dblRange = (dblMax-dblMin);
	if dblRange == 0
		dblRange = 1;
		warning([mfilename ':ZeroVar'],'Input data has zero variance');
	end
	vecMeanTrace = (vecMeanTrace-dblMin)./dblRange;
	vecThisFrac = cumsum(vecMeanTrace) / sum(vecMeanTrace);
	
	%get linear fractions
	vecThisFracLinear = linspace(mean(vecMeanTrace),sum(vecMeanTrace),numel(vecMeanTrace))' / sum(vecMeanTrace);
	
	%assign data
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

