function [vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT] = ...
		getTraceOffsetSR(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur)
	%getTraceOffsetSR Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT,vecMeanTrace] = ...
	%	getTraceOffsetSR(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur)
	%
	%This is a subfunction for zetatstest().
	
	%% prepare
	vecEventStartT = sort(vecEventStartT);
	%pre-allocate
	intTrialNum = numel(vecEventStartT);
	intTimeNum = numel(vecTimestamps);
	%build common timeframe
	cellRefT = cell(1,intTrialNum);
	cellVals = cell(1,intTrialNum);
	for intTrial=1:intTrialNum
		% get original times
		dblStartT = vecEventStartT(intTrial);
		dblStopT = dblStartT+dblUseMaxDur;
		intStartT = max([1 find(vecTimestamps > dblStartT,1) - 1]);
		intStopT = min([intTimeNum find(vecTimestamps > dblStopT,1)]);
		vecSelectSamples = intStartT:intStopT;
		
		%% get data
		cellRefT{intTrial} = vecTimestamps(vecSelectSamples)-dblStartT;
		cellVals{intTrial} = vecData(vecSelectSamples);
	end
	
	%get data
	%build interpolated data
	[vecRefT,vecReorder] = sort(cell2vec(cellRefT));
	vecVals = cell2vec(cellVals);
	vecMeanTrace = vecVals(vecReorder);
	indRemPoints = vecRefT<0 | vecRefT>dblUseMaxDur;
	vecMeanTrace(indRemPoints) = [];
	vecRefT(indRemPoints) = [];
	vecThisFrac = cumsum(vecMeanTrace) / sum(vecMeanTrace);
	
	%get linear fractions
	vecThisFracLinear = linspace(mean(vecMeanTrace),sum(vecMeanTrace),numel(vecMeanTrace))' / sum(vecMeanTrace);
	
	%assign data
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

