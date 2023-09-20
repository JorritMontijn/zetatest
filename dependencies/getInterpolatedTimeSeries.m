function matTracePerTrial = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT)
	%getInterpolatedTimeSeries Builds common timeframe
	%syntax: [vecRefT,matTracePerTrial] = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT)
	%	input:
	%	- vecTimestamps; times (s)
	%	- vecData; data
	%	- vecEventStartT: trial start times (s)
	%	- vecRefT: reference time
	%
	%Version history:
	%1.0 - June 26 2019
	%	Created by Jorrit Montijn
	%2.0 - September 20 2023
	%	Various speed-ups and syntax changes
	
	%% check which interp1 bypass to use
    strNew=which('matlab.internal.math.interp1');
    if ~isempty(strNew)
        boolUseNew=true;
    else
        boolUseNew=false;
	end
	
	%% assign data
	matTracePerTrial = nan(numel(vecEventStartT),numel(vecRefT));
	for intTrial=1:numel(vecEventStartT)
		%% get original times
		dblStartT = vecEventStartT(intTrial);
		intStartT = max([1 find(vecTimestamps > (dblStartT + vecRefT(1)),1) - 1]);
		intStopT = min([numel(vecTimestamps) find(vecTimestamps > (dblStartT + vecRefT(end)),1) + 1]);
		vecSelectSamples = intStartT:intStopT;
		
		%% get data
		vecUseTimes = vecTimestamps(vecSelectSamples);
		vecUseTrace = vecData(vecSelectSamples);
		
		%% interpolate
		vecUseInterpT = vecRefT+dblStartT;
		
		%get real fractions for training set
        %vecInterpTrace = interp1(vecUseTimes,vecUseTrace,vecUseInterpT);
		
        %bypass time-consuming checks and use direct method
        if boolUseNew
            vecInterpTrace = matlab.internal.math.interp1(vecUseTimes,vecUseTrace,'linear','linear',vecUseInterpT);
        else
            F = griddedInterpolant(vecUseTimes,vecUseTrace,'linear');
            vecInterpTrace = F(vecUseInterpT);
        end
		matTracePerTrial(intTrial,:) = vecInterpTrace;
	end
end