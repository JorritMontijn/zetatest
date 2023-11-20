function matTracePerTrial = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT)
	%getInterpolatedTimeSeries Builds common timeframe
	%syntax: matTracePerTrial = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT)
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
	vecEventStartT = sort(vecEventStartT);
    strNew=which('matlab.internal.math.interp1');
    if ~isempty(strNew)
        boolUseNew=true;
    else
        boolUseNew=false;
	end
	boolDebug = false;
	
	if nargout > 1
		error;
	end
	
	%% assign data
	intSampleNum = numel(vecTimestamps);
	dblCritDiff = max(diff(vecRefT))/100;
	matTracePerTrial = nan(numel(vecEventStartT),numel(vecRefT));
	intCurrIdx=1;
	boolSkipToStop = min(diff(vecEventStartT)) > range(vecRefT);
	for intTrial=1:numel(vecEventStartT)
		%% get original times
		dblStartT = vecEventStartT(intTrial);
		%intStartT = max([1 find(vecTimestamps >= (dblStartT + vecRefT(1)),1) - 1]);
		%intStopT = min([numel(vecTimestamps) find(vecTimestamps > (dblStartT + vecRefT(end)),1)]);
		dblLocalStartT=dblStartT+vecRefT(1);
		for intCurrIdx=intCurrIdx:intSampleNum
			if vecTimestamps(intCurrIdx)>dblLocalStartT, break, end
		end
		intStartT = max([1 intCurrIdx-1]);
		intStartIdx = intCurrIdx;
		dblLocalStopT=dblStartT+vecRefT(end);
		for intCurrIdx=intCurrIdx:intSampleNum
			if vecTimestamps(intCurrIdx)>dblLocalStopT, break, end
		end
		intStopT = min([intSampleNum intCurrIdx]);
		if ~boolSkipToStop
			intCurrIdx = intStartIdx+1;
		end
		vecSelectSamples = intStartT:intStopT;
		
		%% get data
		vecUseTimes = vecTimestamps(vecSelectSamples);
		vecUseTrace = vecData(vecSelectSamples);
		
		%% interpolate
		vecUseInterpT = vecRefT+dblStartT;
		
		if boolDebug
			%get real fractions for training set
			vecInterpTrace1 = interp1(vecUseTimes,vecUseTrace,vecUseInterpT);
			if ~all((vecInterpTrace1==vecInterpTrace2) | (isnan(vecInterpTrace1) & isnan(vecInterpTrace2)))
				error
			end
		end
		
		%check if we can perform direct assignment
		boolDirect = false;
		if ndims(vecUseInterpT) == ndims(vecUseTimes)
			for intIdx=1:numel(vecUseTimes)
				if vecUseTimes(intIdx) >= (vecUseInterpT(1)-dblCritDiff)
					intAlignBy = intIdx;
					break
				end
			end
			vecAlignedSamples = (1:numel(vecUseInterpT))+intAlignBy-1;
			if ~any(abs(vecUseTimes(vecAlignedSamples([1 end]))-vecUseInterpT([1 end])) > dblCritDiff)
				boolDirect = true;
			end
		end
		%bypass time-consuming checks and use direct method if possible
		if boolDirect
			vecInterpTrace2 = vecUseTrace(vecAlignedSamples);
		else
			if boolUseNew
				vecInterpTrace2 = matlab.internal.math.interp1(vecUseTimes,vecUseTrace,'linear','none',vecUseInterpT);
			else
				F = griddedInterpolant(vecUseTimes,vecUseTrace,'linear','none');
				vecInterpTrace2 = F(vecUseInterpT);
			end
		end
		
		%assign output
		matTracePerTrial(intTrial,:) = vecInterpTrace2;
	end
end