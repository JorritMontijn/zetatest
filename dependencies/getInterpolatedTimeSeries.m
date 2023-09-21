function matTracePerTrial = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT,cellSampleAssignments)
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
	boolDebug = false;
	
	if nargout > 1
		error;
	end
	
	%% assign data
	matTracePerTrial = nan(numel(vecEventStartT),numel(vecRefT));
	for intTrial=1:numel(vecEventStartT)
		%% get original times
		dblStartT = vecEventStartT(intTrial);
		vecSelectSamples = cellSampleAssignments{intTrial};
		
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
		
		%bypass time-consuming checks and use direct method
        if boolUseNew
            vecInterpTrace2 = matlab.internal.math.interp1(vecUseTimes,vecUseTrace,'linear','none',vecUseInterpT);
        else
            F = griddedInterpolant(vecUseTimes,vecUseTrace,'linear','none');
            vecInterpTrace2 = F(vecUseInterpT);
		end
		
		%assign output
		matTracePerTrial(intTrial,:) = vecInterpTrace2;
	end
end