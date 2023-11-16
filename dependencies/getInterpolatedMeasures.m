function matCond_ext = getInterpolatedMeasures(vecRefT,vecT,matCond)
	%getInterpolatedTimeSeries Builds common timeframe
	%syntax: matCond_ext = getInterpolatedMeasures(vecRefT,vecT,matCond)
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
	intTotN = size(matCond,1);
	matCond_ext = nan(intTotN,numel(vecRefT));
	for intN=1:intTotN
		%% get data
		vecUseTimes = vecT;
		vecUseTrace = matCond(intN,:);
		
		%% interpolate
		%bypass time-consuming checks and use direct method
		if boolUseNew
			vecInterpTrace = matlab.internal.math.interp1(vecUseTimes,vecUseTrace,'linear','none',vecUseInterpT);
		else
			vecInterpTrace = interp1(vecUseTimes,vecUseTrace,vecRefT,'linear',nan);
		end
		
		%assign output
		matCond_ext(intN,:) = vecInterpTrace;
	end
end