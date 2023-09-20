function matTracePerTrial = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT,cellSampleAssignments)
	%getTraceInTrial Builds common timeframe
	%syntax: matTracePerTrial = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,vecRefT,cellSampleAssignments)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%
	%Version history:
	%1.0 - June 26 2019
	%	Created by Jorrit Montijn
	
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
		vecSelectSamples = cellSampleAssignments{intTrial};
		
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
            F = griddedInterpolant(X,V,'linear');
            vecInterpTrace = F(Xqcol);
        end
		matTracePerTrial(intTrial,:) = vecInterpTrace;
	end
end