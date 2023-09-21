function vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur)
	%getTsRefT Get reference time vector
	%   vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur)
	
	%pre-allocate
	vecEventStartT = sort(vecEventStartT);
	intTrialNum = numel(vecEventStartT);
	intTimeNum = numel(vecTimestamps);
	%build common timeframe
	cellRefT = cell(1,intTrialNum);
	for intTrial=1:intTrialNum
		% get original times
		dblStartT = vecEventStartT(intTrial);
		dblStopT = dblStartT+dblUseMaxDur;
		intStartT = max(1,find(vecTimestamps > dblStartT,1) - 1);
		intStopT = min(intTimeNum,find(vecTimestamps > dblStopT,1));
		vecSelectSamples = intStartT:intStopT;
		
		%% get data
		cellRefT{intTrial} = vecTimestamps(vecSelectSamples)-dblStartT;
	end
	
	%set tol
	dblSampInterval = median(diff(vecTimestamps));
	dblTol = dblSampInterval/100;
	vecRefT = uniquetol(sort(cell2vec(cellRefT)),dblTol);
end

