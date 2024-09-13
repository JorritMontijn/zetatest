function vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur,dblSuperResFactor)
	%getTsRefT Get reference time vector
	%   vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur,dblSuperResFactor)
	
	%default
	if ~exist('dblSuperResFactor','var') || isempty(dblSuperResFactor) || dblSuperResFactor < 1 || dblSuperResFactor > 1e6
		dblSuperResFactor = 1;
	end
	
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
	if dblSuperResFactor == 1
		[dummy,intUseEntry]=max(cellfun(@numel,cellRefT));
		vecRefT = cellRefT{intUseEntry};
		dblMedDiff = median(diff(vecRefT));
		vecRefT = round(10*(vecRefT/dblMedDiff))/(10/dblMedDiff);
	else
		dblSampInterval = median(diff(vecTimestamps));
		dblTol = dblSampInterval/dblSuperResFactor;
		% vecRefT = uniquetol(sort(cell2vec(cellRefT)),dblTol);
        vecVals = sort(cell2vec(cellRefT));
		vecRefT = myUniquetol(vecVals,dblTol);
	end
end

function vecData = myUniquetol(vecData,dblTol)
% uniquetol function to match equivalent used in python
vecData = unique(floor(vecData/dblTol)*dblTol);
end
