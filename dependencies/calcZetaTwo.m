function [vecSpikeT,vecRealDiff,vecRealFrac1,vecRealFrac2,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcZetaTwo(vecSpikeTimes1,vecEventStarts1,vecSpikeTimes2,vecEventStarts2,dblUseMaxDur,intResampNum,boolDirectQuantile,boolUseParallel)
	%calcZetaTwo Calculates neuronal responsiveness difference
	%[vecSpikeT,vecRealDiff,vecRealFrac1,vecRealFrac2,vecRealFracLinear,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcZetaTwo(vecSpikeTimes1,vecEventStarts1,vecSpikeTimes2,vecEventStarts2,dblUseMaxDur,intResampNum,boolDirectQuantile,boolUseParallel)
	
	%% check inputs and pre-allocate error output
	vecSpikeT = [];
	vecRealDiff = [];
	vecRealFrac1 = [];
	vecRealFrac2 = [];
	cellRandT = [];
	cellRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	if ~exist('boolUseParallel','var') || isempty(boolUseParallel)
        objPool = gcp('nocreate');
        if isempty(objPool) || ~isprop(objPool,'NumWorkers') || objPool.NumWorkers < 4
            boolUseParallel = false;
        else
            boolUseParallel = true;
        end
	end
	
	%% reduce spikes
	if size(vecEventStarts1,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	vecEventT1 = vecEventStarts1(:,1);
	dblStartT = max([vecSpikeTimes1(1) min(vecEventT1)-dblUseMaxDur]);
	dblStopT = max(vecEventT1)+dblUseMaxDur*2;
	vecSpikeTimes1(vecSpikeTimes1 < dblStartT | vecSpikeTimes1 > dblStopT) = [];
	vecEventT2 = vecEventStarts2(:,1);
	dblStartT = max([vecSpikeTimes2(1) min(vecEventT2)-dblUseMaxDur]);
	dblStopT = max(vecEventT2)+dblUseMaxDur*2;
	vecSpikeTimes2(vecSpikeTimes2 < dblStartT | vecSpikeTimes2 > dblStopT) = [];

	
	%% get spikes per trial
	[cellTrialPerSpike1,cellTimePerSpike1] = getSpikesInTrial(vecSpikeTimes1,vecEventT1,dblUseMaxDur,true);
	[cellTrialPerSpike2,cellTimePerSpike2] = getSpikesInTrial(vecSpikeTimes2,vecEventT2,dblUseMaxDur,true);
	
	%% run normal
	%normalize to cumsum(v1)+cumsum(v2) = 1
	%take difference
	%mean-subtract
	
	%get difference
	[vecSpikeT,vecRealDiff,vecRealFrac1,vecThisSpikeTimes1,vecRealFrac2,vecThisSpikeTimes2] = ...
		getTempOffsetTwo(cellTimePerSpike1,cellTimePerSpike2,dblUseMaxDur);
	if numel(vecRealDiff) < 2
		return
	end
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	
	%% run bootstraps; try parallel, otherwise run normal loop
	%repeat procedure, but swap trials randomly in each resampling
	cellRandT = cell(1,intResampNum);
	cellRandDiff = cell(1,intResampNum);
	vecMaxRandD = nan(1,intResampNum);
	cellAggregateTrials = cat(1,cellTimePerSpike1,cellTimePerSpike2);
	intT1 = numel(cellTimePerSpike1);
	intT2 = numel(cellTimePerSpike2);
	intTotT = intT1+intT2;
	if boolUseParallel
		parfor intResampling=1:intResampNum
			%% get random subsample
			%if cond1 has 10 trials, and cond2 has 100, then:
			%for shuffle of cond1: take 10 trials from set of 110
			%for shuffle of cond2: take 100 trials from set of 110
			
			vecUseRand1 = randperm(intTotT,intT1);
			cellTimePerSpike1_Rand = cellAggregateTrials(vecUseRand1);
			
			vecUseRand2 = randperm(intTotT,intT2);
			cellTimePerSpike2_Rand = cellAggregateTrials(vecUseRand2);
			if sum(cellfun(@numel,cellTimePerSpike2_Rand)) == 0 && sum(cellfun(@numel,cellTimePerSpike1_Rand)) == 0
				continue;
			end
			
			%get difference
			[vecRandSpikeT,vecRandDiff] = ...
				getTempOffsetTwo(cellTimePerSpike1_Rand,cellTimePerSpike2_Rand,dblUseMaxDur);
			
			%assign data
			cellRandT{intResampling} = vecRandSpikeT;
			cellRandDiff{intResampling} = vecRandDiff;
			vecMaxRandD(intResampling) = max(abs(vecRandDiff));
		end
	else
		for intResampling=1:intResampNum
			%% get random subsample
			%if cond1 has 10 trials, and cond2 has 100, then:
			%for shuffle of cond1: take 10 trials from set of 110
			%for shuffle of cond2: take 100 trials from set of 110
			
			vecUseRand1 = randperm(intTotT,intT1);
			cellTimePerSpike1_Rand = cellAggregateTrials(vecUseRand1);
			
			vecUseRand2 = randperm(intTotT,intT2);
			cellTimePerSpike2_Rand = cellAggregateTrials(vecUseRand2);
			if sum(cellfun(@numel,cellTimePerSpike2_Rand)) == 0 && sum(cellfun(@numel,cellTimePerSpike1_Rand)) == 0
				continue;
			end
			
			%get difference
			[vecRandSpikeT,vecRandDiff] = ...
				getTempOffsetTwo(cellTimePerSpike1_Rand,cellTimePerSpike2_Rand,dblUseMaxDur);
			
			%assign data
			cellRandT{intResampling} = vecRandSpikeT;
			cellRandDiff{intResampling} = vecRandDiff;
			vecMaxRandD(intResampling) = max(abs(vecRandDiff));
		end
	end
	
	%% calculate p
	%take max-dev: zeta_raw
	vecMaxRandD(isnan(vecMaxRandD))=dblMaxD;
	
	%get p-value
	[dblZetaP,dblZETA] = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
end
