function [vecRefT,vecRealDiff,vecRealFrac1,vecRealFrac2,matRandDiff,dblZetaP,dblZETA,intZETALoc,matTracePerTrial1,matTracePerTrial2] = ...
		calcTsZetaTwo(vecTraceT1,vecTraceAct1,vecEventStarts1,vecTraceT2,vecTraceAct2,vecEventStarts2,dblSuperResFactor,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize,boolUseParallel)
	%calcTsZetaTwo Calculates neuronal responsiveness index zeta
	%[vecRefT,vecRealDiff,vecRealFrac1,vecRealFrac2,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcTsZetaTwo(vecTraceT1,vecTraceAct1,vecEventStarts1,vecTraceT2,vecTraceAct2,vecEventStarts2,dblSuperResFactor,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize,boolUseParallel)
	
    %% check inputs and pre-allocate error output
	vecRefT = [];
	vecRealDiff = [];
	vecRealFrac1 = [];
	vecRealFrac2 = [];
	matRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	if ~exist('boolUseParallel','var') || isempty(boolUseParallel)
        objPool = gcp('nocreate');
        if isempty(objPool) || ~isprop(objPool,'NumWorkers') || objPool.NumWorkers < 4
            boolUseParallel = false;
        else
            boolUseParallel = false; %seems to always be slower than non-parallel
        end
	end
	
	%% reduce data
	if size(vecEventStarts1,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	%discard leading/lagging data1
	vecEventStarts1 = vecEventStarts1(:,1);
	dblPreUse = -dblUseMaxDur;
	dblPostUse = dblUseMaxDur*2;
	dblStartT1 = min(vecEventStarts1) + dblPreUse*2;
	dblStopT1 = max(vecEventStarts1) + dblPostUse*2;
	indRemoveEntries1 = (vecTraceT1 < dblStartT1) | (vecTraceT1 > dblStopT1);
	vecTraceT1(indRemoveEntries1) = [];
	vecTraceAct1(indRemoveEntries1) = [];
	
	%discard leading/lagging data2
	vecEventStarts2 = vecEventStarts2(:,1);
	dblStartT2 = min(vecEventStarts2) + dblPreUse*2;
	dblStopT2 = max(vecEventStarts2) + dblPostUse*2;
	indRemoveEntries2 = (vecTraceT2 < dblStartT2) | (vecTraceT2 > dblStopT2);
	vecTraceT2(indRemoveEntries2) = [];
	vecTraceAct2(indRemoveEntries2) = [];
	
	%rescale
	dblMin = min(min(vecTraceAct1),min(vecTraceAct2));
	dblMax = max(max(vecTraceAct1),max(vecTraceAct2));
	dblRange = (dblMax-dblMin);
	if dblRange == 0
		dblRange = 1;
		warning([mfilename ':ZeroVar'],'Input data has zero variance');
	end
	vecTraceAct1 = (vecTraceAct1-dblMin)./dblRange;
	vecTraceAct2 = (vecTraceAct2-dblMin)./dblRange;
	
	%% build reference time and matrices
	%time
	vecRefT1 = getTsRefT(vecTraceT1,vecEventStarts1,dblUseMaxDur);
	vecRefT2 = getTsRefT(vecTraceT2,vecEventStarts2,dblUseMaxDur);
	%set tol
	dblSampInterval = (median(diff(vecRefT1)) + median(diff(vecRefT2)))/2;
	dblTol = dblSampInterval/dblSuperResFactor;
	vecRefT = uniquetol(cat(1,vecRefT1(:),vecRefT2(:)),dblTol);
	intT = numel(vecRefT);
	
	%matrices
	matTracePerTrial1 = getInterpolatedTimeSeries(vecTraceT1,vecTraceAct1,vecEventStarts1,vecRefT);
	matTracePerTrial2 = getInterpolatedTimeSeries(vecTraceT2,vecTraceAct2,vecEventStarts2,vecRefT);
	
	%% get trial responses
	[vecRealDiff,vecRealFrac1,vecRealFrac2] = getTraceOffsetTwo(matTracePerTrial1,matTracePerTrial2);
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	
	%repeat procedure, but swap trials randomly in each resampling
	matRandDiff = nan(intResampNum,intT);
	vecMaxRandD = nan(1,intResampNum);
	matAggregateTrials = cat(1,matTracePerTrial1,matTracePerTrial2);
	intTrials1 = size(matTracePerTrial1,1);
	intTrials2 = size(matTracePerTrial2,1);
	intTotTrials = intTrials1+intTrials2;
	if boolUseParallel
		parfor intResampling=1:intResampNum
			%% get random subsample
			%if cond1 has 10 trials, and cond2 has 100, then:
			%for shuffle of cond1: take 10 trials from set of 110
			%for shuffle of cond2: take 100 trials from set of 110
			vecUseRand1 = randi(intTotTrials,[1,intTrials1]);
			vecUseRand2 = randi(intTotTrials,[1,intTrials2]);
			
			matTrace1_Rand = matAggregateTrials(vecUseRand1,:);
			matTrace2_Rand = matAggregateTrials(vecUseRand2,:);
			
			%get difference
			vecRandDiff = getTraceOffsetTwo(matTrace1_Rand,matTrace2_Rand);
			
			%assign data
			matRandDiff(intResampling,:) = vecRandDiff;
			vecMaxRandD(intResampling) = max(abs(vecRandDiff));
		end
	else
		for intResampling=1:intResampNum
			%% get random subsample
			%if cond1 has 10 trials, and cond2 has 100, then:
			%for shuffle of cond1: take 10 trials from set of 110
			%for shuffle of cond2: take 100 trials from set of 110
			vecUseRand1 = randi(intTotTrials,[1,intTrials1]);
			vecUseRand2 = randi(intTotTrials,[1,intTrials2]);

			matTrace1_Rand = matAggregateTrials(vecUseRand1,:);
			matTrace2_Rand = matAggregateTrials(vecUseRand2,:);
			
			%get difference
			vecRandDiff = getTraceOffsetTwo(matTrace1_Rand,matTrace2_Rand);
			
			%assign data
			matRandDiff(intResampling,:) = vecRandDiff;
			vecMaxRandD(intResampling) = max(abs(vecRandDiff));
		end
	end
	
	%% calculate p
	%take max-dev: zeta_raw
	vecMaxRandD(isnan(vecMaxRandD))=dblMaxD;
	
	%get p-value
	[dblZetaP,dblZETA] = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
end

