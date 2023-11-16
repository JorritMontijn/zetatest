function [vecRefT,vecRealDiff,vecRealFrac1,vecRealFrac2,matRandDiff,dblZetaP,dblZETA,intZetaIdx,matCond1_ext,matCond2_ext] = ...
			calcRmZetaTwo(vecT1,matCond1,vecT2,matCond2,intResampNum,boolDirectQuantile)
	%calcRmZetaTwo Calculates neuronal responsiveness index zeta
	%[vecRefT,vecRealDiff,vecRealFrac1,vecRealFrac2,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcRmZetaTwo(vecTraceT1,vecTraceAct1,vecEventStarts1,vecTraceT2,vecTraceAct2,vecEventStarts2,dblSuperResFactor,dblUseMaxDur,intResampNum,boolDirectQuantile,boolUseParallel)
	
    %% check inputs and pre-allocate error output
	vecRefT = [];
	vecRealDiff = [];
	vecRealFrac1 = [];
	vecRealFrac2 = [];
	matRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZetaIdx = nan;
	if ~exist('boolUseParallel','var') || isempty(boolUseParallel)
        objPool = gcp('nocreate');
        if isempty(objPool) || ~isprop(objPool,'NumWorkers') || objPool.NumWorkers < 4
            boolUseParallel = false;
        else
            boolUseParallel = false; %seems to always be slower than non-parallel
        end
	end
	
	%% reduce data
	%rescale
	dblMin = min(min(matCond1(:)),min(matCond2(:)));
	dblMax = max(max(matCond1(:)),max(matCond2(:)));
	dblRange = (dblMax-dblMin);
	if dblRange == 0
		dblRange = 1;
		warning([mfilename ':ZeroVar'],'Input data has zero variance');
	end
	matCond1 = (matCond1-dblMin)./dblRange;
	matCond2 = (matCond2-dblMin)./dblRange;
	
	%% build reference time and matrices
	vecRefT = sort(unique(cat(1,vecT1(:),vecT2(:))));
	indRemT = vecRefT<min(vecT1(:)) | vecRefT<min(vecT2(:)) | vecRefT>max(vecT1(:)) | vecRefT>max(vecT2(:));
	vecRefT(indRemT) = [];
	intT = numel(vecRefT);
	
	%% assign data
	matCond1_ext = getInterpolatedMeasures(vecRefT,vecT1,matCond1);
	matCond2_ext = getInterpolatedMeasures(vecRefT,vecT2,matCond2);

	%% get trial responses
	[vecRealDiff,vecRealFrac1,vecRealFrac2] = getTraceOffsetTwo(matCond1_ext,matCond2_ext);
	[dblMaxD,intZetaIdx]= max(abs(vecRealDiff));
	
	%repeat procedure, but swap trials randomly in each resampling
	matRandDiff = nan(intResampNum,intT);
	vecMaxRandD = nan(1,intResampNum);
	matAggregateTrials = cat(1,matCond1_ext,matCond2_ext);
	intTrials1 = size(matCond1_ext,1);
	intTrials2 = size(matCond2_ext,1);
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

