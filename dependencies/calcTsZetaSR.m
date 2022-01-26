function [vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcTsZetaSR(vecTraceT,vecTraceAct,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	%calcTraceZeta Calculates neuronal responsiveness index zeta
	%[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,vecMeanTrace,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcTraceZeta(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingInterval,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	
	%% check inputs and pre-allocate error output
	vecRefT = [];
	vecRealDiff = [];
	vecRealFrac = [];
	vecRealFracLinear = [];
	cellRandT = [];
	cellRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	
	%% reduce data
	if size(vecEventStarts,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	%discard leading/lagging data
	vecEventStarts = vecEventStarts(:,1);
	dblPreUse = -dblUseMaxDur*dblJitterSize;
	dblPostUse = dblUseMaxDur*(dblJitterSize+1);
	
	dblStartT = min(vecEventStarts) + dblPreUse*2;
	dblStopT = max(vecEventStarts) + dblPostUse*2;
	indRemoveEntries = (vecTraceT < dblStartT) | (vecTraceT > dblStopT);
	vecTraceT(indRemoveEntries) = [];
	vecTraceAct(indRemoveEntries) = [];
	
	%stitch trials
	[vecPseudoT,vecPseudoTrace,vecPseudoStartT] = getPseudoTimeSeries(vecTraceT,vecTraceAct,vecEventStarts,dblUseMaxDur);
	vecPseudoTrace = vecPseudoTrace - min(vecPseudoTrace(:));
	if numel(vecPseudoT) < 3
		return;
	end
	
	
	%% get trial responses
	[vecRealDiff,vecRealFrac,vecRealFracLinear,vecRefT] = ...
		getTraceOffsetSR(vecPseudoT,vecPseudoTrace,vecPseudoStartT',dblUseMaxDur);
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	intSamples = numel(vecRealDiff);
	intTrials = numel(vecPseudoStartT);
	
	%% run bootstraps; try parallel, otherwise run normal loop
	if intResampNum < 1
		%% get parameters
		dblConvAt = intResampNum; %convergence fraction
		boolUseParallel = true;
		try
			objP = gcp('nocreate');
			if isempty(objP)
				%create pool
				objP=gcp;
				%test parfor
				parfor testi=1:10
					x=rand(10);
				end
			end
			intBatchSize = objP.NumWorkers; %size of pool
		catch
			boolUseParallel = false;
			intBatchSize = 20; %size of pool
		end
		
		%% run initial set
		[vecMaxRandD,cellRandDiff,cellRandT] = getRandTsZetaBatchSR(vecPseudoT,vecPseudoTrace,vecPseudoStartT,dblUseMaxDur,dblJitterSize,intBatchSize,boolUseParallel);
		
		%calculate initial p
		dblOldZetaP = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
		intCurrIter = intBatchSize;
		
		%% continue until converged
		intMaxIters = 1000;
		boolConverged = false;
		while ~boolConverged && intCurrIter < intMaxIters
			[vecBatchMaxRandD,cellBatchRandDiff,cellBatchRandT] = getRandTsZetaBatchSR(vecPseudoT,vecPseudoTrace,vecPseudoStartT,dblUseMaxDur,dblJitterSize,intBatchSize,boolUseParallel);
			intCurrIter = intCurrIter + intBatchSize;
			vecMaxRandD = cat(2,vecMaxRandD,vecBatchMaxRandD);
			cellRandDiff = cat(2,cellRandDiff,cellBatchRandDiff);
			cellRandT = cat(2,cellRandT,cellBatchRandT);
			%check if p is stable
			[dblZetaP,dblZETA] = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
			dblChange = abs(1 - (dblZetaP / dblOldZetaP));
			dblOldZetaP = dblZetaP;
			if dblChange < dblConvAt || (dblZetaP == 0 && dblOldZetaP == 0)
				boolConverged = true;
			end
		end
	else
		% run pre-set number of iterations
		cellRandT = cell(1,intResampNum);
		cellRandDiff = cell(1,intResampNum);
		vecMaxRandD = nan(1,intResampNum);
		vecStartOnly = vecPseudoStartT(:);
		%vecJitterPerTrial = dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2); %same as original zeta
		vecJitterPerTrial = dblJitterSize*linspace(dblUseMaxDur/intTrials,dblUseMaxDur,intTrials)';
		matJitterPerTrial = nan(intTrials,intResampNum);
		for intResampling=1:intResampNum
			matJitterPerTrial(:,intResampling) = vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
		end
		try
			parfor intResampling=1:intResampNum
				%% get random subsample
				vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
				
				%get temp offset
				[vecRandDiff,vecThisFrac,vecThisFracLinear,vecRandT] = getTraceOffsetSR(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,dblUseMaxDur);
			
				%assign data
				cellRandT{intResampling} = vecRandT;
				cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
				vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
			end
		catch
			for intResampling=1:intResampNum
				%% get random subsample
				vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
				
				%get temp offset
				[vecRandDiff,vecThisFrac,vecThisFracLinear,vecRandT] = getTraceOffsetSR(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,dblUseMaxDur);
				
				%assign data
				cellRandT{intResampling} = vecRandT;
				cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
				vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
			end
		end
		
		%% calculate significance
		[dblZetaP,dblZETA] = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
	end
end

