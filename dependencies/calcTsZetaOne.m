function [vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcTsZetaOne(vecTraceT,vecTraceAct,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	%calcTsZeta Calculates neuronal responsiveness index zeta for timeseries data
	%[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcTsZetaOne(vecTraceT,vecTraceAct,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	
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
		getTraceOffsetOne(vecPseudoT,vecPseudoTrace,vecPseudoStartT',vecRefT,dblUseMaxDur);
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	intSamples = numel(vecRealDiff);
	intTrials = numel(vecPseudoStartT);
	
	%% run bootstraps; try parallel, otherwise run normal loop
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
			[vecRandDiff,vecThisFrac,vecThisFracLinear,vecRandT] = getTraceOffsetOne(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,vecRefT,dblUseMaxDur);
			
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
			[vecRandDiff,vecThisFrac,vecThisFracLinear,vecRandT] = getTraceOffsetOne(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,vecRefT,dblUseMaxDur);
			
			%assign data
			cellRandT{intResampling} = vecRandT;
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	end
	
	%% calculate significance
	[dblZetaP,dblZETA] = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
end

