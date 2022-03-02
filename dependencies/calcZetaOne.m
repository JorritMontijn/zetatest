function [vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = calcZetaOne(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize,boolStitch)
	%calcZetaOne Calculates neuronal responsiveness index zeta
	%[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcZetaOne(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize,boolStitch)
	
	%% check inputs and pre-allocate error output
	vecSpikeT = [];
	vecRealDiff = [];
	vecRealFrac = [];
	vecRealFracLinear = [];
	cellRandT = [];
	cellRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end
	if ~exist('boolStitch','var') || isempty(boolStitch)
		boolStitch = true;
	end
	
	%% reduce spikes
	if size(vecEventStarts,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	vecEventT = vecEventStarts(:,1);
	dblStartT = max([vecSpikeTimes(1) min(vecEventT)-dblUseMaxDur*5*dblJitterSize]);
	dblStopT = max(vecEventT)+dblUseMaxDur*5*dblJitterSize;
	vecSpikeTimes(vecSpikeTimes < dblStartT | vecSpikeTimes > dblStopT) = [];
	if numel(vecSpikeTimes) < 3
		return;
	end
	
	%% build pseudo data, stitching stimulus periods
	if boolStitch
		[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(vecSpikeTimes,vecEventT,dblUseMaxDur);
	else
		vecPseudoSpikeTimes = vecSpikeTimes;
		vecPseudoEventT = vecEventT;
	end
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFrac,vecRealFracLinear,vecSpikeT] = ...
		getTempOffsetOne(vecPseudoSpikeTimes,vecPseudoEventT,dblUseMaxDur);
	if numel(vecRealDiff) < 3
		return
	end
	vecRealDiff = vecRealDiff - mean(vecRealDiff);
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	
	%% run bootstraps; try parallel, otherwise run normal loop
	% run pre-set number of iterations
	cellRandT = cell(1,intResampNum);
	cellRandDiff = cell(1,intResampNum);
	vecMaxRandD = nan(1,intResampNum);
	vecStartOnly = vecPseudoEventT(:);
	intTrials = numel(vecStartOnly);
	%vecJitterPerTrial = dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2); %original zeta
	vecJitterPerTrial = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials)'; %new
	matJitterPerTrial = nan(intTrials,intResampNum);
	for intResampling=1:intResampNum
		matJitterPerTrial(:,intResampling) = vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
	end
	try
		parfor intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			[vecRandDiff,vecThisSpikeFracs,vecThisFracLinear,vecThisSpikeTimes] = ...
				getTempOffsetOne(vecPseudoSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
			
			%assign data
			cellRandT{intResampling} = vecThisSpikeTimes;
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	catch
		for intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			[vecRandDiff,vecThisSpikeFracs,vecThisFracLinear,vecThisSpikeTimes] = ...
				getTempOffsetOne(vecPseudoSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
			
			%assign data
			cellRandT{intResampling} = vecThisSpikeTimes;
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	end
	
	%% calculate significance
	[dblZetaP,dblZETA] = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile);
end

