function [vecSpikeT,vecRealDiff,vecRealFrac1,vecRealFrac2,vecRealFracLinear,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcZetaDiff(vecSpikeTimes1,vecEventStarts1,vecSpikeTimes2,vecEventStarts2,boolPairedTest,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	%calcZetaDiff Calculates neuronal responsiveness difference
	%[vecSpikeT,vecRealDiff,vecRealFrac1,vecRealFrac2,vecRealFracLinear,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcZetaDiff(vecSpikeTimes1,vecEventStarts1,vecSpikeTimes2,vecEventStarts2,boolPairedTest,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	
	%% check inputs and pre-allocate error output
	vecSpikeT = [];
	vecRealDiff = [];
	vecRealFrac1 = [];
	vecRealFrac2 = [];
	vecRealFracLinear = [];
	cellRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	
	%% reduce spikes
	if size(vecEventStarts1,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	vecEventT1 = vecEventStarts1(:,1);
	dblStartT = max([vecSpikeTimes1(1) min(vecEventT1)-dblUseMaxDur*5*dblJitterSize]);
	dblStopT = max(vecEventT1)+dblUseMaxDur*5*dblJitterSize;
	vecSpikeTimes1(vecSpikeTimes1 < dblStartT | vecSpikeTimes1 > dblStopT) = [];
	if numel(vecSpikeTimes1) < 3
		return;
	end
	vecEventT2 = vecEventStarts2(:,1);
	dblStartT = max([vecSpikeTimes2(1) min(vecEventT2)-dblUseMaxDur*5*dblJitterSize]);
	dblStopT = max(vecEventT2)+dblUseMaxDur*5*dblJitterSize;
	vecSpikeTimes2(vecSpikeTimes2 < dblStartT | vecSpikeTimes2 > dblStopT) = [];
	if numel(vecSpikeTimes2) < 3
		return;
	end
	
	%% build pseudo data, stitching stimulus periods
	[vecPseudoSpikeTimes1,vecPseudoEventT1] = getPseudoSpikeVectors(vecSpikeTimes1,vecEventT1,dblUseMaxDur);
	[vecPseudoSpikeTimes2,vecPseudoEventT2] = getPseudoSpikeVectors(vecSpikeTimes2,vecEventT2,dblUseMaxDur);
	
	%% combine spikes
	%get reference vector
	vecSpikeT1 = getSpikeT(vecPseudoSpikeTimes1,vecPseudoEventT1,dblUseMaxDur);
	vecSpikeT2 = getSpikeT(vecPseudoSpikeTimes2,vecPseudoEventT2,dblUseMaxDur);
	vecSpikeT = sort(cat(1,vecSpikeT1,vecSpikeT2));
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFracLinear,vecRealFrac1,vecThisSpikeTimes1,vecRealFrac2,vecThisSpikeTimes2] = ...
		getTempOffsetDiff(vecSpikeT,vecPseudoSpikeTimes1,vecPseudoEventT1,vecPseudoSpikeTimes2,vecPseudoEventT2,dblUseMaxDur);
	if numel(vecRealDiff) < 3
		return
	end
	vecRealDiff = vecRealDiff - mean(vecRealDiff);
	
	%% run bootstraps; try parallel, otherwise run normal loop
	cellRandDiff = cell(1,intResampNum);
	vecMaxRandD = nan(1,intResampNum);
	intTrials1 = numel(vecPseudoEventT1);
	intTrials2 = numel(vecPseudoEventT2);
	%vecJitterPerTrial = dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2); %original zeta
	vecJitterPerTrial1 = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials1)'; %new
	matJitterPerTrial1 = nan(intTrials1,intResampNum);
	vecJitterPerTrial2 = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials2)'; %new
	matJitterPerTrial2 = nan(intTrials2,intResampNum);
	for intResampling=1:intResampNum
		matJitterPerTrial1(:,intResampling) = vecJitterPerTrial1(randperm(numel(vecJitterPerTrial1)));
		matJitterPerTrial2(:,intResampling) = vecJitterPerTrial2(randperm(numel(vecJitterPerTrial2)));
	end
	if boolPairedTest
		matJitterPerTrial2 = matJitterPerTrial1;
	end
	try
		parfor intResampling=1:intResampNum
			%% get random subsample
			vecPseudoRandEventT1 = vecPseudoEventT1 + matJitterPerTrial1(:,intResampling);
			vecPseudoRandEventT2 = vecPseudoEventT2 + matJitterPerTrial2(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTempOffsetDiff(vecSpikeT,vecPseudoSpikeTimes1,vecPseudoRandEventT1,vecPseudoSpikeTimes2,vecPseudoRandEventT2,dblUseMaxDur);
			
			%assign data
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	catch
		for intResampling=1:intResampNum
			%% get random subsample
			vecPseudoRandEventT1 = vecPseudoEventT1 + matJitterPerTrial1(:,intResampling);
			vecPseudoRandEventT2 = vecPseudoEventT2 + matJitterPerTrial2(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTempOffsetDiff(vecSpikeT,vecPseudoSpikeTimes1,vecPseudoRandEventT1,vecPseudoSpikeTimes2,vecPseudoRandEventT2,dblUseMaxDur);
			
			%assign data
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	end
	
	%% calculate measure of effect size (for equal n, d' equals Cohen's d)
	%find highest peak and retrieve value
	dblRandMu = mean(vecMaxRandD);
	dblRandVar = var(vecMaxRandD);
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	
	if boolDirectQuantile
		%calculate statistical significance using empirical quantiles
		%define p-value
		dblZetaP = 1 - (sum(dblMaxD>vecMaxRandD)/(1+numel(vecMaxRandD)));
		
		%transform to output z-score
		dblZETA = -norminv(dblZetaP/2);
	else
		%calculate statistical significance using Gumbel distribution
		[dblZetaP,dblZETA] = getGumbel(dblRandMu,dblRandVar,dblMaxD);
		%fprintf('Pre-correction d=%.3f,post-correction z=%.3f (p=%.3f)\n',dblD,dblZETA,dblP);
	end
end

