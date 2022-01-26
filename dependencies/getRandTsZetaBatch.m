function [vecMaxRandD,cellRandDiff] = getRandTsZetaBatch(vecPseudoT,vecPseudoTrace,vecEventT,dblSamplingInterval,dblUseMaxDur,dblJitterSize,intBatchSize,boolUseParallel)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%calculate jitter positions
	vecEventT = vecEventT(:);
	cellRandDiff = cell(1,intBatchSize);
	vecMaxRandD = nan(1,intBatchSize);
	intTrials = numel(vecEventT);
	%vecJitterPerTrial = dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2); %original zeta
	vecJitterPerTrial = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials)'; %new
	matJitterPerTrial = nan(intTrials,intBatchSize);
	for intResampling=1:intBatchSize
		matJitterPerTrial(:,intResampling) = vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
	end
	if boolUseParallel
		parfor intResampling=1:intBatchSize
			%% get random subsample
			vecStimUseOnTime = vecEventT + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTraceOffset(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,dblSamplingInterval,dblUseMaxDur);
			
			%assign data
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	else
		for intResampling=1:intBatchSize
			%% get random subsample
			vecStimUseOnTime = vecEventT + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTraceOffset(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,dblSamplingInterval,dblUseMaxDur);
			
			%assign data
			cellRandDiff{intResampling} = vecRandDiff - mean(vecRandDiff);
			vecMaxRandD(intResampling) = max(abs(cellRandDiff{intResampling}));
		end
	end
end

