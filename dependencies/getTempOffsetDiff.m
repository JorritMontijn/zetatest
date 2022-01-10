function [vecThisDiff,vecThisFracLinear,vecThisFrac1,vecThisSpikeTimes1,vecThisFrac2,vecThisSpikeTimes2] = ...
			getTempOffsetDiff(vecSpikeT,vecSpikeTimes1,vecStimUseOnTime1,vecSpikeTimes2,vecStimUseOnTime2,dblUseMaxDur)
	%getTempOffset Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisFracLinear,vecThisFrac1,vecThisSpikeTimes1,vecThisFrac2,vecThisSpikeTimes2] = ...
	%	getTempOffsetDiff(vecSpikeT,vecSpikeTimes1,vecStimUseOnTime1,vecSpikeTimes2,vecStimUseOnTime2,dblUseMaxDur)
	%
	%This is a subfunction for zetatest2().
	
	%% get linear fractions
	vecThisFracLinear = (vecSpikeT./dblUseMaxDur);
	
	%% get temp diff vector 1
	%get reference vector
	vecSpikeT1 = getSpikeT(vecSpikeTimes1,vecStimUseOnTime1,dblUseMaxDur);
	vecThisSpikeTimes1 = unique(vecSpikeT1);
	
	%pre-allocate
	vecThisSpikeFracs1 = linspace(1/numel(vecThisSpikeTimes1),1,numel(vecThisSpikeTimes1))';
	vecThisFrac1 = interp1(vecThisSpikeTimes1,vecThisSpikeFracs1,vecSpikeT);
	
	%calc deviation
	vecDeviation1 = vecThisFrac1 - vecThisFracLinear;
	
	%% get temp diff vector 2
	%get reference vector
	vecSpikeT2 = getSpikeT(vecSpikeTimes2,vecStimUseOnTime2,dblUseMaxDur);
	vecThisSpikeTimes2 = unique(vecSpikeT2);
	
	%pre-allocate
	vecThisSpikeFracs2 = linspace(1/numel(vecThisSpikeTimes2),1,numel(vecThisSpikeTimes2))';
	vecThisFrac2 = interp1(vecThisSpikeTimes2,vecThisSpikeFracs2,vecSpikeT);
	
	%calc deviation
	vecDeviation2 = vecThisFrac2 - vecThisFracLinear;
	
	%% calc diff
	vecThisDiff = vecDeviation1 - vecDeviation2;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

