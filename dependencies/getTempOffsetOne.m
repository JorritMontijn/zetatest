function [vecThisDiff,vecThisSpikeFracs,vecThisFracLinear,vecThisSpikeTimes] = ...
			getTempOffsetOne(vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur)
	%getTempOffset Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisSpikeFracs,vecThisFracLinear,vecThisSpikeTimes] = ...
	%	getTempOffsetOne(vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur)
	%
	%This is a subfunction for zetatest()/calcZetaOne().
	
	%% get temp diff vector
	%pre-allocate
	vecSpikesInTrial = getSpikeT(vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
	vecThisSpikeTimes = getUniqueSpikes(vecSpikesInTrial);
	
	%flat fracs
	vecThisSpikeFracs = linspace(1/numel(vecThisSpikeTimes),1,numel(vecThisSpikeTimes))';
	
	%get linear fractions
	vecThisFracLinear = (vecThisSpikeTimes./dblUseMaxDur);
	
	%calc difference
	vecThisDiff = vecThisSpikeFracs - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

