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
	[vecThisSpikeTimes,ia,ic] = unique(vecSpikesInTrial);
	%introduce minimum jitter to identical spikes
	vecNotUnique = vecSpikesInTrial(ia(diff(ia)>1));
	if ~isempty(vecNotUnique)
		dblUniqueOffset = max(eps(vecSpikesInTrial));
		for intNotUnique=1:numel(vecNotUnique)
			vecIdx = find(vecNotUnique(intNotUnique)==vecSpikesInTrial);
			vecSpikesInTrial(vecIdx) = vecSpikesInTrial(vecIdx) + dblUniqueOffset*((1:numel(vecIdx))'-numel(vecIdx)/2);
		end
		[vecThisSpikeTimes,ia,ic] = unique(vecSpikesInTrial);
	end
	vecThisSpikeFracs = linspace(1/numel(vecThisSpikeTimes),1,numel(vecThisSpikeTimes))';
	
	%get linear fractions
	vecThisFracLinear = (vecThisSpikeTimes./dblUseMaxDur);
	
	%calc difference
	vecThisDiff = vecThisSpikeFracs - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

