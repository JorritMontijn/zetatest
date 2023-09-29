function vecUniqueSpikeTimes = getUniqueSpikes(vecSpikeTimes)
	[vecUniqueSpikeTimes,ia,ic] = unique(vecSpikeTimes);
	%introduce minimum jitter to identical spikes
	vecNotUnique = vecSpikeTimes(ia(diff(ia)>1));
	if ~isempty(vecNotUnique)
		dblUniqueOffset = max(eps(vecSpikeTimes));
		for intNotUnique=1:numel(vecNotUnique)
			vecIdx = find(vecNotUnique(intNotUnique)==vecSpikeTimes);
			vecSpikeTimes(vecIdx) = vecSpikeTimes(vecIdx) + dblUniqueOffset*((1:numel(vecIdx))'-numel(vecIdx)/2);
		end
		[vecUniqueSpikeTimes,ia,ic] = unique(vecSpikeTimes);
	end
end