function vecSpikeTimes = getUniqueSpikes(vecSpikeTimes)
	%introduce minimum jitter to identical spikes
	vecSpikeTimes = sort(vecSpikeTimes);
	dblUniqueOffset = max(eps(vecSpikeTimes));
	indDuplicates = [false;diff(vecSpikeTimes)<dblUniqueOffset];
	vecNotUnique = vecSpikeTimes(indDuplicates);
	if ~isempty(vecNotUnique)
		vecJitter = cat(1,1+4*rand([numel(vecNotUnique),1]),-1-4*rand([numel(vecNotUnique),1]));
		vecJitter = dblUniqueOffset*vecJitter(randperm(numel(vecJitter),numel(vecNotUnique)));
		vecSpikeTimes(indDuplicates) = vecSpikeTimes(indDuplicates) + vecJitter;
	end
end