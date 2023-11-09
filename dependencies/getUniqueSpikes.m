function vecSpikeTimes = getUniqueSpikes(vecSpikeTimes)
	%introduce minimum jitter to identical spikes
	vecSpikeTimes = sort(vecSpikeTimes);
	dblUniqueOffset = max(eps(vecSpikeTimes));
	indDuplicates = [false;diff(vecSpikeTimes)<dblUniqueOffset];
	while any(indDuplicates)
		vecNotUnique = vecSpikeTimes(indDuplicates);
		vecJitter = cat(1,1+9*rand([numel(vecNotUnique),1]),-1-9*rand([numel(vecNotUnique),1]));
		vecJitter = dblUniqueOffset*vecJitter(randperm(numel(vecJitter),numel(vecNotUnique)));
		vecSpikeTimes(indDuplicates) = vecSpikeTimes(indDuplicates) + vecJitter;
		vecSpikeTimes = sort(vecSpikeTimes);
		indDuplicates = [false;diff(vecSpikeTimes)<dblUniqueOffset];
	end
end