function vecSpikeTimes = getUniqueSpikes(vecSpikeTimes)
	%introduce minimum jitter to identical spikes
	vecSpikeTimes = sort(vecSpikeTimes);
	%dblUniqueOffset = max(eps(vecSpikeTimes));
	dblUniqueOffset = eps(class(vecSpikeTimes)); % to match Python version
	indDuplicates = [false;diff(vecSpikeTimes)<dblUniqueOffset];
	while any(indDuplicates)
		vecNotUnique = vecSpikeTimes(indDuplicates);
		vecJitter = cat(1,1+9*rand([numel(vecNotUnique),1]),-1-9*rand([numel(vecNotUnique),1]));
		vecJitter = dblUniqueOffset*vecJitter(my_randperm(numel(vecJitter),numel(vecNotUnique)));
		vecSpikeTimes(indDuplicates) = vecSpikeTimes(indDuplicates) + vecJitter;
		vecSpikeTimes = sort(vecSpikeTimes);
		indDuplicates = [false;diff(vecSpikeTimes)<dblUniqueOffset];
        dblUniqueOffset = dblUniqueOffset * 2; % to avoid endless loop if vecJitter is too small
    end
end

function x = my_randperm(n,k)
% randperm introduced to make results reproducable between python and
% MATLAB implementation
[~,ind] = sort(rand(n,1));
x = ind(1:k);
end