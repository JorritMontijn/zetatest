function spikeTimes = getDistinctSpikes(spikeTimes)
%introduce minimal jitter to repeating spike times (if any)

% 17 October 2025
% - introduced my_randperm() to ensure consistency between MATLAB and Python

spikeTimes = sort(spikeTimes);
uniqueOffset = max(eps(spikeTimes));
idxRepeat = [false;diff(spikeTimes)<uniqueOffset];
while any(idxRepeat)
    notUnique = spikeTimes(idxRepeat);
    addJitter = cat(1,1+9*rand([numel(notUnique),1]),-1-9*rand([numel(notUnique),1]));
    addJitter = uniqueOffset*addJitter(my_randperm(numel(addJitter),numel(notUnique)));
    spikeTimes(idxRepeat) = spikeTimes(idxRepeat)+addJitter;
    spikeTimes = sort(spikeTimes);
    idxRepeat = [false;diff(spikeTimes)<uniqueOffset];
end
