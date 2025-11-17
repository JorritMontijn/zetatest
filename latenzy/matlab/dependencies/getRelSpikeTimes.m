function [relSpikeTimes,spikesPerEvent] = getRelSpikeTimes(spikeTimes,eventTimes,useDur,addArtifSpikes)
% create a vector of spike times relative to event times, syntax:
%   [relSpikeTimes,spikesPerEvent] = getRelSpikeTimes(spikeTimes,eventTimes,useDur,addArtifSpks)
%   inputs:
%   - spikeTimes [S x 1]: spike times (s)
%   - eventTimes [T x 1]: event (start) times (s)
%   - useDur: scalar or [pre post], time to include after/around event times (s)
%   - addArtifSpikes: boolean, add artificial spikes at beginning and end of epoch (default: true)
%
%   outputs:
%   - relSpikeTimes: spike times relative to events (s), sorted
%   - spikesPerEvent: relative spike times per event (s), sorted
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025

%% prep
%ensure correct orientation
spikeTimes = spikeTimes(:);
eventTimes = eventTimes(:);

%check inputs
if ~exist('useDur','var') || isempty(useDur)
    eventTimes = sort(eventTimes);
    useDur = min(diff(eventTimes));
end

if isscalar(useDur), useDur = sort([0 useDur]); end
assert(useDur(2)>useDur(1),[mfilename ':WrongMaxDurInput'],...
    sprintf('The second element of useDur must be larger than the first element, you requested [%.3f %.3f]',...
    useDur(1),useDur(2)));

if ~exist('addArtifSpikes','var') || isempty(addArtifSpikes)
    addArtifSpikes = false;
end

%% compute relative spike times
spikesPerEvent = arrayfun(@(x) spikeTimes(spikeTimes > (x+useDur(1)) & spikeTimes < (x+useDur(2)))-x, ...
    eventTimes,'UniformOutput',false);

%concatenate
relSpikeTimes = sort(cell2vec(spikesPerEvent));

%% if requested, add artificial spikes to cover full epoch
if addArtifSpikes && ~isempty(relSpikeTimes)
    relSpikeTimes = unique([useDur(1); relSpikeTimes; useDur(2)]);
end

