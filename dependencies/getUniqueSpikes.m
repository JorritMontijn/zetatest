function vecSpikeTimes = getUniqueSpikes(vecSpikeTimes)
	%introduce minimum jitter to identical spikes
	vecSpikeTimes = sort(vecSpikeTimes);
	dblUniqueOffset = max(eps(vecSpikeTimes));
	indDuplicates = [false;diff(vecSpikeTimes)<dblUniqueOffset];
	vecNotUnique = vecSpikeTimes(indDuplicates);
	if ~isempty(vecNotUnique)
		vecJitter = cat(1,1+4*rand([numel(vecNotUnique),1]),-1-4*rand([numel(vecNotUnique),1]));
		vecJitter = vecJitter(randperm(numel(vecJitter),numel(vecNotUnique)));
		vecSpikeTimes(indDuplicates) = vecSpikeTimes(indDuplicates) + vecJitter;
	end
end
%{
function vecUniqueSpikeTimes = getUniqueSpikes(vecSpikeTimes)
	%introduce minimum jitter to identical spikes
	vecSpikeTimes = sort(vecSpikeTimes);
	[vecUniqueSpikeTimes,ia,ic] = unique(vecSpikeTimes);
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

runZeta2Steinmetz_202311091263.783 s1.270 sanovan555143.708 s0.272 sanovan>dofit2220140.656 s140.610 szetatest277675.790 s0.161 scalcZetaTwo77675.016 s6.631 sgetTempOffsetTwo15711767.864 s6.159 sopthist55541.771 s0.109 sopthist>findC253641.326 s0.079 ssshist2298141.234 s7.556 sinterp131423425.474 s21.844 ssshist>sscost68943024.579 s15.695 sgetUniqueSpikes31423423.540 s8.357 sunique46788415.349 s6.253 s

%}