function vecThisFrac = fillnans(vecThisFrac,intSp,intT)
	indIsNan = isnan(vecThisFrac);
	if isnan(vecThisFrac(1)) || isnan(vecThisFrac(end))
		indNanDiff = diff(indIsNan);
	end
	if isnan(vecThisFrac(1))
		intLeadingNanLength1 = find(indNanDiff,1,'first');
		vecThisFrac(1:intLeadingNanLength1) = 1/intT;
	end
	if isnan(vecThisFrac(end))
		intLaggingNanLength1 = find(indNanDiff(end:-1:1),1,'first');
		vecThisFrac((end-intLaggingNanLength1+1):end) = intSp/intT;
	end
end