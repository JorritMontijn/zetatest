function [dblZetaDev,dblZetaDev_Inv,intZetaDevIdx,intZetaDevIdx_Inv] = getZetaDev(vecIn)
	%getZetaDev Returns +ZETA and -ZETA
	%   [dblZetaDev,dblZetaDev_Inv,intZetaDevIdx,intZetaDevIdx_Inv] = getZetaDev(vecIn)
	
	[dummy,intZetaDevIdx]=max(abs(vecIn));
	dblZetaDev = vecIn(intZetaDevIdx);
	[dblZetaDev_Inv,intZetaDevIdx_Inv]=max(vecIn*-sign(dblZetaDev));
end

