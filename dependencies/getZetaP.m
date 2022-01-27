function [dblZetaP,dblZETA] = getZetaP(dblMaxD,vecMaxRandD,boolDirectQuantile)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%% calculate significance
	%find highest peak and retrieve value
	vecMaxRandD = sort(unique(vecMaxRandD));
	dblRandMu = mean(vecMaxRandD);
	dblRandVar = var(vecMaxRandD);
	
	if boolDirectQuantile
		%calculate statistical significance using empirical quantiles
		%define p-value
		dblZetaP = nan(size(dblMaxD));
		for i=1:numel(dblMaxD)
			if dblMaxD < min(vecMaxRandD) || isnan(dblMaxD)
				dblValue = 0;
			elseif dblMaxD > max(vecMaxRandD) || isinf(dblMaxD)
				dblValue = numel(vecMaxRandD);
			else
				dblValue = interp1(vecMaxRandD,1:numel(vecMaxRandD),dblMaxD);
			end
			dblZetaP(i) = 1 - (dblValue/(1+numel(vecMaxRandD)));
		end
		
		%transform to output z-score
		dblZETA = -norminv(dblZetaP/2);
	else
		%calculate statistical significance using Gumbel distribution
		[dblZetaP,dblZETA] = getGumbel(dblRandMu,dblRandVar,dblMaxD);
		%fprintf('Pre-correction d=%.3f,post-correction z=%.3f (p=%.3f)\n',dblD,dblZETA,dblP);
	end
end

