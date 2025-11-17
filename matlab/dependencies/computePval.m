function [pVals,Z] = computePval(maxD,maxRandD,useDirectQuant)
% compute p-values and z-scores for statistical significance testing, syntax:
%   [pVals,Z] = computePval(maxD,maxRandD,useDirectQuant)
%
% inputs:
%   - maxD: [N x 1] observed maximum values.
%   - maxRandD: [M x 1] randomized maximum values (used as null distribution).
%   - useDirectQuant: boolean, use empirical quantiles (default: false).
%
% outputs:
%   - pVals: p-values for the observed values
%   - Z: Z-scores for the observed values
%
% description:
%   this function computes p-values and z-scores for a set of observed
%   values (`maxD`) compared to a randomized null distribution (`maxRandD`), 
%   supports two methods:
%     1. empirical quantile-based p-values (if `useDirectQuant = true`).
%     2. Gumbel distribution fitting for statistical testing (default).
%
% history:
%   6 January 2025 - v0.9
%   - created by Robin Haak based on code by Jorrit Montijn
%   v1.0 - 30 June 2025

%% check inputs
if ~exist('useDirectQuant','var') || isempty(useDirectQuant)
    useDirectQuant = false;
end

%ensure `maxRandD` is sorted and unique
maxRandD = sort(unique(maxRandD));
if isempty(maxRandD)
    error('computePval:EmptyMaxRandD', 'Input maxRandD must not be empty.');
end

%% calculate significance
if useDirectQuant
    %empirical quantile-based p-values
    pVals = nan(size(maxD));
    for i = 1:numel(maxD)
        %handle edge cases
        if maxD(i) < min(maxRandD) || isnan(maxD(i))
            pVals(i) = 1; % p-value of 1 for low or invalid values
        elseif maxD(i) > max(maxRandD) || isinf(maxD(i))
            pVals(i) = 1/(1+numel(maxRandD)); % small p-value for high values
        else
            %interpolate empirical quantile
            value = interp1(maxRandD, 1:numel(maxRandD),maxD(i),'linear','extrap');
            pVals(i) = 1-(value/(1+numel(maxRandD)));
        end
    end
    
    %transform to z-scores
    Z = -norminv(pVals / 2);
else
    %Gumbel distribution-based p-values and z-scores
    randMu = mean(maxRandD);
    randVar = var(maxRandD);
    [pVals,Z] = getGumbel(randMu,randVar,maxD);
end
end
