function [latency,sLatenzy] = latenzy(spikeTimes,eventTimes,useDur,resampNum,jitterSize,peakAlpha,doStitch,useParPool,useDirectQuant,restrictNeg,makePlots)
% get event-related spiking latency, syntax:
%   [latency,sLatenzy] = latenzy(spikeTimes,eventTimes,useDur,resampNum,jitterSize,minPeakZ,doStitch,useParPool,useDirectQuant,restrictNeg,makePlots)
%   
%   inputs:
%   - spikeTimes: [S x 1]: spike times (s)
%   - eventTimes: [T x 1]: event times (s)
%   - useDur: scalar or [N x 2], time to include after/around event times (s)
%       - if a scalar is provided, it is interpreted as the duration after each event, with the start time set automatically to zero
%       - default: [0, min(diff(eventtimes))], i.e., from the event time to the minimum interval between events
%   - resampNum: integer, number of resamples (default: 100)
%   - jitterSize: scalar, temporal jitter window relative to useDur (s) (default: 2)
%   - peakAlpha: scalar, significance threshold (default: 0.05)
%   - doStitch: boolean flag, perform data stitching, highly recommended! (default: true)
%   - useParPool: boolean flag, use parallel pool for resamples (default: true, but only when parallel pool is already active!)
%   - useDirectQuant: boolean flag, use the empirical null-distribution rather than the Gumbel approximation (default: false)
%   - restrictNeg: boolean flag, restrict negative latencies (default: true)
%   - makePlots: integer, plotting switch (0=none, 1=raster+traces, 2=traces only, default: 0)
%
%   outputs:
%   - latency: response latency (s) (NaN when no latency could be estimated)
%   - sLatenzy: structure with fields:
%       - latency: response latency (s)
%       - peakTimes: detected peak times, one per iter (s)
%       - peakVals: detected peak values, one per iter
%       - realFrac: see plotting function for details
%       - fracLin: idem
%       - realDiff: idem
%       - realTime: idem
%       - meanRealDiff: idem
%       - randDiff: idem
%       - randTime: idem
%       - meanRandDiff: idem
%       - pVals: p-values for the observed peak maxima
%       - peakZ: two-tailed z-scores corresponding to the p-values
%       - latenzyIdx: use to index arrays above
%       - handleFigs: figure handles
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025

%% prep
%ensure correct orientation
spikeTimes = spikeTimes(:);
eventTimes = eventTimes(:);

%get useDur
if ~exist('useDur','var') || isempty(useDur)
    eventTimes = sort(eventTimes);
    useDur = min(diff(eventTimes));
end
if isscalar(useDur)
    useDur = sort([0 useDur]);
elseif numel(useDur)~=2
    error([mfilename ':WrongMaxDurInput'],'useDur must be a scalar or a two-element array');
end

%check useDur
if useDur(2)>0
    assert(useDur(1)<=0,[mfilename ':WrongMaxDurInput'],...
        sprintf('When useDur(2) > 0, useDur(1) must be a negative scalar or 0, you requested [%.3f %.3f]',useDur(1),useDur(2)));
elseif useDur(2)==0
    assert(useDur(1)<0,[mfilename ':WrongMaxDurInput'],...
        sprintf('When useDur(2) is 0, useDur(1) must be a negative scalar, you requested [%.3f %.3f]',useDur(1),useDur(2)));
elseif useDur(2)<0
        error([mfilename ':WrongMaxDurInput'],'useDur(2) cannot be negative when useDur(1) is negative!');
end

%get resampNum
if ~exist('resampNum','var') || isempty(resampNum)
    resampNum = 100;
end

%get jitterSize
if ~exist('jitterSize','var') || isempty(jitterSize)
    jitterSize = 2;
end
assert(jitterSize>0,[mfilename ':WrongJitterInput'], ...
    sprintf('jitterSize must be a positive scalar, you requested %.3f',jitterSize));

%get peakAlpha
if ~exist('peakAlpha','var') || isempty(peakAlpha)
    peakAlpha = 0.05;
end

%get doStitch
if ~exist('doStitch','var') || isempty(doStitch)
    doStitch = true;
end

%get useParPool
if ~exist('useParPool','var') || isempty(useParPool)
    try
        objPool = gcp('nocreate'); %get current parallel pool (no creation)
        if isempty(objPool) || ~isprop(objPool,'NumWorkers') || objPool.NumWorkers < 4
            useParPool = false; 
        else
            useParPool = true;
        end
    catch
        useParPool = false;
    end
end

%useDirectQuant
if ~exist('useDirectQuant','var') || isempty(useDirectQuant)
    useDirectQuant = false;
end

%restricNeg
if ~exist('restrictNeg','var') || isempty(restrictNeg)
    restrictNeg = false;
end

%get makePlots
if ~exist('makePlots','var') || isempty(makePlots)
    makePlots = 0;
end

%enable warning when estimates are > T/2
giveLateWarn = false;

%% MAIN
%pre-allocate
%#ok<*AGROW>
%#ok<*UNRCH>
latency = nan;
peakTimesAgg = [];
peakValsAgg = [];
realFracAgg = {};
fracLinAgg = {};
realDiffAgg = {};
realTimeAgg = {};
meanRealDiffAgg = [];
randDiffAgg = {};
randTimeAgg = {};
meanRandDiffAgg = [];
pValPeakAgg = [];
peakZAgg = [];
keepPeaks = [];
thisMaxDur = useDur;
doContinue = true;
thisIter = 0;
sLatenzy = struct;

%check if negative latencies are restricted
minLatency = useDur(1);
if restrictNeg
    minLatency = 0;
end

%run
while doContinue
    thisIter = thisIter+1;

    %perform data-stitching
    if doStitch
        discardEdges = true;
        [pseudoSpikeTimes,pseudoEventTimes] = getPseudoTimes(spikeTimes,eventTimes,thisMaxDur,discardEdges);
    else
        pseudoSpikeTimes = spikeTimes;
        pseudoEventTimes = eventTimes;
    end

    %get temporal deviation
    [realDiff,realTime,spikeFrac,fracLinear] = calcTempDiff(pseudoSpikeTimes,pseudoEventTimes,thisMaxDur);
    if numel(realDiff) < 3
        return
    end
    
    %get largest deviation
    [maxDiff,maxIdx] = max(realDiff);
    [minDiff,minIdx] = min(realDiff);
    if abs(minDiff) >= abs(maxDiff)
        realMaxD = minDiff;
        realMaxIdx = minIdx;
    else
        realMaxD = maxDiff;
        realMaxIdx = maxIdx;
    end
    realPeakT = realTime(realMaxIdx);
    realPeakSub = realMaxD-mean(realDiff);

    %run bootstraps
    [peaksRand,randDiff,randTime] = runJitterBootstraps(pseudoSpikeTimes,pseudoEventTimes,...
        thisMaxDur,resampNum,jitterSize,useParPool);
    meanRandDiff = cellfun(@(x)mean(x),randDiff);
    peaksRandSub = peaksRand-meanRandDiff;

    %compute significance
    [pValPeak,peakZ] = computePval(abs(realPeakSub),peaksRandSub(~isnan(peaksRandSub)),useDirectQuant);

    %store
    if ~isnan(realPeakT)
        peakValsAgg(1,thisIter) = realMaxD;
        peakTimesAgg(1,thisIter) = realPeakT;
        realFracAgg{1,thisIter} = spikeFrac;
        fracLinAgg{1,thisIter} = fracLinear;
        realDiffAgg{1,thisIter} = realDiff;
        realTimeAgg{1,thisIter} = realTime;
        meanRealDiffAgg(1,thisIter) = mean(realDiff);
        randDiffAgg(:,thisIter) = randDiff;
        randTimeAgg(:,thisIter) = randTime;
        meanRandDiffAgg(:,thisIter) = meanRandDiff;
        pValPeakAgg(1,thisIter) = pValPeak;
        peakZAgg(1,thisIter) = abs(peakZ);
    end

    %check whether to continue
    if realPeakT > minLatency && pValPeak < peakAlpha
        keepPeaks(thisIter) = true;
        thisMaxDur(2) = realPeakT;
    else
        doContinue = false;
        keepPeaks(thisIter) = false;
    end
end

%get latency
thesePeakTimes = peakTimesAgg(logical(keepPeaks));
if ~isempty(thesePeakTimes)
    latency = thesePeakTimes(end);

    %warning
    if latency > (useDur(1)+sum(abs(useDur))/2) && giveLateWarn
        warning('Estimated latency is quite late in the window (>T/2). Consider plotting for visual check and/or adjusting window');
    end
else
    return
end

%build output
sLatenzy.latency = latency;
sLatenzy.peakTimes = peakTimesAgg;
sLatenzy.peakVals = peakValsAgg;
sLatenzy.realFrac = realFracAgg;
sLatenzy.fracLin = fracLinAgg;
sLatenzy.realDiff = realDiffAgg;
sLatenzy.realTime = realTimeAgg;
sLatenzy.meanRealDiff = meanRealDiffAgg;
sLatenzy.randDiff = randDiffAgg;
sLatenzy.randTime = randTimeAgg;
sLatenzy.meanRandDiff = meanRandDiffAgg;
sLatenzy.pValsPeak = pValPeakAgg;
sLatenzy.peakZ = peakZAgg;
sLatenzy.latenzyIdx = peakTimesAgg==latency;

%plot, optional
if makePlots>0, sLatenzy.figHandles = makeLatenzyFigs(sLatenzy,spikeTimes,eventTimes,useDur,makePlots); end

end
