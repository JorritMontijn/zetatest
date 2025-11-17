function [latency,sLatenzy2] = latenzy2(spikeTimes1,eventTimes1,spikeTimes2,eventTimes2,useDur,resampNum,peakAlpha,useParPool,useDirectQuant,restrictNeg,makePlots)
% get latency of spiking difference between two conditions, syntax:
% [latency,sLatenzy2] = latenzy2(spikeTimes1,eventTimes1,spikeTimes2,eventTimes2,useDur,resampNum,minPeakZ,useParPool,useDirectQuant,restrictNeg,makePlots)  
%   inputs:
%   - spikeTimes1: [S x 1]: spike times condition 1 (s)
%   - eventTimes1: [T x 1]: event times condition 1 (s)
%   - spikeTimes2: [S x 1]: spike times condition 2 (s) (default: spikeTimes2 = spikeTimes1);
%   - eventTimes2: [T x 1]: event times condition 2 (s)
%   - useDur: scalar or [N x 2], time to include after/around event times (s)
%       - if a scalar is provided, it is interpreted as the duration after each event, with the start time set automatically to zero
%       - default: [0, min(diff(eventtimes))], i.e., from the event time to the minimum interval between events%   - resampNum: integer, number of resamples (default: 100)
%   - jitterSize: scalar, temporal jitter window relative to useDur (s) (default: 2)
%   - peakAlpha: scalar, significance threshold (default: 0.05)
%   - doStitch: boolean flag, perform data stitching, highly recommended! (default: true)
%   - useParPool: boolean flag, use parallel pool for resamples (default: true, but only when parallel pool is already active!)
%   - useDirectQuant: boolean flag, use the empirical null-distribution rather than the Gumbel approximation (default: false)
%   - restrictNeg: boolean flag, restrict negative latencies (default: true)
%   - makePlots: integer, plotting switch (0=none, 1=raster+traces, 2=traces only, default: 0)
%
%   *alternative input for latenzy2() is possible:
%   - spikeTimes1 and spikeTimes2 should be cell arrays, where every cell contains the aligned(!) spikes for event.
%   - set eventTimes1 and eventTimes 2 to [].
%   - set useDur to match (or use a smaller window) (default: [0    max(spikeTimes)
%
%   outputs:
%   - latency: response latency (s) (NaN when no latency could be estimated)
%   - sLatenzy2: structure with fields:
%       - latency: response latency (s)
%       - peakTimes: detected peak times, one per iter (s)
%       - peakVals: detected peak values, one per iter
%       - realFrac: see plotting function for details
%       - tempDiffUnSub: idem
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
%   v0.9 - 18 February 2025
%   - created by Robin Haak
%   v1.0 - 30 June 2025

%% prep
%check for alternative input
altInput = false;
if isempty(eventTimes1) && isempty(eventTimes2) && iscell(spikeTimes1) && iscell(spikeTimes2)
    altInput = true;
    %ensure correct orientation
    spikeTimes1 = spikeTimes1(:);
    spikeTimes2 = spikeTimes2(:);
else
    %ensure correct orientation
    spikeTimes1 = spikeTimes1(:);
    if isempty(spikeTimes2)
        spikeTimes2 = spikeTimes1;
    else
        spikeTimes2 = spikeTimes2(:);
    end
    eventTimes1 = eventTimes1(:);
    eventTimes2 = eventTimes2(:);
end

%get useDur
if ~exist('useDur','var') || isempty(useDur)
    if altInput
        useDur = min([max(vertcat(spikeTimes1{:})) max(vertcat(spikeTimes2{:}))]);
    else
        eventTimes1 = sort(eventTimes1);
        eventTimes2 = sort(eventTimes2);
        useDur = min([min(diff(eventTimes1)) min(diff(eventTimes2))]);
    end
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
    resampNum = 250;
end

%get peakAlpha
if ~exist('peakAlpha','var') || isempty(peakAlpha)
    peakAlpha = 0.05;
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

%allowNegative
if ~exist('restrictNeg','var') || isempty(restrictNeg)
    restrictNeg = false;
end

%get makePlots
if ~exist('makePlots','var') || isempty(makePlots)
    makePlots = 0;
end

%enable warning for late
giveLateWarn = false;

%% MAIN
%pre-allocate
%#ok<*AGROW>
%#ok<*UNRCH>
latency = nan;
peakTimesAgg = [];
peakValsAgg = [];
realFracAgg = {};
tempDiffUnSubAgg = {};
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
sLatenzy2 = struct;

%check if negative latencies are restricted
minLatency = useDur(1);
if restrictNeg
    minLatency = 0;
end

%check if we can use the fast interpolation
useFastInterp = false;
try
    vecTest = lininterp1f([0;0.25;1],[0;0.5;1],[-1;0;0.1;0.5;1],nan);
    if isnan(vecTest(1)) && all(vecTest(2:5)==[0 0.2 2/3 1])
        useFastInterp = true;
    else
        error('lininterp1f gives incorrect output');
    end
catch
    try
        mex('lininterp1f.c');
        vecTest = lininterp1f([0;0.25;1],[0;0.5;1],[-1;0;0.1;0.5;1],nan);
        if isnan(vecTest(1)) && all(vecTest(2:5)==[0 0.2 2/3 1])
            useFastInterp = true;
        else
            error('lininterp1f gives incorrect output');
        end
    catch
        useParPool = false;
    end
end

%run
while doContinue
    thisIter = thisIter+1;

    if altInput
        %keep spikes in window
        spikesPerEvent1 = cellfun(@(x) x(x>thisMaxDur(1) & x<thisMaxDur(2)),spikeTimes1,'UniformOutput',false);
        spikesPerEvent2 = cellfun(@(x) x(x>thisMaxDur(1) & x<thisMaxDur(2)),spikeTimes2,'UniformOutput',false);
    else
        %get spikes per event, for both conditions
        [~,spikesPerEvent1] = getRelSpikeTimes(spikeTimes1,eventTimes1,thisMaxDur);
        [~,spikesPerEvent2] = getRelSpikeTimes(spikeTimes2,eventTimes2,thisMaxDur);
    end

    %get temporal difference
    [realDiff,realTime,spikeFrac1,~,spikeFrac2,~,tempDiffUnSub,fracLinear] = ...
        calcTempDiff2(spikesPerEvent1,spikesPerEvent2,thisMaxDur,useFastInterp);
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
    [peaksRand,randDiff,randTime] = runSwapBootstraps(spikesPerEvent1,spikesPerEvent2,...
        thisMaxDur,resampNum,useParPool,useFastInterp);
    meanRandDiff = cellfun(@(x)mean(x),randDiff);
    peaksRandSub = peaksRand-meanRandDiff;

    %compute significance
    [pValPeak,peakZ] = computePval(abs(realPeakSub),peaksRandSub(~isnan(peaksRandSub)),useDirectQuant);

    %store
    if ~isnan(realPeakT)
        peakValsAgg(1,thisIter) = realMaxD;
        peakTimesAgg(1,thisIter) = realPeakT;
        realFracAgg{1,thisIter} = spikeFrac1;
        realFracAgg{2,thisIter} = spikeFrac2;
        tempDiffUnSubAgg{1,thisIter} = tempDiffUnSub;
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
        warning('Estimated latency is quite late in the window (>T/2), consider plotting and/or adjusting window');
    end
else
    return
end

%build output
sLatenzy2.latency = latency;
sLatenzy2.peakTimes = peakTimesAgg;
sLatenzy2.peakVals = peakValsAgg;
sLatenzy2.realFrac = realFracAgg;
sLatenzy2.diffUnSub = tempDiffUnSubAgg;
sLatenzy2.fracLin = fracLinAgg;
sLatenzy2.realDiff = realDiffAgg;
sLatenzy2.realTime = realTimeAgg;
sLatenzy2.meanRealDiff = meanRealDiffAgg;
sLatenzy2.randDiff = randDiffAgg;
sLatenzy2.randTime = randTimeAgg;
sLatenzy2.meanRandDiff = meanRandDiffAgg;
sLatenzy2.pValsPeak = pValPeakAgg;
sLatenzy2.peakZ = peakZAgg;
sLatenzy2.latenzyIdx = peakTimesAgg==latency;

%plot, optional
if makePlots>0
    if altInput,makePlots=9;end
    sLatenzy2.figHandles = makeLatenzy2Figs(sLatenzy2,spikeTimes1,eventTimes1,spikeTimes2,eventTimes2,useDur,makePlots);
end

end
