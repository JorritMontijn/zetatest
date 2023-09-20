function [vecRefT,cellSampleAssignments] = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur)
    %getTsRefT Get reference time vector and sample assignments
    %   [vecRefT,cellSampleAssignments] = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur)

    %pre-allocate
    vecEventStartT = sort(vecEventStartT);
    intTrialNum = numel(vecEventStartT);
    intTimeNum = numel(vecTimestamps);
    if size(vecTimestamps,1)==1,vecTimestamps=vecTimestamps';end
    
    %build common timeframe
    cellRefT = cell(1,intTrialNum);
    cellSampleAssignments = cell(1,intTrialNum);
    try
        [vecStart,vecStop] = findgtentries2_mex(vecTimestamps,vecEventStartT,vecEventStartT+dblUseMaxDur);
    catch
        [vecStart,vecStop] = findgtentries2(vecTimestamps,vecEventStartT,vecEventStartT+dblUseMaxDur);
    end
    for intTrial=1:intTrialNum
        vecSelectSamples = vecStart(intTrial):vecStop(intTrial);
        cellSampleAssignments{intTrial} = vecSelectSamples;
        cellRefT{intTrial} = vecTimestamps(vecSelectSamples)-vecEventStartT(intTrial);
    end

%     %old and slow
%     for intTrial=1:intTrialNum
%         % get original times
%         dblStartT = vecEventStartT(intTrial);
%         dblStopT = dblStartT+dblUseMaxDur;
%         intStartT = max([1,findfirstgreaterthan(vecTimestamps(intStartT:end),dblStartT) + intStartT - 2]);
%         intStopT = min([intTimeNum,findfirstgreaterthan(vecTimestamps(intStopT:end),dblStopT)+ intStopT - 1]);
%         %intStartT = max([1 find(vecTimestamps > dblStartT,1) - 1]);
%         %intStopT = min([intTimeNum find(vecTimestamps > dblStopT,1)]);
%         vecSelectSamples = intStartT:intStopT;
%         cellSampleAssignments{intTrial} = vecSelectSamples;
% 
%         %% get data
%         cellRefT{intTrial} = vecTimestamps(vecSelectSamples)-dblStartT;
%     end

    %set tol
    dblSampInterval = median(diff(vecTimestamps));
    dblTol = dblSampInterval/100;
    vecRefT = uniquetol(sort(cell2vec(cellRefT)),dblTol);
end
function [vecCompIdx1,vecCompIdx2] = findgtentries2(vecIn,vecComp1,vecComp2)
    intVecNum = numel(vecIn);
    intCompNum1 = numel(vecComp1);
    intCompNum2 = numel(vecComp1);
    vecCompIdx1 = intVecNum*ones(intCompNum1,1);
    vecCompIdx2 = intVecNum*ones(intCompNum2,1);
    intCompIdx1=1;
    intCompIdx2=1;
    dblCompVal1=vecComp1(intCompIdx1);
    dblCompVal2=vecComp2(intCompIdx2);
    boolComp1Done = false;
    boolComp2Done = false;
    for intIdx=1:intVecNum
        if boolComp1Done && boolComp2Done
            return;
        end
        dblCurrVal = vecIn(intIdx);
        if ~boolComp1Done && dblCurrVal > dblCompVal1
            vecCompIdx1(intCompIdx1) = intIdx;
            intCompIdx1 = intCompIdx1 + 1;
            boolComp1Done = intCompIdx1 > intCompNum1;
            if ~boolComp1Done
                dblCompVal1 = vecComp1(intCompIdx1);
            end
        end
        if ~boolComp2Done && dblCurrVal > dblCompVal2
            vecCompIdx2(intCompIdx2) = intIdx;
            intCompIdx2 = intCompIdx2 + 1;
            boolComp2Done = intCompIdx2 > intCompNum2;
            if ~boolComp2Done
                dblCompVal2 = vecComp2(intCompIdx2);
            end
        end
    end
end
function intIdx = findfirstgreaterthan(vecIn,dblCompVal)
    for intIdx=1:numel(vecIn)
        if vecIn(intIdx) > dblCompVal
            return;
        end
    end
end
