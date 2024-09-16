%testZetaTest Test for the ZETA-test implementation
%
% To run use runtests('testZetatest')
%
% 2024, Alexander Heimel, based on runExampleZETA

%% zetatest_default
disp('zetatest_default, expected to take about 2 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
rng(1,'twister'); 
dblZetaP = zetatest(sLoad.vecSpikeTimes1,sLoad.vecStimulusStartTimes);
assert(abs(dblZetaP - 7.702804163978172e-05)<1E-6)

%% zetatest_specified
disp('zetatest_specified, expected to take about 2 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
rng(1,'twister'); 
intPlot = 0;
vecRestrictRange = [0 inf];
boolDirectQuantile = false; 
boolStitch = true; 
[dblZetaP,sZETA,sRate,sLatencies] = zetatest(sLoad.vecSpikeTimes1,sLoad.matEventTimes,...
    sLoad.dblUseMaxDur,sLoad.intResampNum,intPlot,vecRestrictRange,boolDirectQuantile,...
    sLoad.dblJitterSize,boolStitch);
assert(abs(dblZetaP - 1.353659556455611e-04)<1E-6)

%% zetatstest_default
disp('zetatstest_default, expected to take about 6 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
rng(1,'twister');
dblZetaP = zetatstest(sLoad.vecTimestamps, sLoad.vecData1, sLoad.matEventTimesTs(:,1),[],sLoad.intResampNum);
assert(abs(dblZetaP-0.027278506931302)<1E-6)

%% zetatstest_specified
disp('zetatstest_specified, expected to take about 6 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
rng(1,'twister');
intPlot = 0;
vecRestrictRange = [0 inf];
boolDirectQuantile = false; 
boolStitch = true; 
[dblZetaP, sZeta] = zetatstest(sLoad.vecTimestamps, sLoad.vecData1, sLoad.matEventTimesTs,...
    sLoad.dblUseMaxDur, sLoad.intResampNum, intPlot, boolDirectQuantile,...
    sLoad.dblJitterSize,boolStitch);
assert(abs(dblZetaP - 0.027276318742224)<1E-6)


%% zetatest2_neurons
disp('zetatest2_two_neurons, expected to take about 0.7 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
rng(1,'twister');
intTrials = 240; %that's already more than enough
intResampNum = [];
intPlot = 0;
[dblZetaP,sZeta] = zetatest2(sLoad.vecSpikeTimes1,sLoad.matEventTimes(1:intTrials,:),...
    sLoad.vecSpikeTimes2,sLoad.matEventTimes(1:intTrials,:),...
    sLoad.dblUseMaxDur,[],intPlot);
assert(abs(dblZetaP - 8.685750946035853e-06)<1E-6)

%% zetatest2_stimuli
disp('zetatest2_two_stimuli, expected to take about 0.1 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
vecTrials1 = sLoad.vecStimulusOrientation==0;
vecTrials2 = sLoad.vecStimulusOrientation==90;
intResampNum = [];
intPlot = 0;
rng(1,'twister');
[dblZetaP,sZeta] = zetatest2(sLoad.vecSpikeTimes1,sLoad.matEventTimes(vecTrials1,:),...
    sLoad.vecSpikeTimes1,sLoad.matEventTimes(vecTrials2,:),sLoad.dblUseMaxDur,intResampNum,intPlot);
assert(abs(dblZetaP - 0.010218584062786)<1E-6)

%% zetatstest2_neurons
disp('zetatstest2_neurons, expected to take about 0.2 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
rng(1,'twister');
intTrials = 240; %that's already more than enough
intResampNum = [];
intPlot = 0;
[dblZetaP,sZeta] = zetatstest2(sLoad.vecTimestamps,sLoad.vecData1,sLoad.matEventTimesTs(1:intTrials,:),...
    sLoad.vecTimestamps,sLoad.vecData2,sLoad.matEventTimesTs(1:intTrials,:),...
    sLoad.dblUseMaxDur,intResampNum,intPlot);
assert(abs(dblZetaP - 7.201222648300920e-06)<1E-6)

%% zetatstest2_stimuli
disp('zetatstest2_stimuli, expected to take about 0.04 s')
pth = fileparts(which('testZetatest'));
sLoad = load(fullfile(pth,'testZetaTestData.mat'));
rng(1,'twister');
intResampNum = [];
intPlot = 0;
vecTrials1 = sLoad.vecStimulusOrientation(1:480)==0;
vecTrials2 = sLoad.vecStimulusOrientation(1:480)==90;
[dblZetaP,sZeta] = zetatstest2(sLoad.vecTimestamps,sLoad.vecData1,sLoad.matEventTimesTs(vecTrials1,:),...
    sLoad.vecTimestamps,sLoad.vecData1,sLoad.matEventTimesTs(vecTrials2,:),...
    sLoad.dblUseMaxDur,intResampNum,intPlot);
assert(abs(dblZetaP - 0.033518074062296)<1E-6)

% 
% sLoad = load('testZetaTestData.mat');
% vecSpikeTimes1 = sLoad.vecSpikeTimes1;
% vecSpikeTimes2 = sLoad.vecSpikeTimes2;
% vecStimulusOrientation = sLoad.vecStimulusOrientation;
% vecStimulusStartTimes = sLoad.vecStimulusStartTimes;
% vecStimulusStopTimes= sLoad.vecStimulusStopTimes;
% vecData1 = sLoad.vecData1;
% vecData2 = sLoad.vecData2;
% vecTimestamps = sLoad.vecTimestamps;
% intUseTrialNum = sLoad.intUseTrialNum;
% matEventTimesTs = sLoad.matEventTimesTs;
% intResampNum = sLoad.intResampNum;
% 
% 
% save('testZetaTestData.mat','vecSpikeTimes1','vecSpikeTimes2',...
%     'vecStimulusOrientation','vecStimulusStartTimes','vecStimulusStopTimes', ...
%     'matEventTimes','dblUseMaxDur','intResampNum','dblJitterSize',...
%     'vecData1','vecData2','vecTimestamps','matEventTimesTs'); %#ok<UNRCH>
