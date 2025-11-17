# *latenZy* — Preprocessing guidelines 
**Contents:** 
- Preparing trial-aligned data for ***latenZy***
- Using ***latenZy2*** with trial-aligned spike data

## Preparing trial-aligned data for ***`latenZy`***
***`latenZy`*** requires spike and event times as continuous, absolute timestamps. If your data is trial-aligned (e.g., spikes relative to stimulus onset) and you do not have the original event times, you can simulate them by assigning large, fixed offsets between repetitions (e.g., 100s apart):

**Python example:**
```python
import numpy as np

#trial-aligned spike times for 5 trials (relative to stimulus onset)
aligned_spikes = [
    [-0.3, 0.2, 0.3, 0.7],
    [-0.4, 0.1, 0.5],
    [-0.6, 0.3, 0.6],
    [-0.2, 0.1, 0.4, 0.9],
    [-0.5, 0.6]
]

#simulated absolute event (stimulus) times, spaced 100s apart
event_times = np.arange(100, 100 * len(aligned_spikes) + 1, 100)

#offset each trial's spikes by its simulated global event time
new_spikes = []
for i, t_event in enumerate(event_times):
    new_spikes.append(np.array(aligned_spikes[i]) + t_event)

#flatten all spikes into a global spike time array
spike_times = np.concatenate(new_spikes)
```

**MATLAB example:**
```matlab
%trial-aligned spike times for 5 trials (relative to stimulus onset)
alignedSpikes = {
    [-0.3, 0.2, 0.3, 0.7];
    [-0.4, 0.1, 0.5];
    [-0.6, 0.3, 0.6];
    [-0.2, 0.1, 0.4, 0.9];
    [-0.5, 0.6]
};

%simulated absolute event (stimulus) times, spaced 100s apart
eventTimes = 100:100:(100 * numel(alignedSpikes));

%offset each trial's spikes by its simulated global event time
newSpikes = cell(size(alignedSpikes));
for i = 1:numel(alignedSpikes)
    newSpikes{i} = alignedSpikes{i} + eventTimes(i);
end

%flatten all spikes into a global spike time array
spikeTimes = [newSpikes{:}];
```

This produces pseudo-global spike and event times compatible with *`latenZy`*. 

In the first step of the algorithm, data is stitched across repetitions by removing spikes outside the event window `use_dur`/`useDur`. The excluded intervals between event repetiations are substracted from all subsequent times, creating a continuous timeline of only event-related activity for statistics. Make sure that `do_stitch`/`doStitch` is set to **True** (it is by default).
> ⚠️ **Important:** Make sure `use_dur`/`useDur` does **not** include periods without spikes, as silent intervals distort the stitched data and bias latency estimates.


## Using ***`latenZy2`*** with trial-aligned spike data
***`latenZy2`*** can directly analyze trial-aligned spike data without needing event times (because the statistics do not require jittering event times).
To do this:
- Pass the spike times for each trial in lists (Python) or cell arrays (MATLAB).
- Set event time inputs to empty lists/arrays ([ ]).
- Define the event window relative to stimulus (e.g., [0, 1] seconds).

**Python example:**
```python
t, s_latenzy2 = latenzy2(spike_times1, [], spiketimes2, [], use_dur=[0, 1])
```

**MATLAB example:**
```matlab
[t, sLatenzy2] = latenzy2(spikeTimes1, [], spikeTimes2, [], [0 1]);
```


