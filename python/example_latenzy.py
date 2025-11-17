# -*- coding: utf-8 -*-
"""
example_latenzy

Example code using the latenzy function to compute the latency of a spiking 
response.

2025, Alexander Heimel, Robin Haak
"""

import numpy as np
import scipy.io as sio
import os
from latenzy import latenzy

np.random.seed(1)

data_path = os.path.dirname(os.path.dirname(__file__))
mat_data = sio.loadmat(os.path.join(data_path, 'example_data', 'Topo2_20220126_AP.mat'),
                       struct_as_record=False, squeeze_me=True)
sAP = mat_data['sAP']

# Apply inclusion criteria
clusters = sAP.sCluster
is_good = np.array([c.KilosortGood for c in clusters])
low_contam = np.array([c.Contamination for c in clusters]) < 0.1
in_primary_visual = np.array(['primary visual' in c.Area.lower() for c in clusters])
idx_incl = (is_good | low_contam) & in_primary_visual
spike_times_agg = [c.SpikeTimes for i, c in enumerate(clusters) if idx_incl[i]]
event_times = sAP.cellBlock[3].vecStimOnTime
spike_times = spike_times_agg[15]  # 15

# Compute latency for a single cluster
result, s_latenzy = latenzy(
    spike_times,
    event_times,
    use_dur=1,
    resamp_num=100,
    jitter_size=2,
    peak_alpha=0.05,
    do_stitch=True,
    use_par_pool=False,
    use_direct_quant=False,
    restrict_neg=True,
    make_plots=1
)

print(result)
print(s_latenzy)
