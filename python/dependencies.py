# -*- coding: utf-8 -*-
"""
dependencies.py

helper functions for latenzy and latenzy2

2025, Alexander Heimel and Robin Haak, translation from MATLAB version by Robin Haak
"""
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
from scipy.stats import norm
from scipy.interpolate import interp1d
from matplotlib import gridspec


def get_pseudo_times(spike_times, event_times, use_dur=None, discard_edges=False):
    """
    Perform data-stitching.

    Parameters:
        spike_times (np.ndarray): Spike times (s)
        event_times (np.ndarray): Event times (s)
        use_dur (tuple or float): Time window around event (start, end)
        discard_edges (bool): Whether to discard spikes outside event-aligned window

    Returns:
        pseudo_spike_times (np.ndarray): Stitched spike times
        pseudo_event_times (np.ndarray): Stitched event times
    """

    # Ensure correct orientation and sorting
    spike_times = np.sort(np.ravel(spike_times))
    event_times = np.sort(np.ravel(event_times))

    # Check and format use_dur
    if use_dur is None or len(np.atleast_1d(use_dur)) == 0:
        if len(event_times) < 2:
            raise ValueError("Cannot infer use_dur: need at least two event times.")
        use_dur = np.min(np.diff(event_times))
    use_dur = np.atleast_1d(use_dur)
    if np.isscalar(use_dur) or use_dur.size == 1:
        use_dur = np.sort([0, float(use_dur)])
    else:
        use_dur = np.sort(use_dur)

    assert use_dur[1] > use_dur[0], f"The second element of use_dur must be greater than the first: {use_dur}"

    # Initialization
    sample_num = spike_times.size
    event_num = event_times.size
    pseudo_spike_t = []
    pseudo_event_times = np.full(event_num, np.nan)
    duration = use_dur[1] - use_dur[0]
    pseudo_event_t = 0
    last_used_sample = -1
    first_sample = None
    pseudo_t0 = 0

    # Loop over each event
    for this_event in range(event_num):
        event_t = event_times[this_event] + use_dur[0]

        # Find start and end indices of spikes in this window
        start_idx_candidates = np.where(spike_times >= event_t)[0]
        end_idx_candidates = np.where(spike_times < (event_t + duration))[0]

        start_idx = start_idx_candidates[0] if start_idx_candidates.size > 0 else None
        end_idx = end_idx_candidates[-1] if end_idx_candidates.size > 0 else None

        # Handle edge cases for missing spikes
        if start_idx is None:
            start_idx = sample_num
        if end_idx is None or (start_idx is not None and start_idx > end_idx):
            start_idx = None
            end_idx = sample_num - 1

        # Build sample range
        if start_idx is not None:
            use_samples = np.arange(start_idx, end_idx + 1)
        else:
            use_samples = np.array([], dtype=int)

        # Handle edge conditions
        if use_samples.size > 0:
            if this_event == 0 and not discard_edges:
                use_samples = np.arange(0, use_samples[-1] + 1)
            elif this_event == event_num - 1 and not discard_edges:
                use_samples = np.arange(use_samples[0], sample_num)

        # Spike times to add
        add_t = spike_times[use_samples] if use_samples.size > 0 else np.array([])
        samples_overlap = use_samples <= last_used_sample

        # Overlap and pseudo-event time adjustment
        if this_event == 0:
            pseudo_event_t = 0
        elif duration > (event_t - event_times[this_event - 1]):
            use_samples = use_samples[~samples_overlap]
            add_t = spike_times[use_samples]
            pseudo_event_t += event_t - event_times[this_event - 1]
        else:
            pseudo_event_t += duration

        # Local pseudo times
        if use_samples.size == 0:
            local_pseudo_t = np.array([])
        else:
            last_used_sample = use_samples[-1]
            local_pseudo_t = add_t - event_t + pseudo_event_t

        # Record first sample
        if first_sample is None and use_samples.size > 0:
            first_sample = use_samples[0]
            pseudo_t0 = pseudo_event_t

        pseudo_spike_t.append(local_pseudo_t)
        pseudo_event_times[this_event] = pseudo_event_t

    # Add beginning spikes 
    if not discard_edges and first_sample is not None and first_sample > 0:
        step_begin = spike_times[first_sample] - spike_times[first_sample - 1]
        samp_add_beginning = np.arange(0, first_sample)
        if samp_add_beginning.size > 0:
            prepend = spike_times[samp_add_beginning]
            prepend_shifted = prepend - prepend[0] + pseudo_t0 - step_begin - np.ptp(prepend)
            pseudo_spike_t.insert(0, prepend_shifted)

    # Add end spikes
    if not discard_edges:
        indices = np.where(spike_times > (event_times[-1] + duration))[0]
        if indices.size == 0:
            start_idx = np.where(spike_times >= event_times[-1])[0][0]
            samp_add_end = np.arange(start_idx, sample_num)
        else:
            samp_add_end = np.arange(indices[0], sample_num)

        if samp_add_end.size > 0:
            append = spike_times[samp_add_end] - event_times[-1] + pseudo_event_t
            pseudo_spike_t.append(append)

    # Recombine
    valid_arrays = [p for p in pseudo_spike_t if p.size > 0]
    if valid_arrays:
        pseudo_spike_times = np.concatenate(valid_arrays)
    else:
        pseudo_spike_times = np.array([])
        
    pseudo_event_times = pseudo_event_times + abs(use_dur[0])

    return pseudo_spike_times, pseudo_event_times


def get_rel_spike_times(spike_times, event_times, use_dur=None, add_artif_spikes=False):
    """
    Create a vector of spike times relative to event times.

    Parameters:
        spike_times (np.ndarray): Array of spike times (shape: [S])
        event_times (np.ndarray): Array of event times (shape: [T])
        use_dur (float or tuple): Scalar or (pre, post) interval around event (in seconds)
        add_artif_spikes (bool): If True, add artificial spikes at epoch edges

    Returns:
        rel_spike_times (np.ndarray): Sorted vector of spike times relative to events
        spikes_per_event (list of np.ndarray): Spike times per event, relative to that event
    """
    spike_times = np.sort(spike_times).reshape(-1)
    event_times = np.sort(event_times).reshape(-1)

    # Default window duration if not specified
    if use_dur is None:
        if len(event_times) < 2:
            raise ValueError("use_dur must be provided if fewer than 2 events.")
        use_dur = (0, np.min(np.diff(event_times)))

    if np.isscalar(use_dur):
        use_dur = (0, use_dur)

    use_dur = tuple(sorted(use_dur))
    if use_dur[1] <= use_dur[0]:
        raise ValueError(f"Second element of use_dur must be greater than first: got {use_dur}")

    # Compute spike times relative to each event
    spikes_per_event = []
    for t in event_times:
        mask = (spike_times > t + use_dur[0]) & (spike_times < t + use_dur[1])
        rel_spikes = spike_times[mask] - t
        spikes_per_event.append(rel_spikes)

    rel_spike_times = np.sort(np.concatenate(spikes_per_event)) if spikes_per_event else np.array([])
    
    if add_artif_spikes and rel_spike_times.size > 0:
        rel_spike_times = np.unique(np.concatenate(([use_dur[0]], rel_spike_times, [use_dur[1]])))

    return rel_spike_times, spikes_per_event


def calc_temp_diff(spike_times, event_times, use_dur):
    """
    Compute temporal offset vector.

    Parameters:
        spike_times (np.ndarray): Spike timestamps
        event_times (np.ndarray): Event timestamps
        use_dur (tuple or float): Time window to consider around each event

    Returns:
        temp_diff (np.ndarray): Difference between spike fraction and linear fraction
        rel_spike_times (np.ndarray): Spike times relative to event
        spike_frac (np.ndarray): Fractional spike indices
        frac_linear (np.ndarray): Linear fraction (uniformly spaced)
    """
    temp_diff = np.array([])
    spike_frac = np.array([])
    frac_linear = np.array([])

    rel_spike_times,_ = get_rel_spike_times(spike_times, event_times, use_dur, add_artif_spikes=True)
    if rel_spike_times is None or len(rel_spike_times) == 0:
        return temp_diff, rel_spike_times, spike_frac, frac_linear

    rel_spike_times = get_distinct_spikes(rel_spike_times)

    num_spikes = len(rel_spike_times)
    spike_frac = np.linspace(1/num_spikes, 1, num_spikes)

    frac_linear = (rel_spike_times - rel_spike_times[0]) / (rel_spike_times[-1] - rel_spike_times[0])
    temp_diff = spike_frac - frac_linear

    return temp_diff, rel_spike_times, spike_frac, frac_linear



def get_distinct_spikes(spike_times):
    """
    Introduce minimal jitter to repeating spike times to ensure all are unique.

    Parameters:
        spike_times (np.ndarray): Array of spike times

    Returns:
        np.ndarray: Array with minimally jittered spike times (sorted and unique)
    """
    spike_times = np.sort(spike_times)
    unique_offset = np.max(np.finfo(spike_times.dtype).eps * np.abs(spike_times))

    # Identify repeated or nearly repeated spikes
    diffs = np.diff(spike_times)
    idx_repeat = np.concatenate(([False], diffs < unique_offset))

    while np.any(idx_repeat):
        not_unique = spike_times[idx_repeat]
        n = len(not_unique)

        # Random jitter between ±1–10× unique_offset
        jitter = np.concatenate([
            1 + 9 * np.random.rand(n),      # positive: [1, 10]
            -1 - 9 * np.random.rand(n)      # negative: [-10, -1]
        ])
        perm_idx = my_randperm(len(jitter), n)
        jitter = jitter[perm_idx] * unique_offset

        spike_times[idx_repeat] = not_unique + jitter
        spike_times = np.sort(spike_times)

        diffs = np.diff(spike_times)
        idx_repeat = np.concatenate(([False], diffs < unique_offset))

    return spike_times


def my_randperm(n, k=None):
    """
    MATLAB-like randperm function.
    
    Parameters:
        n (int): Total number of elements
        k (int, optional): Number of elements to select. Defaults to n.
    
    Returns:
        np.ndarray: Random permutation of indices
    """
    if k is None:
        k = n
    return np.random.permutation(n)[:k]



def run_jitter_bootstraps(spike_times, event_times, use_dur, resamp_num,
                          jitter_size, use_par_pool=False):
    """
    Run jittered bootstrap analysis.

    Parameters:
        spike_times (np.ndarray): Spike timestamps
        event_times (np.ndarray): Event timestamps
        use_dur (tuple): (start, end) time window around events
        resamp_num (int): Number of resampling iterations
        jitter_size (float): Jitter multiplier
        use_par_pool (bool): Whether to run in parallel

    Returns:
        peaks_rand_d (np.ndarray): Peak deviations from null iterations
        resamp_d (list of np.ndarray): Jittered temporal deviation values
        resamp_t (list of np.ndarray): Corresponding time vectors
    """
    spike_times = np.sort(np.asarray(spike_times).flatten())
    event_times = np.sort(np.asarray(event_times).flatten())
    event_num = len(event_times)
    full_duration = use_dur[1] - use_dur[0]

    # Generate jitter matrix: eventNum x resampNum
    jitter_per_trial = jitter_size * full_duration * (np.random.rand(event_num, resamp_num) - 0.5) * 2

    peaks_rand_d = np.full(resamp_num, np.nan)
    resamp_d = [None] * resamp_num
    resamp_t = [None] * resamp_num

    def single_resample(resamp_idx):
        rand_event_t = event_times + jitter_per_trial[:, resamp_idx]
        rand_d, rand_t, _, _ = calc_temp_diff(spike_times, rand_event_t, use_dur)

        if len(rand_d) == 0:
            return None, None, np.nan

        max_val = np.max(rand_d)
        min_val = np.min(rand_d)
        max_rand_d = min_val if abs(min_val) >= abs(max_val) else max_val

        return rand_d, rand_t, max_rand_d

    if use_par_pool:
        with ThreadPoolExecutor() as executor:
            results = list(executor.map(single_resample, range(resamp_num)))
    else:
        results = [single_resample(i) for i in range(resamp_num)]

    for i, (rd, rt, peak) in enumerate(results):
        resamp_d[i] = rd
        resamp_t[i] = rt
        if not np.isnan(peak):
            peaks_rand_d[i] = peak

    return peaks_rand_d, resamp_d, resamp_t


def compute_pval(max_d, max_rand_d, use_direct_quant=False):
    """
    Compute p-values and z-scores for observed max deviations.

    Parameters:
        max_d (array-like): Observed maximum values
        max_rand_d (array-like): Randomized max values (null distribution)
        use_direct_quant (bool): If True, use empirical quantiles

    Returns:
        p_vals (np.ndarray): p-values
        z_scores (np.ndarray): z-scores
    """
    max_d = np.atleast_1d(max_d)
    max_rand_d = np.unique(np.sort(np.asarray(max_rand_d)))
    if max_rand_d.size == 0:
        raise ValueError("max_rand_d must not be empty.")

    p_vals = np.full(max_d.shape, np.nan)

    if use_direct_quant:
        N = len(max_rand_d)
        for i, val in enumerate(max_d):
            if np.isnan(val) or val < max_rand_d[0]:
                p_vals[i] = 1
            elif np.isinf(val) or val > max_rand_d[-1]:
                p_vals[i] = 1 / (1 + N)
            else:
                # Empirical rank-based interpolation
                interp = interp1d(max_rand_d, np.arange(1, N + 1), kind='linear', fill_value='extrapolate')
                value_rank = interp(val)
                p_vals[i] = 1 - (value_rank / (1 + N))

        z_scores = -norm.ppf(p_vals / 2)
    else:
        rand_mu = np.mean(max_rand_d)
        rand_var = np.var(max_rand_d)
        p_vals, z_scores, _, _ = get_gumbel(rand_mu, rand_var, max_d)

    return p_vals, z_scores


def get_gumbel(mean_val, var_val, x_vals):
    """
    Compute p-values and z-scores for values under a fitted Gumbel distribution.

    Parameters:
        mean_val (float): Mean of the maximum value distribution
        var_val (float): Variance of the maximum value distribution
        x_vals (array-like): Observed maximum value(s)

    Returns:
        p_vals (np.ndarray): p-values under the Gumbel distribution
        z_scores (np.ndarray): Corresponding z-scores
        mode (float): Gumbel mode parameter
        beta (float): Gumbel scale parameter
    """
    x_vals = np.asarray(x_vals)
    euler_mascheroni = 0.5772156649015329  # Euler–Mascheroni constant

    # Estimate Gumbel parameters from mean and variance
    beta = np.sqrt(6 * var_val) / np.pi
    mode = mean_val - beta * euler_mascheroni

    # Gumbel CDF
    gumbel_cdf = np.exp(-np.exp(-(x_vals - mode) / beta))
    p_vals = 1 - gumbel_cdf

    # Initial z-score conversion
    with np.errstate(divide='ignore'):
        z_scores = -norm.ppf(p_vals / 2)

    # Handle infinite z-scores by re-approximating p
    inf_mask = np.isinf(z_scores)
    if np.any(inf_mask):
        p_vals[inf_mask] = np.exp((mode - x_vals[inf_mask]) / beta)
        z_scores[inf_mask] = -norm.ppf(p_vals[inf_mask] / 2)

    return p_vals, z_scores, mode, beta


def make_latenzy_figs(s_latenzy, spike_times, event_times, use_dur, make_plots):
    """
    Generate figures to visualize latenzy analysis.
    
    Parameters:
        s_latenzy: dict with keys like 'latency', 'peakTimes', etc.
        spike_times: 1D np.ndarray
        event_times: 1D np.ndarray
        use_dur: tuple or list (start, end)
        make_plots: int (1 to make full raster, 0 to skip)

    Returns:
        fig_handles: list of matplotlib Axes handles
    """
    latency = s_latenzy['latency']
    peak_times = s_latenzy['peakTimes']
    peak_vals = s_latenzy['peakVals']
    real_frac = s_latenzy['realFrac']
    frac_lin = s_latenzy['fracLin']
    real_diff = s_latenzy['realDiff']
    real_time = s_latenzy['realTime']
    mean_real_diff = s_latenzy['meanRealDiff']
    rand_diff = s_latenzy['randDiff']
    rand_time = s_latenzy['randTime']
    mean_rand_diff = s_latenzy['meanRandDiff']
    #peak_z = s_latenzy['peakZ']
    pvals_peak = s_latenzy['pValsPeak']
    
    latenzy_idx = s_latenzy['latenzyIdx']

    num_iters = len(peak_times)
    use_colors = plt.cm.Blues(np.linspace(1, 0.5, num_iters))
    line_width = 1.5
    marker_size = 60

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs = axs.flatten()
    fig_handles = axs.tolist()

    if make_plots == 1:
        # Raster plot
        axs[0].cla()
        for i, evt in enumerate(event_times):
            mask = (spike_times >= evt + use_dur[0]) & (spike_times <= evt + use_dur[1])
            rel_times = spike_times[mask] - evt
            axs[0].vlines(rel_times, i + 0.5, i + 1.5, color='black', linewidth=0.5)
        axs[0].axvline(latency, color=(0.8627, 0.0784, 0.2353), linestyle='--', linewidth=line_width)
        axs[0].set_ylim(-0.5, len(event_times) + 0.5)
        axs[0].set_ylabel('Trial')
        axs[0].set_xlim(use_dur)
        axs[0].set_xlabel('Time from event (s)')
        axs[0].set_title('Aligned spikes')
        axs[0].tick_params(direction='out')

    # Cumulative spike fractions
    axs[1].cla()
    legend_lines = []
    legend_labels = []
    for i in range(num_iters):
        axs[1].plot(real_time[i], frac_lin[i], color='gray', linewidth=line_width)
        line = axs[1].plot(real_time[i], real_frac[i], color=use_colors[i], linewidth=line_width)[0]
        legend_lines.append(line)
        legend_labels.append(f'{i + 1}')
    axs[1].set_xlim(use_dur)
    axs[1].set_ylim([0,1])
    axs[1].set_xlabel('Time from event (s)')
    axs[1].set_ylabel('Fractional spike position')
    axs[1].set_title('Cumulative spikes')
    axs[1].legend(legend_lines, legend_labels, title='Iteration', loc='lower right')
    axs[1].tick_params(direction='out')

    # Deviation from linear
    axs[2].cla()
    axs[2].axhline(0, color='gray', linestyle='--', linewidth=line_width)
    for i in range(num_iters):
        axs[2].plot(real_time[i], real_diff[i], color=use_colors[i], linewidth=line_width)
    axs[2].scatter(np.array(peak_times)[~np.array(latenzy_idx)],
               np.array(peak_vals)[~np.array(latenzy_idx)],
               s=marker_size, color='black', marker='x', linewidth=line_width, zorder=3)
    axs[2].scatter(np.array(peak_times)[latenzy_idx],
                   np.array(peak_vals)[latenzy_idx],
                   s=marker_size, color=(0.8627, 0.0784, 0.2353), marker='x',
                   linewidth=line_width, zorder=3)
    axs[2].set_xlim(use_dur)
    axs[2].set_xlabel('Time from event (s)')
    axs[2].set_ylabel('Deviation (Δfraction)')
    axs[2].set_title('Deviation from uniform')
    axs[2].tick_params(direction='out')


    # Real + jittered deviation
    axs[3].cla()
    lat_idx = np.where(latenzy_idx)[0]
    for i in range(len(rand_diff[0])):
        axs[3].plot(rand_time[lat_idx[0]][i],
                    rand_diff[lat_idx[0]][i] - mean_rand_diff[lat_idx[0]][i],
                    color='gray', linewidth=line_width)
    axs[3].plot(real_time[lat_idx[0]],
                real_diff[lat_idx[0]] - mean_real_diff[lat_idx[0]],
                color=use_colors[lat_idx[0]], linewidth=line_width)
    axs[3].scatter(peak_times[lat_idx[0]],
                   peak_vals[lat_idx[0]] - mean_real_diff[lat_idx[0]],
                   s=marker_size, color=(0.8627, 0.0784, 0.2353), marker='x',
                   linewidth=line_width, zorder=3)
    if latenzy_idx[0]:
        axs[3].set_xlim(use_dur)
    else:
        idx = np.where(latenzy_idx)[0][0]
        axs[3].set_xlim([use_dur[0], peak_times[idx - 1]])
    axs[3].set_xlabel('Time from event (s)')
    axs[3].set_ylabel('Deviation (Δfraction)')
    axs[3].set_title(f'Real + jittered data (p={pvals_peak[lat_idx[0]][0]:.4f})')
    axs[3].tick_params(direction='out')

    fig.suptitle(f'latenZy estimate = {latency:.4f}s', fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    return fig_handles



def calc_temp_diff2(spikes_per_trial1, spikes_per_trial2, use_dur, use_fast_interp=False):
    """
    Compute temporal offset between spike trains from two conditions.

    Parameters:
        spikes_per_trial1: list of 1D arrays (condition 1)
        spikes_per_trial2: list of 1D arrays (condition 2)
        use_dur: tuple/list (start, end)
        use_fast_interp: bool, skip (currently same interp used)

    Returns:
        temp_diff: difference of spike fractions after linear subtraction
        rel_spike_times_agg: combined time vector
        spike_frac1: interpolated spike fraction for cond 1
        rel_spike_times1: raw cond 1 spike times
        spike_frac2: interpolated spike fraction for cond 2
        rel_spike_times2: raw cond 2 spike times
        temp_diff_unsub: uncorrected difference
        frac_linear: linear fit baseline
    """
    # Flatten and deduplicate spikes
    rel_spike_times1 = get_distinct_spikes(np.concatenate(spikes_per_trial1))
    rel_spike_times2 = get_distinct_spikes(np.concatenate(spikes_per_trial2))

    # Combine and sort, pad with artificial start/end spikes
    rel_spike_times_agg = np.sort(np.concatenate([rel_spike_times1, rel_spike_times2]))
    if rel_spike_times_agg.size == 0 or np.max(rel_spike_times_agg) < use_dur[1]:
        rel_spike_times_agg = np.concatenate((
            [use_dur[0]],
            rel_spike_times_agg,
            [use_dur[1]]
        ))
    else:
        rel_spike_times_agg = np.concatenate(([use_dur[0]], rel_spike_times_agg))

    # Spike fractions
    num_sp1 = len(rel_spike_times1)
    num_sp2 = len(rel_spike_times2)
    num_ev1 = len(spikes_per_trial1)
    num_ev2 = len(spikes_per_trial2)

    unique_spike_fracs1 = np.arange(1, num_sp1 + 1) / num_ev1
    interp_x1 = np.concatenate(([use_dur[0]], rel_spike_times1, [use_dur[1]]))
    interp_y1 = np.concatenate(([0], unique_spike_fracs1, [num_sp1 / num_ev1]))
    spike_frac1 = np.interp(rel_spike_times_agg, interp_x1, interp_y1, left=np.nan, right=np.nan)
    spike_frac1 = fillnans(spike_frac1, num_sp1, num_ev1)

    unique_spike_fracs2 = np.arange(1, num_sp2 + 1) / num_ev2
    interp_x2 = np.concatenate(([use_dur[0]], rel_spike_times2, [use_dur[1]]))
    interp_y2 = np.concatenate(([0], unique_spike_fracs2, [num_sp2 / num_ev2]))
    spike_frac2 = np.interp(rel_spike_times_agg, interp_x2, interp_y2, left=np.nan, right=np.nan)
    spike_frac2 = fillnans(spike_frac2, num_sp2, num_ev2)

    # Difference
    temp_diff = spike_frac1 - spike_frac2
    temp_diff_unsub = temp_diff.copy()

    # Linear baseline
    frac_linear = temp_diff[0] + (temp_diff[-1] - temp_diff[0]) * \
                  (rel_spike_times_agg - rel_spike_times_agg[0]) / \
                  (rel_spike_times_agg[-1] - rel_spike_times_agg[0])
    temp_diff = temp_diff - frac_linear

    return (temp_diff, rel_spike_times_agg, spike_frac1, rel_spike_times1,
            spike_frac2, rel_spike_times2, temp_diff_unsub, frac_linear)


def fillnans(vec, int_sp, int_t):
    """
    Fill leading and trailing NaNs in spike fraction vector.

    Parameters:
        vec (np.ndarray): The array to modify (1D)
        int_sp (int): Number of spikes
        int_t (int): Number of trials

    Returns:
        np.ndarray: Modified vector with edge NaNs filled
    """
    vec = vec.copy()
    isnan = np.isnan(vec)

    if isnan[0] or isnan[-1]:
        nan_diff = np.diff(isnan.astype(int))

    if isnan[0]:
        leading_nan_len = np.argmax(nan_diff == 1) + 1 if 1 in nan_diff else len(vec)
        vec[:leading_nan_len] = 1 / int_t

    if isnan[-1]:
        reversed_diff = np.diff(isnan[::-1].astype(int))
        lagging_nan_len = np.argmax(reversed_diff == 1) + 1 if 1 in reversed_diff else len(vec)
        vec[-lagging_nan_len:] = int_sp / int_t

    return vec


def run_swap_bootstraps(spikes_per_event1, spikes_per_event2, use_dur,
                        resamp_num, use_par_pool=False, use_fast_interp=False):
    """
    Run trial-swapping bootstraps to build null distribution of max deviations.

    Parameters:
        spikes_per_event1: list of np.ndarray, trials for condition 1
        spikes_per_event2: list of np.ndarray, trials for condition 2
        use_dur: (start, end) time window to use
        resamp_num: number of bootstrap iterations
        use_par_pool: run in parallel
        use_fast_interp: if True, enable faster interpolation (optional)

    Returns:
        peaks_rand_d: np.ndarray of max deviations from null
        resamp_d: list of temporal deviation arrays
        resamp_t: list of corresponding time arrays
    """
    num_ev1 = len(spikes_per_event1)
    num_ev2 = len(spikes_per_event2)
    spikes_agg = spikes_per_event1 + spikes_per_event2
    num_ev_tot = num_ev1 + num_ev2

    idx_spikes_empty = np.array([len(sp) == 0 for sp in spikes_agg])
    peaks_rand_d = np.full(resamp_num, np.nan)
    resamp_d = [None] * resamp_num
    resamp_t = [None] * resamp_num

    def single_bootstrap(resamp_idx):
        shuffled_idx = my_randperm(num_ev_tot)
        use_rand1 = shuffled_idx[:num_ev1]
        use_rand2 = shuffled_idx[num_ev1:]

        spikes_rand1 = [spikes_agg[i] for i in use_rand1]
        spikes_rand2 = [spikes_agg[i] for i in use_rand2]

        if all(idx_spikes_empty[i] for i in use_rand1) and all(idx_spikes_empty[i] for i in use_rand2):
            return None, None, np.nan

        rand_d, rand_t, *_ = calc_temp_diff2(spikes_rand1, spikes_rand2, use_dur, use_fast_interp)

        max_val = np.max(rand_d)
        min_val = np.min(rand_d)
        max_rand_d = min_val if abs(min_val) >= abs(max_val) else max_val

        return rand_d, rand_t, max_rand_d

    if use_par_pool:
        with ThreadPoolExecutor() as executor:
            results = list(executor.map(single_bootstrap, range(resamp_num)))
    else:
        results = [single_bootstrap(i) for i in range(resamp_num)]

    for i, (rd, rt, peak) in enumerate(results):
        resamp_d[i] = rd
        resamp_t[i] = rt
        if not np.isnan(peak):
            peaks_rand_d[i] = peak

    return peaks_rand_d, resamp_d, resamp_t


def make_latenzy2_figs(s_latenzy2, spike_times1, event_times1, spike_times2, event_times2, use_dur, make_plots):
    """
    Generate figures to visualize latency analysis for two conditions.
    
    Parameters:
        s_latenzy2: dict with latency analysis results (keys similar to make_latenzy)
        spike_times1: 1D np.ndarray for condition 1 spike times
        event_times1: 1D np.ndarray for condition 1 event times
        spike_times2: 1D np.ndarray for condition 2 spike times
        event_times2: 1D np.ndarray for condition 2 event times
        use_dur: tuple/list (start, end) time window for plotting
        make_plots: int (1 or 9 to make raster plots, else no raster)
    
    Returns:
        fig_handles: list of matplotlib Axes handles
    """
    latency = s_latenzy2['latency']
    peak_times = s_latenzy2['peakTimes']
    peak_vals = s_latenzy2['peakVals']
    real_frac = s_latenzy2['realFrac']
    temp_diff_unsub = s_latenzy2['diffUnSub']
    frac_lin = s_latenzy2['fracLin']
    real_diff = s_latenzy2['realDiff']
    real_time = s_latenzy2['realTime']
    mean_real_diff = s_latenzy2['meanRealDiff']
    rand_diff = s_latenzy2['randDiff']
    rand_time = s_latenzy2['randTime']
    mean_rand_diff = s_latenzy2['meanRandDiff']
    #peak_z = s_latenzy2['peakZ']
    pvals_peak = s_latenzy2['pValsPeak']
    latenzy_idx = s_latenzy2['latenzyIdx']

    num_iters = len(peak_times)
    use_colors = plt.cm.Blues(np.linspace(1, 0.5, num_iters))
    line_width = 1.5
    marker_size = 60

    fig = plt.figure(figsize=(18, 10))
    gs = gridspec.GridSpec(2, 3, figure=fig)
    axs = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]  # 6 axes total

    # Raster plots if requested
    if make_plots == 1 or make_plots == 9:
        axs[0].cla()
        for i, evt in enumerate(event_times1):
            times = spike_times1[(spike_times1 > evt + use_dur[0]) & (spike_times1 < evt + use_dur[1])] - evt
            axs[0].vlines(times, i + 0.5, i + 1.5, color='k')
        axs[0].axvline(x=latency, linestyle='--', linewidth=line_width, color=(0.8627, 0.0784, 0.2353))
        axs[0].set_ylim(-0.5, len(event_times1) + 0.5)
        axs[0].set_xlim(use_dur)
        axs[0].set_title('Aligned spikes condition 1')
        axs[0].set_xlabel('Time from event (s)')
        axs[0].set_ylabel('Trial')

        axs[3].cla()
        for i, evt in enumerate(event_times2):
            times = spike_times2[(spike_times2 > evt + use_dur[0]) & (spike_times2 < evt + use_dur[1])] - evt
            axs[3].vlines(times, i + 0.5, i + 1.5, color='k')
        axs[3].axvline(x=latency, linestyle='--', linewidth=line_width, color=(0.8627, 0.0784, 0.2353))
        axs[3].set_ylim(-0.5, len(event_times2) + 0.5)
        axs[3].set_xlim(use_dur)
        axs[3].set_title('Aligned spikes condition 2')
        axs[3].set_xlabel('Time from event (s)')
        axs[3].set_ylabel('Trial')

    # Cumulative spike count
    axs[1].cla()
    axs[1].plot(real_time[0], real_frac[0][0], label='1', color=(0.8510, 0.4510, 0.1340), linewidth=line_width)
    axs[1].plot(real_time[0], real_frac[0][1], label='2', color=(0.9255, 0.7255, 0.5670), linewidth=line_width)
    axs[1].set_xlim(use_dur)
    axs[1].set_title('Cumulative spikes')
    axs[1].set_xlabel('Time from event (s)')
    axs[1].set_ylabel('Cumulative spike count (norm.)')
    axs[1].legend(title='Condition', loc='lower right')

    # Linear difference
    axs[2].cla()
    axs[2].plot(real_time[0], frac_lin[0], color='gray', linewidth=line_width)
    axs[2].plot(real_time[0], temp_diff_unsub[0], color=use_colors[0], linewidth=line_width)
    axs[2].set_xlim(use_dur)
    axs[2].set_title('Condition 1 - condition 2')
    axs[2].set_xlabel('Time from event (s)')
    axs[2].set_ylabel('Spike count difference')

    # Offset from linear
    axs[4].cla()
    axs[4].axhline(y=0, color='gray', linestyle='--', linewidth=line_width)
    
    line_handles = []
    for iter in range(num_iters):
        line = axs[4].plot(real_time[iter], real_diff[iter],
                           color=use_colors[iter], linewidth=line_width)[0]
        line_handles.append(line)
    axs[4].scatter(np.array(peak_times)[~np.array(latenzy_idx)],
                   np.array(peak_vals)[~np.array(latenzy_idx)],
                   marker='x', s=marker_size, color='k', linewidths=line_width, zorder=3)
    axs[4].scatter(np.array(peak_times)[latenzy_idx],
                   np.array(peak_vals)[latenzy_idx],
                   marker='x', s=marker_size, color=(0.8627, 0.0784, 0.2353),
                   linewidths=line_width, zorder=3)
    axs[4].set_xlim(use_dur)
    axs[4].set_title('Offset from linear')
    axs[4].set_xlabel('Time from event (s)')
    axs[4].set_ylabel('Deviation (Δcount)')
    axs[4].legend(line_handles, [str(i+1) for i in range(num_iters)],
                  title='Iteration', loc='lower right')

    # Real + shuffled deviation, mean-subtracted
    axs[5].cla()
    lat_idx = np.where(latenzy_idx)[0]
    for i in range(len(rand_diff[lat_idx[0]])):
        axs[5].plot(rand_time[lat_idx[0]][i],
                    rand_diff[lat_idx[0]][i] - mean_rand_diff[lat_idx[0]][i],
                    color='gray', linewidth=line_width)
    axs[5].plot(real_time[lat_idx[0]],
                real_diff[lat_idx[0]] - mean_real_diff[lat_idx[0]],
                color=use_colors[lat_idx[0]], linewidth=line_width)
    axs[5].scatter(peak_times[lat_idx[0]],
                   peak_vals[lat_idx[0]] - mean_real_diff[lat_idx[0]],
                   s=marker_size, color=(0.8627, 0.0784, 0.2353), marker='x',
                   linewidth=line_width, zorder=3)

    if latenzy_idx[0]:
        axs[5].set_xlim(use_dur)
    else:
        idx = lat_idx[0]
        axs[5].set_xlim([use_dur[0], peak_times[idx - 1]])
    
    axs[5].set_xlabel('Time from event (s)')
    axs[5].set_ylabel('Deviation (Δcount)')
    axs[5].set_title(f'Real + shuffled data (p={pvals_peak[lat_idx[0]][0]:.4f})')
    axs[5].tick_params(direction='out')

    fig.suptitle(f'latenZy2 estimate = {latency:.4f}s', fontweight='bold')
    plt.tight_layout()
    plt.show()

    return axs


