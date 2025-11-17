# -*- coding: utf-8 -*-
"""
latenzy.py

Contains latenzy and latenzy2 to compute latencies for spiking responses
See Haak et al. 2025

2025, Alexander Heimel, translated from MATLAB version by Robin Haak
"""

import numpy as np

from dependencies import (get_pseudo_times,calc_temp_diff,
                          run_jitter_bootstraps,compute_pval,make_latenzy_figs,
                          get_rel_spike_times,calc_temp_diff2,
                          run_swap_bootstraps,make_latenzy2_figs)


def latenzy(spike_times, 
            event_times, 
            use_dur=None, 
            resamp_num=100, 
            jitter_size=2,
            peak_alpha=0.05, 
            do_stitch=True, 
            use_par_pool=False,
            use_direct_quant=False, 
            restrict_neg=False, 
            make_plots=0):
    """
    Compute event-related spiking latency.

    Parameters
    ----------
    spike_times : array-like
        Spike times (seconds).
    event_times : array-like
        Event times (seconds).
    use_dur : scalar or list of two floats, optional
        Time window to include after/around event times (seconds). If scalar, interpreted as duration after each event, with start time set to zero. Default: [0, min(diff(event_times))].
    resamp_num : int, optional
        Number of resamples for bootstrapping (default: 100).
    jitter_size : float, optional
        Temporal jitter window relative to use_dur (seconds, default: 2).
    peak_alpha : float, optional
        Significance threshold for peak detection (default: 0.05).
    do_stitch : bool, optional
        Perform data stitching, highly recommended! (default: True).
    use_par_pool : bool, optional
        Use parallel pool for resamples (default: False).
    use_direct_quant : bool, optional
        Use empirical null-distribution rather than Gumbel approximation (default: False).
    restrict_neg : bool, optional
        Restrict negative latencies (default: False).
    make_plots : int, optional
        Plotting switch (0=none, 1=raster+traces, 2=traces only, default: 0).

    Returns
    -------
    latency : float
        Response latency (seconds), NaN if no latency could be estimated.
    s_latenzy : dict
        Dictionary with fields:
            - latency: response latency (s)
            - peakTimes: detected peak times, one per iteration (s)
            - peakVals: detected peak values, one per iteration
            - realFrac: see plotting function for details
            - fracLin: idem
            - realDiff: idem
            - realTime: idem
            - meanRealDiff: idem
            - randDiff: idem
            - randTime: idem
            - meanRandDiff: idem
            - pVals: p-values for the observed peak maxima
            - peakZ: two-tailed z-scores corresponding to the p-values
            - latenzyIdx: use to index arrays above
            - figHandles: figure handles (if make_plots > 0)
    """

    spike_times = np.asarray(spike_times).flatten()
    event_times = np.asarray(event_times).flatten()

    # Set defailt use_dur
    event_times.sort()
    if use_dur is None:
        use_dur = [0, np.min(np.diff(event_times))]
    elif np.isscalar(use_dur):
        use_dur = sorted([0, use_dur])
    elif len(use_dur) != 2:
        raise ValueError("useMaxDur must be a scalar or a two-element list")

    # Validate use_dur
    if use_dur[1] > 0:
        assert use_dur[0] <= 0, f"When useMaxDur[1] > 0, useMaxDur[0] must be ≤ 0, got {use_dur}"
    elif use_dur[1] == 0:
        assert use_dur[0] < 0, f"When useMaxDur[1] = 0, useMaxDur[0] must be < 0, got {use_dur}"
    else:
        raise ValueError("useMaxDur[1] cannot be negative when useMaxDur[0] is negative!")

    # Main loop
    latency = np.nan
    peak_times_agg = []
    peak_vals_agg = []
    real_frac_agg = []
    frac_lin_agg = []
    real_diff_agg = []
    real_time_agg = []
    mean_real_diff_agg = []
    rand_diff_agg = []
    rand_time_agg = []
    mean_rand_diff_agg = []
    p_val_peak_agg = []
    peak_z_agg = []
    keep_peaks = []
    this_max_dur = use_dur.copy()
    do_continue = True
    this_iter = 0
    give_late_warn = False

    min_latency = 0 if restrict_neg else use_dur[0]

    while do_continue:
        this_iter += 1

        if do_stitch:
            discard_edges = True
            pseudo_spike_times, pseudo_event_times = get_pseudo_times(spike_times, event_times, this_max_dur, discard_edges)
        else:
            pseudo_spike_times = spike_times
            pseudo_event_times = event_times

        real_diff, real_time, spike_frac, frac_linear = calc_temp_diff(pseudo_spike_times, pseudo_event_times, this_max_dur)
        if len(real_diff) < 3:
            return np.nan, {}

        # Peak detection
        max_diff = np.max(real_diff)
        min_diff = np.min(real_diff)
        if abs(min_diff) >= abs(max_diff):
            real_max_d = min_diff
            real_max_idx = np.argmin(real_diff)
        else:
            real_max_d = max_diff
            real_max_idx = np.argmax(real_diff)

        real_peak_t = real_time[real_max_idx]
        real_peak_sub = real_max_d - np.mean(real_diff)

        # Bootstrapping
        peaks_rand, rand_diff, rand_time = run_jitter_bootstraps(
            pseudo_spike_times, pseudo_event_times, this_max_dur, resamp_num,
            jitter_size, use_par_pool)

        mean_rand_diff = np.array([np.mean(rd) if rd is not None and len(rd) > 0 else np.nan for rd in rand_diff])

        peaks_rand_sub = peaks_rand - mean_rand_diff

        # Compute significance
        p_val_peak, peak_z = compute_pval(np.abs(real_peak_sub), peaks_rand_sub[~np.isnan(peaks_rand_sub)], use_direct_quant)

        # Store iteration
        if not np.isnan(real_peak_t):
            peak_vals_agg.append(real_max_d)
            peak_times_agg.append(real_peak_t)
            real_frac_agg.append(spike_frac)
            frac_lin_agg.append(frac_linear)
            real_diff_agg.append(real_diff)
            real_time_agg.append(real_time)
            mean_real_diff_agg.append(np.mean(real_diff))
            rand_diff_agg.append(rand_diff)
            rand_time_agg.append(rand_time)
            mean_rand_diff_agg.append(mean_rand_diff)
            p_val_peak_agg.append(p_val_peak)
            peak_z_agg.append(abs(peak_z))

        if real_peak_t > min_latency and p_val_peak < peak_alpha:
            keep_peaks.append(True)
            this_max_dur[1] = real_peak_t
        else:
            do_continue = False
            keep_peaks.append(False)

    # Aggregate results
    keep_peaks = np.array(keep_peaks, dtype=bool)
    peak_times_arr = np.array(peak_times_agg)

    if np.any(keep_peaks):
        latency = peak_times_arr[keep_peaks][-1]
        if latency > (use_dur[0] + sum(np.abs(use_dur)) / 2) and give_late_warn:
            print("Warning: Estimated latency is late in the window (>T/2). Consider plotting for visual check and/or adjusting window.")
    else:
        return np.nan, {}

    s_latenzy = {
        'latency': latency,
        'peakTimes': peak_times_agg,
        'peakVals': peak_vals_agg,
        'realFrac': real_frac_agg,
        'fracLin': frac_lin_agg,
        'realDiff': real_diff_agg,
        'realTime': real_time_agg,
        'meanRealDiff': mean_real_diff_agg,
        'randDiff': rand_diff_agg,
        'randTime': rand_time_agg,
        'meanRandDiff': mean_rand_diff_agg,
        'pValsPeak': p_val_peak_agg,
        'peakZ': peak_z_agg,
        'latenzyIdx': peak_times_arr == latency
    }

    if make_plots > 0:
        s_latenzy['figHandles'] = make_latenzy_figs(s_latenzy, spike_times, event_times, use_dur, make_plots)

    return latency, s_latenzy

def latenzy2(
    spike_times1,
    event_times1,
    spike_times2=None,
    event_times2=None,
    use_dur=None,
    resamp_num=250,
    peak_alpha=0.05,
    use_par_pool=False,
    use_direct_quant=False,
    restrict_neg=False,
    make_plots=0
):
    """
    Compute latency of spiking difference between two conditions.

    Parameters
    ----------
    spike_times1 : array-like or list of arrays
        Spike times for condition 1 (seconds), or list of aligned spikes per event.
    event_times1 : array-like or None
        Event times for condition 1 (seconds), or None if using aligned input.
    spike_times2 : array-like or list of arrays, optional
        Spike times for condition 2 (seconds), or list of aligned spikes per event. Default: same as spike_times1.
    event_times2 : array-like or None, optional
        Event times for condition 2 (seconds), or None if using aligned input.
    use_dur : scalar or list of two floats, optional
        Time window to include after/around event times (seconds). If scalar, interpreted as duration after each event, with start time set to zero. Default: [0, min(diff(event_times1)), min(diff(event_times2))].
    resamp_num : int, optional
        Number of resamples for bootstrapping (default: 250).
    peak_alpha : float, optional
        Significance threshold for peak detection (default: 0.05).
    use_par_pool : bool, optional
        Use parallel pool for resamples (default: False).
    use_direct_quant : bool, optional
        Use empirical null-distribution rather than Gumbel approximation (default: False).
    restrict_neg : bool, optional
        Restrict negative latencies (default: False).
    make_plots : int, optional
        Plotting switch (0=none, 1=raster+traces, 2=traces only, default: 0).

    Returns
    -------
    latency : float
        Response latency (seconds), NaN if no latency could be estimated.
    s_latenzy2 : dict
        Dictionary with fields:
            - latency: response latency (s)
            - peakTimes: detected peak times, one per iteration (s)
            - peakVals: detected peak values, one per iteration
            - realFrac: see plotting function for details
            - diffUnSub: idem
            - fracLin: idem
            - realDiff: idem
            - realTime: idem
            - meanRealDiff: idem
            - randDiff: idem
            - randTime: idem
            - meanRandDiff: idem
            - pVals: p-values for the observed peak maxima
            - peakZ: two-tailed z-scores corresponding to the p-values
            - latenzyIdx: use to index arrays above
            - figHandles: figure handles (if make_plots > 0)
    """

    # ALT INPUT CASE: aligned spikes in cells
    alt_input = (
        isinstance(spike_times1, list)
        and isinstance(spike_times2, list)
        and event_times1 is None
        and event_times2 is None
    )

    if alt_input:
        spike_times1 = [np.asarray(s) for s in spike_times1]
        spike_times2 = [np.asarray(s) for s in spike_times2]
    else:
        spike_times1 = np.asarray(spike_times1).flatten()
        event_times1 = np.asarray(event_times1).flatten()
        spike_times2 = np.asarray(spike_times2).flatten() if spike_times2 is not None else spike_times1
        event_times2 = np.asarray(event_times2).flatten()

    # Set default use_dur
    if use_dur is None:
        if alt_input:
            use_dur = [0., min(np.max(np.concatenate(spike_times1)), np.max(np.concatenate(spike_times2)))]
        else:
            use_dur = [0., min(np.min(np.diff(event_times1)), np.min(np.diff(event_times2)))]

    if np.isscalar(use_dur):
        use_dur = sorted([0., use_dur])
    elif len(use_dur) != 2:
        raise ValueError("use_dur must be a scalar or two-element list")

    if use_dur[1] <= 0 or (use_dur[1] == 0 and use_dur[0] >= 0):
        raise ValueError("Invalid use_dur: must include a positive duration")

    # Main loop
    latency = np.nan
    peak_times_agg = []
    peak_vals_agg = []
    real_frac_agg = []
    temp_diff_unsub_agg = []
    frac_lin_agg = []
    real_diff_agg = []
    real_time_agg = []
    mean_real_diff_agg = []
    rand_diff_agg = []
    rand_time_agg = []
    mean_rand_diff_agg = []
    p_val_peak_agg = []
    peak_z_agg = []
    keep_peaks = []
    this_max_dur = use_dur.copy()
    do_continue = True
    this_iter = 0

    min_latency = 0. if restrict_neg else use_dur[0]

    while do_continue:
        this_iter += 1

        if alt_input:
            spikes1 = [s[(s > this_max_dur[0]) & (s < this_max_dur[1])] for s in spike_times1]
            spikes2 = [s[(s > this_max_dur[0]) & (s < this_max_dur[1])] for s in spike_times2]
        else:
            _, spikes1 = get_rel_spike_times(spike_times1, event_times1, this_max_dur)
            _, spikes2 = get_rel_spike_times(spike_times2, event_times2, this_max_dur)

        real_diff, real_time, spike_frac1, _, spike_frac2, _, temp_diff_unsub, frac_linear = \
            calc_temp_diff2(spikes1, spikes2, this_max_dur)

        if len(real_diff) < 3:
            return np.nan, {}

        # Peak detection
        max_val = np.max(real_diff)
        min_val = np.min(real_diff)
        real_max_d = min_val if abs(min_val) >= abs(max_val) else max_val
        real_peak_idx = np.argmin(real_diff) if abs(min_val) >= abs(max_val) else np.argmax(real_diff)
        real_peak_t = real_time[real_peak_idx]
        real_peak_sub = real_max_d - np.mean(real_diff)

        # Bootstrapping
        peaks_rand, rand_diff, rand_time = run_swap_bootstraps(
            spikes1, spikes2, this_max_dur, resamp_num, use_par_pool
        )
        mean_rand_diff = np.array([np.mean(rd) if rd is not None and len(rd) > 0 else np.nan for rd in rand_diff])
        peaks_rand_sub = peaks_rand - mean_rand_diff

        # Compute significance
        p_val_peak, peak_z = compute_pval(abs(real_peak_sub), peaks_rand_sub[~np.isnan(peaks_rand_sub)], use_direct_quant)

        # Store iteration
        if not np.isnan(real_peak_t):
            peak_vals_agg.append(real_max_d)
            peak_times_agg.append(real_peak_t)
            real_frac_agg.append([spike_frac1, spike_frac2])
            temp_diff_unsub_agg.append(temp_diff_unsub)
            frac_lin_agg.append(frac_linear)
            real_diff_agg.append(real_diff)
            real_time_agg.append(real_time)
            mean_real_diff_agg.append(np.mean(real_diff))
            rand_diff_agg.append(rand_diff)
            rand_time_agg.append(rand_time)
            mean_rand_diff_agg.append(mean_rand_diff)
            p_val_peak_agg.append(p_val_peak)
            peak_z_agg.append(abs(peak_z))

        if real_peak_t > min_latency and p_val_peak < peak_alpha:
            keep_peaks.append(True)
            this_max_dur[1] = real_peak_t
        else:
            do_continue = False
            keep_peaks.append(False)

    # Aggregate results
    keep_peaks = np.array(keep_peaks)
    these_peak_times = np.array(peak_times_agg)[keep_peaks]

    if len(these_peak_times) == 0:
        return np.nan, {}

    latency = these_peak_times[-1]

    s_latenzy2 = {
        'latency': latency,
        'peakTimes': peak_times_agg,
        'peakVals': peak_vals_agg,
        'realFrac': real_frac_agg,
        'diffUnSub': temp_diff_unsub_agg,
        'fracLin': frac_lin_agg,
        'realDiff': real_diff_agg,
        'realTime': real_time_agg,
        'meanRealDiff': mean_real_diff_agg,
        'randDiff': rand_diff_agg,
        'randTime': rand_time_agg,
        'meanRandDiff': mean_rand_diff_agg,
        'pValsPeak': p_val_peak_agg,
        'peakZ': peak_z_agg,
        'latenzyIdx': np.array(peak_times_agg) == latency
    }

    if make_plots > 0:
        s_latenzy2['figHandles'] = make_latenzy2_figs(
            s_latenzy2, spike_times1, event_times1, spike_times2, event_times2, use_dur, make_plots
        )

    return latency, s_latenzy2

##

