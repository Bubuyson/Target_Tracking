function fa_times = calculateFATimes(start_time, end_time, nfa, lambda)
    if nfa == 0
        fa_times = [];
        return
    end
    interval_dur = end_time - start_time;
    inter_arrival_times = -log(1 - rand(nfa, 1)) / lambda;
    normalized_times = cumsum(inter_arrival_times)';
    fa_times = start_time + (normalized_times / sum(inter_arrival_times)) * interval_dur;
    last_diff = fa_times(end) - fa_times(end-1);
    last_diff_updated = randf([0 last_diff]);
    fa_times(end) = fa_times(end - 1) + last_diff_updated;
end