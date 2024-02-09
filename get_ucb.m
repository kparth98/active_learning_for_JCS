function [ucb_idx, ucb_grid] = get_ucb(avg_rewards,visit_counts, beta, delta)
%GET_UCB Summary of this function goes here
%   Detailed explanation goes here
ucb = zeros(length(visit_counts), 1);
for i = 1:16
    if visit_counts(i) == 0
        ucb(i) = 1e8;
    else
        ucb(i) = avg_rewards(i) + beta * sqrt(log(1/delta) / visit_counts(i));
    end
end
[maxs, idxs] = max(ucb);
ucb_idx = idxs;
ucb_grid = ucb;
end

