function T_sparse = sparsify_delay_doppler(T, l_i, k_i, N, M)
    % Compute significant bins from the channel response
    delay_bins = unique(round(mod(l_i, N)) + 1); % Delay bins
    doppler_bins = unique(round(mod(k_i, M)) + 1); % Doppler bins

    % Ensure indices are valid (inside [1, N] and [1, M])
    delay_bins = delay_bins(delay_bins >= 1 & delay_bins <= N);
    doppler_bins = doppler_bins(doppler_bins >= 1 & doppler_bins <= M);

    % Create a mask preserving only dominant delay-Doppler bins
    mask = zeros(N, M);
    for d = delay_bins
        for k = doppler_bins
            disp([d, k]); % Debugging output
            mask(d, k) = 1;
        end
    end

    % Apply the sparsification mask
    T_sparse = T .* mask;
end