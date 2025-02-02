function [DWHT] = dwht_matrix(N, M)
    % DWHT_MATRIX_OTFS computes the DWHT matrix optimized for OTFS.
    %
    % INPUT:
    %   N - Length of the signal (must be divisible by M and a power of 2).
    %   M - Number of delay taps.
    %
    % OUTPUT:
    %   DWHT - DWHT matrix optimized for OTFS of size N x N.

    % Ensure N is divisible by M and is a power of 2
    if mod(N, M) ~= 0 || mod(log2(N), 1) ~= 0
        error('N must be divisible by M and a power of 2.');
    end

    % Number of Doppler bins
    L = N / M;

    % Initialize the DWHT matrix
    DWHT = zeros(N, N);

    % Populate the DWHT matrix
    for k = 0:L-1 % Doppler index
        for m = 0:M-1 % Delay index
            for n = 0:N-1 % Time index
                % Check if the time index belongs to the current delay segment
                if mod(n, M) == m
                    DWHT(k + m * L + 1, n + 1) = hadamard(M)(mod(n, M) + 1, k + 1);
                end
            end
        end
    end

    % Normalize (optional)
    DWHT = DWHT / sqrt(N);
end
