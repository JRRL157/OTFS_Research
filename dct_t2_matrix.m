function [DCT] = dct_t2_matrix(N, M)
    % DCT_MATRIX_OTFS computes the DCT Type-2 matrix optimized for OTFS.
    %
    % INPUT:
    %   N - Length of the signal (must be divisible by M).
    %   M - Number of delay taps.
    %
    % OUTPUT:
    %   DCT - DCT Type-2 matrix optimized for OTFS of size N x N.

    % Ensure N is divisible by M
    if mod(N, M) ~= 0
        error('N must be divisible by M.');
    end

    % Number of Doppler bins
    L = N / M;

    % Initialize the DCT matrix
    DCT = zeros(N, N);

    % Populate the DCT matrix
    for k = 0:L-1 % Doppler index
        for m = 0:M-1 % Delay index
            for n = 0:N-1 % Time index
                % Check if the time index belongs to the current delay segment
                if mod(n, M) == m
                    alpha = sqrt(1 / M) * (k == 0) + sqrt(2 / M) * (k > 0);
                    DCT(k + m * L + 1, n + 1) = alpha * cos(pi * (2 * mod(n, M) + 1) * k / (2 * M));
                end
            end
        end
    end
end

