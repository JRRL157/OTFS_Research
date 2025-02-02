function [DFT] = dft_matrix(N, M)
    % DFT_MATRIX_OTFS computes the DFT matrix optimized for OTFS.
    %
    % INPUT:
    %   N - Length of the signal (must be divisible by M).
    %   M - Number of delay taps.
    %
    % OUTPUT:
    %   DFT - DFT matrix optimized for OTFS of size N x N.

    % Ensure N is divisible by M
    if mod(N, M) ~= 0
        error('N must be divisible by M.');
    end

    % Number of Doppler bins
    L = N / M;

    % Initialize the DFT matrix
    DFT = zeros(N, N); % Complex entries

    % Populate the DFT matrix
    for k = 0:L-1 % Doppler index
        for m = 0:M-1 % Delay index
            for n = 0:N-1 % Time index
                % Check if the time index belongs to the current delay segment
                if mod(n, M) == m
                    DFT(k + m * L + 1, n + 1) = exp(-1j * 2 * pi * k * n / M);
                end
            end
        end
    end
end
