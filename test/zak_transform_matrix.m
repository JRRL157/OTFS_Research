function [Z] = zak_transform_matrix(N, M)
    % ZAK_TRANSFORM_MATRIX computes the Zak transform matrix Z for a signal
    % of length N with M delay taps.
    %
    % INPUT:
    %   N - Length of the signal (must be divisible by M).
    %   M - Number of delay taps.
    %
    % OUTPUT:
    %   Z - Zak transform matrix of size N x N.

    % Ensure N is divisible by M
    if mod(N, M) ~= 0
        error('N must be divisible by M.');
    end

    % Number of Doppler bins
    L = N / M;

    % Initialize the Zak transform matrix
    Z = zeros(N, N); % Complex entries, automatically cast when adding exp terms

    % Populate the Zak transform matrix
    for k = 0:L-1 % Doppler index
        for m = 0:M-1 % Delay index
            for n = 0:N-1 % Time index
                % Check if the time index belongs to the current delay segment
                if mod(n, M) == m
                    Z(k + m * L + 1, n + 1) = exp(-1j * 2 * pi * k * n / M);
                end
            end
        end
    end
end

