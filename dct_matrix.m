function [DCT] = dct_matrix(N)
    % DCT_MATRIX computes the Discrete Cosine Transform Type-2 matrix.
    %
    % INPUT:
    %   N - Size of the DCT matrix (NxN).
    %
    % OUTPUT:
    %   DCT - DCT Type-2 matrix of size N x N.

    DCT = zeros(N, N);

    % Populate the DCT matrix
    for k = 0:N-1
        for n = 0:N-1
            if k == 0
                alpha = sqrt(1 / N);
            else
                alpha = sqrt(2 / N);
            end
            DCT(k+1, n+1) = alpha * cos(pi * (2*n + 1) * k / (2 * N));
        end
    end
end
