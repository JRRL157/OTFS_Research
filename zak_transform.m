function Z = zak_transform(X, M, N)
    % Calcula a transformada de Zak
    Z = zeros(M, N);
    for m = 1:M
        for n = 1:N
            Z(m, n) = sum(X(:, n) .* exp(-1j * 2 * pi * (m-1) * (0:M-1)' / M));
        end
    end
end
