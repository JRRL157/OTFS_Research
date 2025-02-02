function X = zak_inverse(Z, M, N)
    % Calcula a transformada inversa de Zak
    X = zeros(M, N);
    for m = 1:M
        for n = 1:N
            X(m, n) = sum(Z(:, n) .* exp(1j * 2 * pi * (0:M-1)' * (m-1) / M));
        end
    end
    X = X / M; % Normalização
end
