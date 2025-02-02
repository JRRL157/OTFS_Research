function Z_matrix = dzt_matrix(N)
    % Função para calcular a matriz de transformação de Zak de tamanho NxN.
    %
    % Entrada:
    %   N - Tamanho da matriz (inteiro positivo).
    % Saída:
    %   Z_matrix - Matriz de transformação de Zak (NxN).

    % Inicializa a matriz de Zak
    Z_matrix = zeros(N, N);

    % Preenche a matriz de Zak
    for m = 0:N-1
        for n = 0:N-1
            % Índices ajustados para formato de Zak
            % Multiplicação por exp(2*pi*i*m*n/N) faz a modulação
            Z_matrix(m+1, n+1) = exp(-2 * pi * 1i * (mod(m * n, N)) / N);
        end
    end
end
