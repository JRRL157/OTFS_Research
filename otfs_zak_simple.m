function [x, x_hat2] = otfs_zak_simple(N, M, spd, SNR_db, mod_size)
  %TODO: Please delete these later
  delta_f = 15e3;
  T=1/delta_f;
  fc = 4e9;
  c = 299792458;

  delay_resolution = 1/(M*delta_f);
  doppler_resolution = 1/(N*T);

  Zak = exp(-1i * 2 * pi * (0:M-1)' * (0:N-1) / M);

  % Gerando o quadro OTFS
  N_syms_per_frame = N * M;

  random_syms = randi([0, mod_size-1], N_syms_per_frame, 1);
  tx_info_symbols = qammod(random_syms, mod_size);

  X = reshape(tx_info_symbols, M, N);
  x = reshape(X.', N * M, 1);

  % Modulação OTFS
  P = zeros(N * M, N * M);
  for j = 1:N
    for i = 1:M
      E = zeros(M, N);
      E(i, j) = 1;
      P((j-1)*M+1:j*M, (i-1)*N+1:i*N) = E;
    end
  end

  % Aplicando transformada de Zak
  s = zak_transform(X, M, N);
  s = reshape(s, 1, M * N);

  % Canal simplificado (sem desvanecimento complexo)
  max_ue_spd_mps = spd / 3.6;
  nu_max = (max_ue_spd_mps * fc) / c;
  k_max = nu_max / doppler_resolution;

  % Gerando parâmetros do canal (sem a complexidade de múltiplos taps)
  g_i = sqrt(0.5) * (randn(1, 1) + 1i * randn(1, 1)); % Canal de Rayleigh simples
  k_i = k_max * cos(2 * pi * rand(1, 1)); % Doppler simplificado

  % Geração da matriz do canal G
  G = g_i * exp(1i * 2 * pi / (N * M) * k_i * (0:N * M - 1));

  % Gerar r ao passar o sinal Tx pelo canal
  r = s .* G;

  Es = mean(abs(qammod(0:mod_size - 1, mod_size).^2));
  SNR = 10^(SNR_db / 10);
  sigma_w_2 = Es / SNR;
  noise = sqrt(sigma_w_2 / 2) * (randn(N * M, 1) + 1i * randn(N * M, 1));

  % Canal ideal (apenas AWGN, sem multi-path)
  r = s + noise.';

  % Demodulação OTFS com Zak
  Y_zak = reshape(r, M, N);
  Y = zak_inverse(Y_zak, M, N);

  % Detecção LMMSE
  H = eye(N * M);  % Canal idealizado
  x_hat = (H' * H + sigma_w_2)^(-1) * (H' * Y(:));
  x_hat = qamdemod(x_hat, mod_size);
  x_hat2 = qammod(x_hat, mod_size);
end

