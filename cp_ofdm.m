function [x, x_hat2] = cp_ofdm(N, M, spd, SNR_db, mod_size)
  Fn = dftmtx(N);
  Fn = Fn / norm(Fn);

  % Parâmetros do canal
  delta_f = 15e3;
  T = 1 / delta_f;
  fc = 4e9;
  c = 299792458;

  delay_resolution = 1 / (M * delta_f);
  doppler_resolution = 1 / (N * T);

  % Geração do frame OFDM
  N_syms_per_frame = N * M;

  random_syms = randi([0, mod_size - 1], N_syms_per_frame, 1);
  tx_info_symbols = qammod(random_syms, mod_size);

  X = reshape(tx_info_symbols, M, N);
  x = reshape(X.', N * M, 1);

  % Modulação CP-OFDM
  s = ifft(X, [], 1);
  cp_length = round(0.1 * N);
  s_cp = [s(end-cp_length+1:end, :); s];
  s = s_cp(:);

  % Canal
  max_ue_spd_mps = spd / 3.6;
  nu_max = (max_ue_spd_mps * fc) / c;
  k_max = nu_max / doppler_resolution;

  delays_EVA = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
  pdp_EVA = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];

  pdp_linear = 10.^(pdp_EVA / 10);
  pdp_linear = pdp_linear / sum(pdp_linear);
  taps = length(pdp_linear);

  g_i = sqrt(pdp_linear) .* (sqrt(1/2) * (randn(1, taps) + 1i * randn(1, taps)));
  l_i = round(delays_EVA ./ delay_resolution);
  k_i = (k_max * cos(2 * pi * rand(1, taps)));

  G = zeros(N * M, N * M);
  for q = 0:N*M-1
    for ell = 0:max(l_i)
      if (q >= ell)
        G(q + 1, q - ell + 1) = g_i(ell + 1);
      end
    end
  end

  % Recepção com ruído
  r = conv(s, G(:, 1), 'same');

  Es = mean(abs(qammod(0:mod_size - 1, mod_size).^2));
  SNR = 10.^(SNR_db / 10);
  sigma_w_2 = Es / SNR;
  noise = sqrt(sigma_w_2 / 2) * (randn(size(r)) + 1i * randn(size(r)));

  r = r + noise;

  % Remoção do CP e demodulação
  rx_signal_cp = reshape(r, cp_length + N, []).';
  rx_signal_no_cp = rx_signal_cp(:, cp_length + 1:end);
  Y = fft(rx_signal_no_cp, [], 1);

  % Equalização
  H_channel = fft(G(:, 1), N);
  X_hat = Y ./ H_channel.';

  % Decodificação
  x_hat = qamdemod(X_hat(:), mod_size);
  x_hat2 = qammod(x_hat, mod_size);
end