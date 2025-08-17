function [x, x_hat2] = otfs_wh_rs(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr)
  WHn = hadamard(N);
  WHm = hadamard(M);
  WHn=WHn/norm(WHn);
  WHm=WHm/norm(WHm);

  %TODO: Please delete these later
  T=1/delta_f;
  c = 299792458;

  delay_resolution = 1/(M*delta_f);
  doppler_resolution = 1/(N*T);

  % Generating OTFS frame
  N_syms_per_frame = N*M;

  random_syms = randi([0 mod_size-1], [N_syms_per_frame 1]);

  % RS Encoder
  n_rs = 255; %Codeword length
  k_rs = 223; %Message length
  encoded_msg = rsenc(random_syms, n_rs, k_rs);

  tx_info_symbols = qammod(encoded_msg, mod_size);

  X = reshape(tx_info_symbols, M, N);
  x = reshape(X.', N*M, 1);

  % OTFS Modulation
  Im = eye(M);

  P = zeros(N*M, N*M);
  for j=1:N
    for i=1:M
      E=zeros(M,N);
      E(i,j)=1;
      P((j-1)*M+1:j*M,(i-1)*N+1:i*N)=E;
    end
  end

  X_tf = WHm*X*WHn';
  X_til = WHm' * X_tf;
  s = reshape(X_til, 1, N*M);

  % Channel
  max_ue_spd_mps = spd / 3.6;
  nu_max = (max_ue_spd_mps*fc)/c;
  k_max = nu_max/doppler_resolution;

  % Generate standard channel parameters
  delays = delays_arr;
  pdp=pdp_arr;
  pdp_linear = 10.^(pdp/10);
  pdp_linear = pdp_linear/sum(pdp_linear);
  taps=length(pdp_linear);

  g_i = sqrt(pdp_linear).*(sqrt(1/2)*(randn(1,taps)+1i*randn(1,taps)));
  l_i = round(delays./delay_resolution);
  k_i = (k_max*cos(2*pi*rand(1,taps))); %%%%% --------- %%%%%%

  % Generate discrete delay-time channel coefficients and matrix_type
  z = exp(1i*2*pi/N/M);
  delay_spread = max(l_i);

  gs = zeros(delay_spread+1, N*M);
  for q=0:N*M-1
    for i=1:taps
      gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1) + g_i(i)*z^(k_i(i)*(q-l_i(i)));
    end
  end

  G = zeros(N*M, N*M);
  for q=0:N*M-1
    for ell=0:delay_spread
      if(q >= ell)
        G(q+1,q-ell+1)=gs(ell+1,q+1);
      end
    end
  end

  H_til = P*G*P.';
  H=kron(Im,WHn)*(P' * G * P)*kron(Im,WHn');

  % Generate r by passing the Tx signal through the channel
  r=G*s.';

  % Add AWGN
  Es = mean(abs(x).^2);
  SNR=10.^(SNR_db/10);

  sigma_w_2 = (Es / SNR); % Normalize by number of symbols
  noise = sqrt(sigma_w_2 / 2)*(randn(N*M, 1) + 1i*randn(N*M, 1));

  r = r + noise;

  % OTFS demodulation
  Y_til = reshape(r, M, N);
  Y_tf = WHm * Y_til;
  Y = WHm' * Y_tf * WHn;

  % OTFS delay-doppler LMMSE detection
  y = reshape(Y.', N*M, 1);
  %H = eye(N * M);  % Canal idealizado
  x_hat_symbols = (((H' * H + sigma_w_2*eye(size(H)))^(-1)) * H') * y;
  %x_hat = (G' * G + sigma_w_2)^(-1) * (G' * y); % Aqui ajustamos o canal simplificado
  if isnan(x_hat)
    x_hat = 0;
    x_hat2 = 0;
  else
    received_encoded_msg = qamdemod(x_hat_symbols, mod_size);
    [decoded_msg, num_errors] = rsdec(received_encoded_msg, n_rs, k_rs);
    x_hat = qamdemod(decoded_msg, mod_size);
    x_hat2 = qammod(decoded_msg, mod_size);
  end
end

