function [dataIn, dataOut] = otfs_ce(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr)
  Fn = fft(eye(N));
  Fm = fft(eye(M));
  Fn=Fn/norm(Fn);
  Fm=Fm/norm(Fm);

  %TODO: Please delete these later
  T=1/delta_f;
  c = 299792458;

  delay_resolution = 1/(M*delta_f);
  doppler_resolution = 1/(N*T);

  % --- Define a simple Rate 1/2 Convolutional Code ---
  k = log2(mod_size);
  constrlen = 7; % Standard constraint length
  tPoly = poly2trellis(constrlen, [171 133]); % Standard generator polynomials
  codeRate = 1/2;

  % --- Generate, Encode, and Modulate Data ---
  num_symbols_per_frame = N * M;
  num_encoded_bits = num_symbols_per_frame * k;
  num_message_bits = num_encoded_bits * codeRate;

  % Generate original message bits (as a column vector)
  dataIn = randi([0 1], num_message_bits - (constrlen-1), 1);

  % Add tail bits to flush the encoder (for robust decoding)
  dataIn_terminated = [dataIn; zeros(constrlen-1, 1)];

  % Encode the data
  dataEnc = convenc(dataIn_terminated, tPoly);

  % Group encoded bits into integers
  dataSymbolsIn = bit2int(dataEnc, k);

  % Modulate the integers to create symbols for the OTFS grid
  tx_info_symbols = qammod(dataSymbolsIn, mod_size, 'UnitAveragePower', true);

  X = reshape(tx_info_symbols, M, N);
  x = X(:);

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

  X_tf = Fm*X*Fn';
  X_til = Fm' * X_tf;
  s = X_til(:).';

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
  H=kron(Im,Fn)*(P' * G * P)*kron(Im,Fn');

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
  Y_tf = Fm * Y_til;
  Y = Fm' * Y_tf * Fn;

  % OTFS delay-doppler LMMSE detection
  y = Y(:);
  
  % --- LMMSE Equalizer with CORRECT End-to-End Channel Matrix H ---
  % This is the key fix. We build H that maps X_vec to y_vec.
  H_dd = zeros(N*M, N*M);
  for i = 1:taps
      for m_prime = 1:M
          for n_prime = 1:N
              l_cycl = mod(m_prime - 1 - l_i(i), M);
              H_dd = H_dd + g_i(i) * circshift(diag(Fn(:, n_prime)), [0, l_cycl]) * diag(exp(1i*2*pi*k_i(i)*(0:M-1)/M)) * circshift(diag(Fn(:, n_prime)'), [l_cycl, 0]);
          end
      end
  end
  
  % A simplified but more common H matrix formulation for LMMSE
  % This relates the DD symbols X to the time domain signal r
  % For equalization, we need to bring r back to DD domain first.
  % Let's use the effective channel matrix H_eff that relates X(:) to y(:)
  % H_eff = F_M' * diag(h_circ) * F_M, where h_circ is from G
  % This is complex. Let's use a well-known approximation.
  
  % Let's revert to your original H calculation but fix the transforms
  % This is the most direct fix to your existing code structure

  % 1. Create the permutation matrix P
  P = zeros(N*M, N*M);
  for j=1:N
      for i=1:M
          E=zeros(M,N);
          E(i,j)=1;
          P((j-1)*M+1:j*M,(i-1)*N+1:i*N)=E;
      end
  end

  % 2. Calculate the channel in the Delay-Doppler domain
  H_til = P * G * P';

  % 3. Form the full transform matrix and calculate the final H
  F_MN = kron(Fm, Fn);
  H = F_MN * H_til * F_MN'; % Correct effective channel matrix

  % 4. Perform LMMSE equalization
  % Using '\' is more stable and efficient than inv()
  x_hat_symbols = (H' * H + sigma_w_2 * eye(size(H))) \ (H' * y);

  % --- Demodulate and Decode ---
  if any(isnan(x_hat_symbols))
      dataOut = -1;
  else
    dataSymbOut = qamdemod(x_hat_symbols, mod_size, 'UnitAveragePower', true);
    codedDataOut = int2bit(dataSymbOut, k);
    dataOut = vitdec(codedDataOut, tPoly, traceBack, 'term', 'hard');
    dataOut = dataOut(1:length(dataIn));
  end
end

