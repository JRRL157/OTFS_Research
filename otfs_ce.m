function [dataIn, dataOut] = otfs_ce(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr)
  addpath('utils')
  Fn = fft(eye(N));
  Fm = fft(eye(M));
  Fn=Fn/norm(Fn);
  Fm=Fm/norm(Fm);

  %TODO: Please delete these later
  T=1/delta_f;
  c = 299792458;

  delay_resolution = 1/(M*delta_f);
  doppler_resolution = 1/(N*T);

  bps = log2(mod_size);

  m = bps;
  otfs_frame_length = N * M;
  L_bits = otfs_frame_length * bps;
  t = 40;
  n = 2^m -1;
  k = n - 2*t;
 codeRate = k / n;

    % Encoder Object
  rsEncoder = comm.RSEncoder( ...
      BitInput=true, ...
      CodewordLength=n, ...
      MessageLength=k);

  % Decoder Object
  rsDecoder = comm.RSDecoder( ...
      BitInput=true, ...
      CodewordLength=n, ...
      MessageLength=k);

  % Random message bits
  random_bits_msg = randi([0 1], L_bits, 1);

  dataIn = [];
  encoded_bits = [];
  num_blocks = floor((L_bits)/(k * m));

  for i=1:num_blocks
    msg = random_bits_msg(1 + (i-1)*(k*m): i*(k*m));
    enc_msg = rsEncoder(msg);
    dataIn = [dataIn; msg];
    encoded_bits = [encoded_bits; enc_msg];
  end

  bits2pad = L_bits - length(encoded_bits);
  encoded_bits = [encoded_bits; zeros(bits2pad, 1)];

  % Interleaver Configs
  rng('default'); 
  interleaver_state = randperm(length(encoded_bits));

  % disp(size(interleaver_state.'));
  % disp(size(encoded_bits));

  % Interleaved data
  % interleaved_data = intrlv(encoded_bits, interleaver_state.');
  interleaved_data = turbo_interleaver_3gpp(encoded_bits);

  % Modulate the integers to create symbols for the OTFS grid
  tx_info_symbols = qammod(interleaved_data, mod_size, InputType='bit');

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
  y = reshape(Y.', N*M, 1);
  %H = eye(N * M);  % Canal idealizado
  x_hat = (((H' * H + sigma_w_2*eye(size(H)))^(-1)) * H') * y;
  %x_hat = (G' * G + sigma_w_2)^(-1) * (G' * y); % Aqui ajustamos o canal simplificado
  if isnan(x_hat)
      dataIn = 0;
      dataOut = 0;
  else
    demod_bits = qamdemod(x_hat, mod_size, OutputType='bit');
    
    %De-Interleaver
    % disp(size(demod_bits));
    % disp(size(interleaver_state));
    % deinterleaved_data = deintrlv(demod_bits, interleaver_state.');
    deinterleaved_data = turbo_deinterleaver_3gpp(demod_bits);
    unpadded_bits = deinterleaved_data(1: L_bits - bits2pad);
    decoded_bits = rsDecoder(unpadded_bits);
    dataOut = decoded_bits;
  end
end

