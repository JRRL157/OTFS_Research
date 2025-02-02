function [x, x_hat2] = otfs_simulation(N, M, spd, SNR_db, mod_size)
  Fn = dftmtx(N);
  Fn=Fn/norm(Fn);

  %TODO: Please delete these later
  delta_f = 15e3;
  T=1/delta_f;
  fc = 4e9;
  c = 299792458;

  delay_resolution = 1/(M*delta_f);
  doppler_resolution = 1/(N*T);

  % Generating OTFS frame
  N_syms_per_frame = N*M;

  random_syms = randi([0, mod_size-1], N_syms_per_frame,1);
  tx_info_symbols = qammod(random_syms, mod_size);

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
    endfor
  endfor

  X_til = X*Fn';
  s = reshape(X_til, 1, N*M);

  % Channel
  max_ue_spd_mps = 100 / 3.6;
  nu_max = (max_ue_spd_mps*fc)/c;
  k_max = nu_max/doppler_resolution;

  % 3GPP Standard channel
  delays_EVA = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
  pdp_EVA = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];

  % Generate standard channel parameters
  delays = delays_EVA;
  pdp=pdp_EVA;
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
    endfor
  endfor

  G = zeros(N*M, N*M);
  for q=0:N*M-1
    for ell=0:delay_spread
      if(q >= ell)
        G(q+1,q-ell+1)=gs(ell+1,q+1);
      end
    endfor
  endfor

  H_til = P*G*P.';
  H=kron(Im,Fn)*(P' * G * P)*kron(Im,Fn');

  % Generate r by passing the Tx signal through the channel
  r=G*s.';

  % Add AWGN
  Es = mean(abs(qammod(0:mod_size -1, mod_size).^2));

  SNR=10.^(SNR_db/10);

  sigma_w_2=Es/SNR;
  noise = sqrt(sigma_w_2 / 2)*(randn(N*M, 1) + 1i*randn(N*M, 1));

  r+= noise;

  % OTFS demodulation
  Y_til = reshape(r, M, N);
  Y = Y_til * Fn;

  % OTFS delay-doppler LMMSE detection
  y = reshape(Y.', N*M, 1);
  x_hat = (H' * H + sigma_w_2)^(-1) * (H' * y);
  x_hat = qamdemod(x_hat, mod_size);
  x_hat2 = qammod(x_hat, mod_size);
end
