clc;
clear;

N = 16;
M = 128;
mod_order = 4;
bps = log2(mod_order); %Bits per symbol

% OTFS Frame Length (Decimal)
L = N * M;

% OTFS Frame Length (Binary)
L_bits = L * bps;

% Codeword length (n) is 2^m -1, n = 2^m -1
% Message length (k) is n - m, k = n - m
m = 3;

% Hamming code
[H, G, n, k] = hammgen(m);
codeRate = k / n;
disp('CR = ');
disp(codeRate);

%Simulation Parameters
num_frames = 1000;

% Eb/No Parameters
ebno_vec_db = 0:2:16;

% Results data
ber_uncoded = zeros(size(ebno_vec_db));
ber_coded = zeros(size(ebno_vec_db));

% Error Calculator
error_rate_uncoded = comm.ErrorRate;
error_rate_coded = comm.ErrorRate;

num_blocks = floor(L_bits / n);
num_src_bits = num_blocks * k;

for idx_ebno = 1: length(ebno_vec_db)
  
  ebno_db = ebno_vec_db(idx_ebno);
  fprintf('Simulating Eb/No = %d dB ...\n', ebno_db);
  
  % Reset Error Calculators
  reset(error_rate_uncoded);
  reset(error_rate_coded);
  
  snr_uncoded_db = ebno_db + 10* log10(bps);
  snr_coded_db = ebno_db + 10* log10(bps) + 10*log10(codeRate);
  
  for frame_idx = 1: num_frames

    % Random message bits which will be sent
    random_bits_msg = randi([0 1], num_src_bits, 1);

    % -- UNCODED PATH -- %

    tx_symbols_uncoded = qammod(random_bits_msg, mod_order, InputType='bit', UnitAveragePower=true);
    rx_symbols_uncoded = awgn(tx_symbols_uncoded, snr_uncoded_db, 'measured');
    demod_bits_uncoded = qamdemod(rx_symbols_uncoded, mod_order, OutputType='bit', UnitAveragePower=true);

    ber_stats_uncoded = error_rate_uncoded(random_bits_msg, demod_bits_uncoded);

    % -- HAMMING PATH -- %

    % Encoded bits
    encoded_bits = [];
    for i = 1:num_blocks
      message_block = random_bits_msg((i-1)*k + 1 : i*k);
      cword = mod(G.' * message_block, 2);
      encoded_bits = [encoded_bits; cword];
    end

    % Padding Frame with zeroes
    num_bits_to_pad = L_bits - length(encoded_bits);
    padded_encoded_bits = [encoded_bits; zeros(num_bits_to_pad, 1)];

    tx_symbols_coded = qammod(padded_encoded_bits, mod_order, InputType='bit');
    rx_symbols_coded = awgn(tx_symbols_coded, snr_coded_db, 'measured');
    demod_bits_coded = qamdemod(rx_symbols_coded, mod_order, OutputType='bit', UnitAveragePower=true); 

    decoded_bits = [];
    for i = 1:num_blocks
      rec_cword =  demod_bits_coded((i-1)*n + 1: i*n);
      syndrome = mod(H * rec_cword, 2);
      error_idx = find(ismember(H.', syndrome.', 'rows'));
      if ~isempty(error_idx)
        rec_cword(error_idx) =  mod(rec_cword(error_idx) + 1, 2);
      end
      msg_block = rec_cword(end-k+1:end);
      decoded_bits = [decoded_bits; msg_block];
    end

    ber_stats_coded = error_rate_coded(random_bits_msg, decoded_bits);
  end
  ber_uncoded(idx_ebno) = ber_stats_uncoded(1);
  ber_coded(idx_ebno) = ber_stats_coded(1);
end

% Display Results

fprintf('\n--- Simulation Complete ---\n\n');
fprintf('Eb/No (dB) | Uncoded BER | Coded BER (Hamming)\n');
fprintf('-----------|-------------|---------------------\n');
for i_ebno = 1:length(ebno_vec_db)
    fprintf('%8.1f   | %11.5f | %17.5f\n', ebno_vec_db(i_ebno), ber_uncoded(i_ebno), ber_coded(i_ebno));
end
fprintf('\n');
