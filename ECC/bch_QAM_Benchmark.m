?clc;
clear;

N = 16;
M = 16;
mod_order = 4;
bps = log2(mod_order); %Bits per symbol

% OTFS Frame Length (Decimal)
L = N * M;

% OTFS Frame Length (Binary)
L_bits = L * bps;

% Codeword length (n) is 2^m -1, n = 2^m -1
% Message length (k) is n - m, k = n - m
% Reference: https://uotechnology.edu.iq/dep-eee/lectures/4th/Communication/Information%20theory/2.pdf
% Reference 2: https://users.encs.concordia.ca/~msoleyma/ELEC6131/ELEC6131_2022/Lecture%20Notes/LECTURE%207.pdf

m = 9;
n = 2^m - 1;
k = 259;
codeRate = k / n;

% Hamming code
gp = bchgenpoly(n, k);

% Encoder and Decoder objects
enc = comm.BCHEncoder(n, k, gp);
dec = comm.BCHDecoder(n, k, gp);

t = bchnumerr(n, k);

%Simulation Parameters
num_frames = 1000;

% Eb/No Parameters
ebno_vec_db = 0:1:10;

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
    random_bits_msg_uncoded = randi([0 1], L_bits, 1);

    % -- UNCODED PATH -- %
    tx_symbols_uncoded = qammod(random_bits_msg_uncoded, mod_order, InputType='bit', UnitAveragePower=true);
    rx_symbols_uncoded = awgn(tx_symbols_uncoded, snr_uncoded_db, 'measured');
    demod_bits_uncoded = qamdemod(rx_symbols_uncoded, mod_order, OutputType='bit', UnitAveragePower=true);

    ber_stats_uncoded = error_rate_uncoded(random_bits_msg_uncoded, demod_bits_uncoded);

    % -- BCH PATH -- %
    random_bits_msg_coded = randi([0 1], num_src_bits, 1);

    % Encoded bits
    encoded_bits = [];
    for i = 1:num_blocks
      message_block = random_bits_msg_coded((i-1)*k + 1 : i*k);
      cword = enc(message_block);
      encoded_bits = [encoded_bits; cword];
    end

    % Padding Frame with zeroes
    num_bits_to_pad = L_bits - length(encoded_bits);
    padded_encoded_bits = [encoded_bits; zeros(num_bits_to_pad, 1)];
    
    interleaver_indices = randperm(length(padded_encoded_bits));
    interleaved_data = padded_encoded_bits(interleaver_indices);
    tx_symbols_coded = qammod(interleaved_data, mod_order, InputType='bit');

    rx_symbols_coded = awgn(tx_symbols_coded, snr_coded_db, 'measured');
   
    demod_bits_coded = qamdemod(rx_symbols_coded, mod_order, OutputType='bit', UnitAveragePower=true); 
    deinterleaved_bits = zeros(size(demod_bits_coded));
    deinterleaved_bits(interleaver_indices) = demod_bits_coded;

    decoded_bits = [];
    for i = 1:num_blocks
      rec_cword = deinterleaved_bits(1 + (i-1)*n : i * n);
      msg_block = dec(rec_cword);
      decoded_bits = [decoded_bits; msg_block];
    end

    ber_stats_coded = error_rate_coded(random_bits_msg_coded, decoded_bits);
  end
  ber_uncoded(idx_ebno) = ber_stats_uncoded(1);
  ber_coded(idx_ebno) = ber_stats_coded(1);

end

% Display Results
fprintf('\n--- Simulation Complete ---\n\n');
fprintf('Eb/No (dB) | Uncoded BER | Coded BER (BCH)\n');
fprintf('-----------|-------------|---------------------\n');
for i_ebno = 1:length(ebno_vec_db)
    fprintf('%8.1f   | %11.9f | %17.9f\n', ebno_vec_db(i_ebno), ber_uncoded(i_ebno), ber_coded(i_ebno));
end
fprintf('\n');
