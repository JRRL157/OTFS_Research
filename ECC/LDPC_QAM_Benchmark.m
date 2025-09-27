clc;
clear;
close all;

%==========================================================================
%% 1. Simulation Parameters
%==========================================================================

% --- Modulation Parameters ---
mod_order = 256;
bps = log2(mod_order); % Bits per symbol

% --- Eb/No Parameters ---
ebno_vec_db = 0:1:10; % Eb/No range in dB

% --- LDPC Code Parameters (DVB-S2 Standard) ---
H = dvbs2ldpc(1/2, 'sparse');
cfgLDPCEncoder = ldpcEncoderConfig(H);
cfgLDPCDecoder = ldpcDecoderConfig(H);

k = cfgLDPCEncoder.NumInformationBits; 
n = cfgLDPCEncoder.BlockLength; 
codeRate = k / n;

% --- Simulation Control ---
% We will simulate one LDPC block per frame for simplicity
num_info_bits_per_frame = k;
num_frames_per_point = 1000; % Increase for smoother curves

% --- Results Storage ---
ber_uncoded = zeros(size(ebno_vec_db));
ber_coded = zeros(size(ebno_vec_db));

fprintf('Starting Simulation...\n');
fprintf('Modulation: %d-QAM\n', mod_order);
fprintf('Coding: DVB-S2 LDPC(%d,%d) vs. Uncoded\n', n, k);
fprintf('Info Bits per Frame: %d\n\n', num_info_bits_per_frame);

%==========================================================================
%% 2. Simulation Loop
%==========================================================================

for idx_ebno = 1:length(ebno_vec_db)
    ebno_db = ebno_vec_db(idx_ebno);
    fprintf('Simulating Eb/No = %d dB ...\n', ebno_db);

    % Create error rate calculators for this point
    error_rate_uncoded = comm.ErrorRate;
    error_rate_coded = comm.ErrorRate;

    % Calculate SNR from Eb/No for both systems.
    % The coded system's SNR is lower due to the "price of redundancy".
    snr_uncoded_db = ebno_db + 10*log10(bps);
    snr_coded_db = ebno_db + 10*log10(bps) + 10*log10(codeRate);

    % The noise variance is the reciprocal of the linear SNR
    noiseVar = 10.^(-snr_coded_db/10);

    for frame_idx = 1:num_frames_per_point

        % --- UNCODED PATH ---
        % Transmit the same number of information bits as the coded system.
        uncoded_src_bits = randi([0 1], num_info_bits_per_frame, 1);
        tx_symbols_uncoded = qammod(uncoded_src_bits, mod_order, 'InputType', 'bit', 'UnitAveragePower', true);
        rx_symbols_uncoded = awgn(tx_symbols_uncoded, snr_uncoded_db, 'measured');
        demod_bits_uncoded = qamdemod(rx_symbols_uncoded, mod_order, 'OutputType', 'bit', 'UnitAveragePower', true);
        ber_stats_uncoded = error_rate_uncoded(uncoded_src_bits, demod_bits_uncoded);

        % --- LDPC CODED PATH ---
        coded_src_bits = randi([0 1], k, 1);
        encoded_bits = ldpcEncode(coded_src_bits, cfgLDPCEncoder);
        tx_symbols_coded = qammod(encoded_bits, mod_order, 'InputType', 'bit', 'UnitAveragePower', true);
        % Noise
        rx_symbols_coded = awgn(tx_symbols_coded, snr_coded_db, 'measured');

        % Soft Demodulation 
        % We need LLRs, not hard bits. NoiseVariance is essential.
        demod_llrs = qamdemod(rx_symbols_coded, mod_order, ...
            'OutputType', 'approxllr', ...
            'NoiseVariance', noiseVar, ...
            'UnitAveragePower', true);

        % 6. LDPC Decode
        decoded_bits = double(ldpcDecode(demod_llrs, cfgLDPCDecoder, 50));
        
        % Error rate count
        ber_stats_coded = error_rate_coded(coded_src_bits, decoded_bits);
    end

    % Store final BER for this Eb/No point
    ber_uncoded(idx_ebno) = ber_stats_uncoded(1);
    ber_coded(idx_ebno) = ber_stats_coded(1);
end

%==========================================================================
%% 3. Display and Plot Results
%==========================================================================

fprintf('\n--- Simulation Complete ---\n\n');
fprintf('Eb/No (dB) | Uncoded BER | Coded BER (LDPC)\n');
fprintf('-----------|-------------|---------------------\n');
for i_ebno = 1:length(ebno_vec_db)
    fprintf('%8.1f   | %11.15f | %17.15f\n', ebno_vec_db(i_ebno), ber_uncoded(i_ebno), ber_coded(i_ebno));
end
fprintf('\n');
