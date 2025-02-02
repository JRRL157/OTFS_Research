clc; clear; close all;pkg load communications;

% Parameters
N = 16; % Number of subcarriers
M = 16; % Number of symbols
mod_size = 4; % QPSK
num_iterations = 1000; % Number of frames for averaging
SNR_dB = 0:1:20; % SNR range in dB

% Derived parameters
qam_mod = @(x) qammod(x, mod_size);
qam_demod = @(x) qamdemod(x, mod_size);

% Pre-allocate BER
BER = zeros(size(SNR_dB));

% Loop over SNR values
for snr_idx = 1:length(SNR_dB)
    SNR = 10^(SNR_dB(snr_idx) / 10);
    noise_var = 1 / SNR; % Noise variance for UnitAveragePower
    total_errors = 0;
    total_bits = 0;

    for frame = 1:num_iterations
        % Transmit symbols
        tx_symbols = randi([0, mod_size-1], M * N, 1);
        tx_modulated = qam_mod(tx_symbols);

        % Reshape to OTFS grid
        tx_grid = reshape(tx_modulated, M, N);

        % OTFS Modulation (ISFFT)
        tx_time = sqrt(N) * zak_transform(tx_grid, M, N);
        tx_time = tx_time(:); % Serialize

        % Simplified channel with single tap for debugging
        path_gains = [1 0.8]; % Fixed gains
        delay_taps = [0 1]; % Fixed delays
        doppler_shifts = [0 0.05]; % Fixed Doppler

        % Construct simplified H_eff
        H_eff = zeros(M * N, M * N);
        for p = 1:length(path_gains)
            delay_matrix = circshift(eye(M * N), delay_taps(p), 1);
            doppler_matrix = diag(exp(1j * 2 * pi * doppler_shifts(p) * (0:M * N - 1)));
            H_eff = H_eff + path_gains(p) * delay_matrix * doppler_matrix;
        end

        rx_time = zeros(size(tx_time));
        for p = 1:length(delay_taps)
            delayed_signal = [zeros(delay_taps(p), 1);
            tx_time(1:end-delay_taps(p))];
            doppler_effect = exp(1j * 2 * pi * doppler_shifts(p) * (0:length(tx_time)-1)');
            rx_time = rx_time + path_gains(p) * delayed_signal .* doppler_effect;
        end

        % Add noise
        noise = sqrt(noise_var / 2) * (randn(size(rx_time)) + 1j * randn(size(rx_time)));
        rx_time = rx_time + noise;

        % OTFS Demodulation (SFFT)
        rx_grid = reshape(rx_time, M, N);
        rx_grid = zak_inverse(rx_grid, M, N) / sqrt(N);

        % Detection
        rx_symbols = (H_eff' * H_eff + noise_var * eye(size(H_eff))) \ (H_eff' * rx_grid(:));
        rx_demodulated = qam_demod(rx_symbols);

        % BER Calculation
        total_errors = total_errors + sum(rx_demodulated ~= tx_symbols);
        total_bits = total_bits + numel(tx_symbols) * log2(mod_size);
    end

    % Calculate BER for this SNR
    BER(snr_idx) = total_errors / total_bits;
end

% Plot BER curve
semilogy(SNR_dB, BER, '-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('OTFS BER over Doubly Dispersive Channel');

