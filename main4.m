clc; clear; close all;pkg load communications;

% Parameters
N = 16; % Number of subcarriers
M = 16; % Number of symbols
mod_size = 4; % QPSK
num_frames = 100; % Number of frames for averaging
SNR_dB = 0:5:30; % SNR range in dB

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

    for frame = 1:num_frames
        % Transmit symbols
        tx_symbols = randi([0, mod_size-1], M * N, 1);
        tx_modulated = qam_mod(tx_symbols);

        % Reshape to OTFS grid
        tx_grid = reshape(tx_modulated, M, N);

        % OTFS Modulation (ISFFT)
        tx_time = sqrt(N) * ifft(ifft(tx_grid, [], 2), [], 1);
        tx_time = tx_time(:); % Serialize

        % Channel: Doubly dispersive (simple model)
        delay_taps = [0 2 4];
        doppler_shifts = [0 0.1 0.2]; % Normalized Doppler shifts
        path_gains = [1 0.8 0.6];

        P = length(path_gains); % Number of channel
        H_eff = zeros(M * N, M * N);
        for p = 1:P
            % Delay effect: Circular shift in time
            delay_matrix = circshift(eye(M * N), delay_taps(p), 1);

            % Doppler effect: Diagonal phase shifts
            doppler_matrix = diag(exp(1j * 2 * pi * doppler_shifts(p) * (0:M * N - 1)));

            % Combine delay and Doppler effects with path gain
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
        rx_grid = fft(fft(rx_grid, [], 1), [], 2) / sqrt(N);

        % Detection
        %rx_symbols = rx_grid(:);
        %rx_symbols = (H_eff' * H_eff + regularization * eye(size(H_eff))) \ (H_eff' * rx_grid(:));
        %rx_symbols = pinv(H_sparse) * rx_grid(:);
        regularization = noise_var; % Regularization parameter (depends on noise power)
        rx_symbols = (H_eff' * H_eff + regularization * eye(size(H_eff))) \ (H_eff' * rx_grid(:));
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

