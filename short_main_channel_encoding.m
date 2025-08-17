% Main procedure
% THIS IS A MATLAB SCRIPT!
clc;
clear;

% Prefixo
directory = "images/";

% Simulation type (0 for OFDM, 1 for OTFS)
optimized = 0;

file_suffix = 'simulação.png';

% 3GPP Standard channel
delays_EVA = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
pdp_EVA = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];

delays_EPA = [0, 30, 70, 90, 110, 190, 410]*1e-9;
pdp_EPA = [0.0, -1.0, -2.0, -3.0, -8.0, -17.2, -20.8];

delays_ETU = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000]*1e-9;
pdp_ETU = [-1.0, -1.0, -1.0, 0.0, 0.0, 0.0, -3.0, -5.0, -7.0];

channel_model_name = ["EVA", "EPA", "ETU"];

mod_schemes = {'OTFS-DFT', 'RS OTFS-DFT', 'OTFS-Hadamard', 'RS OTFS-Hadamard'};

MOD_SIZE = 256;

% SNR
SNR_step = 10; % Incremento de SNR em dB
SNR_values = 30:SNR_step:50; % Vetor de valores de SNR

% Number of Iterations
num = 10;

% Initialize BER_values object for all modulation types
BER_values = zeros(length(mod_schemes), length(SNR_values));

%[1]Low-Complexity VLSI Architecture for OTFS Transceiver Under Multipath Fading Channel
%[2]Performance Evaluation of OTFS Over Measured V2V Channels at 60 GHz
%[3]A Robust and Low-Complexity Walsh-Hadamard Modulation for Doubly-Dispersive Channels

% Define selected (idx, N, M, spd, fc, delta_f) tuples
simulation_params = [
    1, 16, 16, 500, 6e9, 30e3; %[1]

    %2, 8, 8, 50, 60e9, 625000;%BW=5MHz, [2]
    %3, 8, 64, 50, 60e9, 78125;%BW=5MHz, [2]
    %4, 8, 8, 50, 60e9, 5e6;%BW=40MHz, [2]
    %5, 8, 64, 50, 60e9, 625e3;%BW=40MHz, [2]
    %6, 8, 8, 50, 60e9, 15e6; %BW=120MHz, [2]
    %7, 8, 64, 50, 60e9, 1.875e6; %BW=120MHz, [2]

    %8, 16, 64, 300, 5.9e9, 78.125e3;%BW=5MHz [3]
    %9, 4, 128, 300, 5.9e9, 78.125e3;%BW=10MHz [3]
];

for idx = 1:size(simulation_params)
    C = simulation_params(idx, 1);  
    N = simulation_params(idx, 2);
    M = simulation_params(idx, 3);
    spd = simulation_params(idx, 4);
    fc = simulation_params(idx, 5);
    delta_f = simulation_params(idx, 6);

    for channel_model_selector = 1:1
        delays_arr = [];
        pdp_arr = [];
        channel_model_name = "";
        switch(channel_model_selector)
            case 1
                delays_arr = delays_EVA;
                pdp_arr = pdp_EVA;
                channel_model_name = "EVA";
            case 2
                delays_arr = delays_EPA;
                pdp_arr = pdp_EPA;
                channel_model_name = "EPA";
            case 3
                delays_arr = delays_ETU;
                pdp_arr = pdp_ETU;
                channel_model_name = "ETU";
        end

        mod_size = MOD_SIZE;
        for type = 1:length(mod_schemes)
          % Loop sobre os valores de SNR
          for snr_idx = 1:length(SNR_values)
              SNR_db = SNR_values(snr_idx);
              total_errors = 0;
              total_bits = 0;

              % Simule para o valor atual de SNR
              for i = 1:num
                  disp(['Type = ', mod_schemes{type}, ', SNR #',num2str(SNR_db), ', Iteration #', num2str(i)]);
                  if (type == 1)
                    [x, x_hat] = otfs(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                  elseif (type == 2)
                    [x, x_hat] = otfs_rs(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                  elseif (type == 3)
                    [x, x_hat] = otfs_wh(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                  elseif (type == 4)
                    [x, x_hat] = otfs_wh_rs(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                  else
                    [x, x_hat] = cp_ofdm(N, M, spd, fc, delta_f, SNR_db, mod_size, optimized);
                  end

                  % Verifique se x e x_hat têm o mesmo comprimento
                  if length(x) ~= length(x_hat)
                      %disp(['Temperature is:' num2str(UU(90)) 'After: ' num2str(timeInMinutes) ' minutes']);
                      disp('x e x_hat devem ter o mesmo comprimento.');
                      disp(['size(x):' num2str(length(x))]);
                      disp(['size(x_hat)' num2str(length(x_hat))]);
                      continue;
                  end

                  % Calcule o número de erros para esta execução
                  errors = sum(x ~= x_hat);
        
                  % Atualize os contadores totais
                  total_errors = total_errors + errors;
                  total_bits = total_bits + length(x);
              end

              % Calcule o BER para o valor atual de SNR
              BER_values(type, snr_idx) = total_errors / total_bits;

              % Write BER values to a file
              filename = 'ber_values.txt';
              file = fopen(filename, 'a'); % Open file in append mode
              if file == -1
                  error('Could not open file for writing.');
              end
              dateAndTime = datestr(now(), 'yyyy_mmmm_dd_HH-MM-SS');
              fprintf(file, 'Date: %s: ,Scenario Index: %d, Type: %s, Modulation: %d-QAM, SNR: %f, BER: %e\n', dateAndTime, C, mod_schemes{type}, mod_size, SNR_db, BER_values(type, snr_idx));
              fclose(file);
          end
        end
    end
end