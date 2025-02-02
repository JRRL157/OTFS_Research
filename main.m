% Main procedure
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

%Multicarrier modulation schemes
otfs_str = "otfs";
otfs_zak_str = "otfs_zak";
otfs_wh_str = "otfs_walsh_hadamard";
otfs_dct_2_str = "otfs_dct_type_2";
ofdm_str = "ofdm";

mod_schemes = {'Classical OTFS', 'OTFS with Zak', 'OTFS with Walsh-Hadamard', 'OTFS with DCT Type II'};

%N = 16;
%M = 16;
MOD_SIZE = 4;

%spd = 100;
%delta_f = 15e3;
%fc = 4e9;

% SNR
SNR_step = 10; % Incremento de SNR em dB
SNR_values = 0:SNR_step:40; % Vetor de valores de SNR

% Number of Iterations
num = 500; % Número de iterações para cada SNR (simulações)

% Initialize BER_values object for all modulation types
BER_values = zeros(length(mod_schemes), length(SNR_values));

% Define selected (N, M, spd, fc, delta_f) tuples
simulation_params = [
    16, 16, 2.5, 3.5e9, 15e3;  % Low mobility, Sub-6 GHz
    16, 64, 2.5, 3.5e9, 15e3;  % Low mobility, Sub-6 GHz
    16, 16, 17, 28e9, 120e3;   % Medium mobility, mmWave
    16, 64, 17, 28e9, 120e3;   % Medium mobility, mmWave
    16, 16, 100, 6e9, 30e3;    % High-speed train, mid-band 5G
    16, 64, 100, 6e9, 30e3;    % High-speed train, mid-band 5G
    64, 64, 100, 6e9, 30e3;    % High-speed train, mid-band 5G
    16, 64, 7500, 20e9, 240e3  % Satellite LEO, Ka-band
    64, 64, 7500, 20e9, 240e3  % Satellite LEO, Ka-band    
];

for idx = 1:size(simulation_params)
    N = simulation_params(idx, 1);
    M = simulation_params(idx, 2);
    spd = simulation_params(idx, 3);
    fc = simulation_params(idx, 4);
    delta_f = simulation_params(idx, 5);
    
    %Título do plot
    plot_title = sprintf('Parâmetros: Num. Iterações=%d, N=%d, M=%d, Speed(m/s)=%d, delta_f(KHz)=%.2f, f_c(MHz)=%.2f', num, N, M, spd, (delta_f/1e3), (fc/1e9));

    for channel_model_selector = 1:3
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
    
        for mod_exp = 1:3
            % Loop sobre as diversas transformadas
            mod_size = MOD_SIZE^mod_exp;
            for type = 1:4
              % Loop sobre os valores de SNR
              for snr_idx = 1:length(SNR_values)
                  SNR_db = SNR_values(snr_idx); % Valor atual de SNR
                  total_errors = 0;
                  total_bits = 0;
            
                  % Simule para o valor atual de SNR
                  for i = 1:num
                      if (type == 1)
                        [x, x_hat] = otfs(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                      elseif (type == 2)
                        [x, x_hat] = otfs_zak(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                      elseif (type == 3)
                        [x, x_hat] = otfs_wh(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                      elseif (type == 4)
                        [x, x_hat] = otfs_dct_type_2(N, M, spd, fc, delta_f, SNR_db, mod_size, delays_arr, pdp_arr);
                      else
                        [x, x_hat] = cp_ofdm(N, M, spd, fc, delta_f, SNR_db, mod_size, optimized);
                      end
            
                      % Verifique se x e x_hat têm o mesmo comprimento
                      if length(x) ~= length(x_hat)
                          error('x e x_hat devem ter o mesmo comprimento.');
                      end
            
                      % Calcule o número de erros para esta execução
                      errors = sum(x ~= x_hat);
            
                      % Atualize os contadores totais
                      total_errors = total_errors + errors;
                      total_bits = total_bits + length(x);
                  end
            
                  % Calcule o BER para o valor atual de SNR
                  BER_values(type, snr_idx) = total_errors / total_bits;
              end
            end
            % Plot do gráfico BER vs. SNR
            fig = figure();
            hold on;
            for type = 1:4
                semilogy(SNR_values, BER_values(type, :), '-o', 'LineWidth', 1.5, 'DisplayName', mod_schemes{type});
            end
            grid on;
            xlabel('SNR (dB)');
            ylabel(sprintf('BER - [%d QAM, %s Channel Model]', mod_size, channel_model_name));
            legend ('show');
            title(plot_title);
            dateAndTime = datestr(now(),'yyyy_mmmm_dd_HH-MM-SS');
            fname = sprintf('%s%s_%s_%d-QAM.png', directory, dateAndTime, channel_model_name, mod_size);
            saveas(fig, fname, 'png');
        end
    end
end