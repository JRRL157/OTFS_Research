clc;
clear;

server_ip = '127.0.0.1';
server_port = 5016;

cmd_code = 0;
M = 32;
CR = 1/2;
K = ceil(M * CR);

data = [cmd_code, M, CR];

Y = NR_LDPC_client(server_ip, server_port, data, 1);

cmd_code = 1;
random_bits_msg = randi([0 1], K, 1).';

data = [cmd_code, K, random_bits_msg];

Y = NR_LDPC_client(server_ip, server_port, data, M);

disp("Original: " + num2str(random_bits_msg));
disp("Encoded: " + num2str(Y));

cmd_code = 2;
qam_symbols = qammod(Y.', 4, 'InputType', 'bit', 'UnitAveragePower', true);

snr_db = 20;
snr_linear = 10^(snr_db/10);
noise_var = 1/snr_linear;
noise = sqrt(noise_var/2) * randn(size(qam_symbols));

rx_symbols = qam_symbols + noise;

demod_llrs = qamdemod(rx_symbols, 4, ...
    'OutputType', 'approxllr', ...
    'NoiseVariance', noise_var, ...
    'UnitAveragePower', true);


data = [cmd_code, 10, length(demod_llrs), demod_llrs.'];

Z = NR_LDPC_client(server_ip, server_port, data, K);

disp("Decoded: " + num2str(Z));
