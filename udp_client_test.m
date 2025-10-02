clc;
clear;

server_ip = '127.0.0.1';
server_port = 5016;

data = [0, 1000, 0.2, ones(1, 10)];

Y = NR_LDPC_client(server_ip, server_port, data, length(data));

disp(Y);