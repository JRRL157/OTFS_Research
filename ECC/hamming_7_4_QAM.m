clc;
clear;

N = 16;
M = 16;
mod_order = 256;
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

num_blocks = floor(L_bits / n);
num_src_bits = num_blocks * k;

% Random message bits which will be sent
random_bits_msg = randi([0 1], num_src_bits, 1);

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

tx_symbols = qammod(padded_encoded_bits, mod_order, InputType='bit');

% Adding Noise
UncodedEbNo = 6;
SNR = convertSNR(UncodedEbNo,"ebno","SNR", ...
    BitsPerSymbol=bps, ...
    CodingRate=codeRate);

rx_symbols = awgn(tx_symbols,SNR);

demod_bits = qamdemod(rx_symbols, mod_order, OutputType='bit'); 

decoded_bits = [];
for i = 1:num_blocks
	rec_cword =  demod_bits((i-1)*n + 1: i*n);
	syndrome = mod(H * rec_cword, 2);
	error_idx = find(ismember(H.', syndrome.', 'rows'));
	if ~isempty(error_idx)
		rec_cword(error_idx) =  mod(rec_cword(error_idx) + 1, 2);
	end
	msg_block = rec_cword(end-k+1:end);
	decoded_bits = [decoded_bits; msg_block];
end

if isequal(random_bits_msg, decoded_bits)
    disp('Decoded message matches the original message.');
else
    disp('Decoded message does NOT match the original message.');
end
