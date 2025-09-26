clear;
clc;

N = 4;
M = 128;

mod_order = 4;
bps = log2(mod_order); %Bits per symbol

% OTFS Frame Length (Decimal)
L = N * M;

% OTFS Frame Length (Binary)
L_bits = L * bps;

% Codeword length (n) is 2^m -1, n = 2^m -1
% Message length (k) is n - m, k = n - m
m = 4;
n = 2^m -1;
k = 5;
codeRate = k/n;

% Generator Polynomial
gp = bchgenpoly(n,k);

% Encoder Object
enc = comm.BCHEncoder(n, k, gp);
% Decoder Object
dec = comm.BCHDecoder(n, k, gp);

% BCH code correction capability
t = bchnumerr(n,k);

num_blocks = floor(L_bits / n);
num_src_bits = num_blocks * k;

% Random message bits which will be sent
random_bits_msg = randi([0 1], num_src_bits, 1);

% Encoded bits
encoded_bits = [];
for i = 1:num_blocks
	message_block = random_bits_msg((i-1)*k + 1 : i*k);
  cword = enc(message_block);
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
  msg_block = dec(rec_cword);
	decoded_bits = [decoded_bits; msg_block];
end

if isequal(random_bits_msg, decoded_bits)
    disp('[SUCCESS] Decoded message matches the original message.');
else
    disp('[FAIL] Decoded message does NOT match the original message.');
end
