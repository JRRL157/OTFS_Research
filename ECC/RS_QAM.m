clear;
clc;

N = 4;
M = 32;

mod_order = 16;
bps = log2(mod_order); %Bits per symbol

% OTFS Frame Length (Decimal)
L = N * M;

% OTFS Frame Length (Binary)
L_bits = L * bps;

m = bps;

% Symbol correction capability 
t = 5;

% Codeword length
n = 2^m -1;

% Message length
k = n - 2 * t;

codeRate = k/n;

% Encoder Object
rsEncoder = comm.RSEncoder( ...
    BitInput=true, ...
    CodewordLength=n, ...
    MessageLength=k);

% Decoder Object
rsDecoder = comm.RSDecoder( ...
    BitInput=true, ...
    CodewordLength=n, ...
    MessageLength=k);

% Random message bits which will be sent
random_bits_msg = randi([0 1],  L * bps, 1);

tx_bits = [];
encoded_bits = [];
num_blocks = floor((L * bps) / (k * m));

for i=1:num_blocks
  msg = random_bits_msg(1 + (i-1)*(k*m): i*(k*m));
  enc_msg = rsEncoder(msg);
  tx_bits = [tx_bits; msg];
  encoded_bits = [encoded_bits; enc_msg];
end
% Encoded bits padding
bits2pad = L * bps - length(encoded_bits);
encoded_bits = [encoded_bits; zeros(bits2pad, 1)];

tx_symbols = qammod(encoded_bits, mod_order, InputType='bit');

% Noise Config
UncodedEbNo = 6;
SNR = convertSNR(UncodedEbNo,"ebno","SNR", ...
    BitsPerSymbol=bps, ...
    CodingRate=codeRate);

% Adding Noise
rx_symbols = awgn(tx_symbols,SNR);

% Demodulation
demod_bits = qamdemod(rx_symbols, mod_order, OutputType='bit'); 

% Unpad bits
unpadded_bits = demod_bits(1: L * bps - bits2pad);

% RS Decoding
rx_bits = rsDecoder(unpadded_bits);

if isequal(tx_bits, rx_bits)
    disp('[SUCCESS] Decoded message matches the original message.');
else
    disp('[FAIL] Decoded message does NOT match the original message.');
end
