clear;
clc;

N = 4;
M = 4;

mod_order = 16;
bps = log2(mod_order); %Bits per symbol

% OTFS Frame Length (Decimal)
L = N * M;

% OTFS Frame Length (Binary)
L_bits = L * bps;

m = bps;

% Symbol correction capability 
t = 6;

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

num_blocks = floor(L_bits / bps);
num_src_bits = num_blocks * k;

% Random message bits which will be sent
random_bits_msg = randi([0 1], num_src_bits, 1);

% Encoded bits
encoded_bits = rsEncoder(random_bits_msg);
% encoded_bits = [];
% for i = 1:num_blocks
%     message_block = random_bits_msg((i-1)* (k*m) + 1 : i * (k * m));
%     disp(message_block);
%     cword = rsEncoder(message_block);
%     encoded_bits = [encoded_bits; cword];
% end

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

decoded_bits = rsDecoder(demod_bits);
% decoded_bits = [];
% for i = 1:num_blocks
% 	rec_cword =  demod_bits((i-1)*k + 1: i*k);
% 	% syndrome = mod(H * rec_cword, 2);
% 	% error_idx = find(ismember(H.', syndrome.', 'rows'));
% 	% if ~isempty(error_idx)
% 	% 	rec_cword(error_idx) =  mod(rec_cword(error_idx) + 1, 2);
% 	% end
% 	% msg_block = rec_cword(end-k+1:end);
%   % msg_block = bchdec(rec_cword, n, k);
%   msg_block = rsDecoder(rec_cword);
% 	decoded_bits = [decoded_bits; msg_block];
% end

if isequal(random_bits_msg, decoded_bits)
    disp('[SUCCESS] Decoded message matches the original message.');
else
    disp('[FAIL] Decoded message does NOT match the original message.');
end
