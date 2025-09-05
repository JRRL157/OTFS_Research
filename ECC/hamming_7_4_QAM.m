clc;
clear;

N = 16;
M = 16;
mod_order = 256;
bps = log2(mod_order); %Bits per symbol

% OTFS Frame Length (Decimal)
L = N * M;

% OTFS Frame Length (Binary)
bL = L * bps;

% Codeword length (n) is 2^m -1, n = 2^m -1
% Message length (k) is n - m, k = n - m
m = floor(log2(bL + 1));

% disp(m);
[H, G, n, k] = hammgen(m);
codeRate = k / n;

msg = randi([0 1], [k 1]);
encoded = mod(G.' * msg, 2);
encoded_padded = [1 encoded.'].';

tx_symbols = qammod(encoded_padded, mod_order, InputType='bit');

% Adding Noise
UncodedEbNo = 6;
SNR = convertSNR(UncodedEbNo,"ebno","SNR", ...
    BitsPerSymbol=bps, ...
    CodingRate=codeRate);

received_signal = awgn(tx_symbols,SNR);

rx_symbols = qamdemod(received_signal, mod_order, OutputType='bit'); 
demod_unpadded = rx_symbols(2: end);
syndrome = mod(H * demod_unpadded, 2);

disp('Syndrome: ');
disp(syndrome);
% disp(H);

% Find the error position by comparing syndrome with columns of H
error_pos = find(ismember(H.', syndrome.', 'rows'));

if ~isempty(error_pos)
    % Correct the error
    demod_unpadded(error_pos) = mod(demod_unpadded(error_pos) + 1, 2);
end

disp('Decoded: ');
decoded = demod_unpadded(end-k+1:end);
disp(num2str(decoded.'));

if isequal(msg, decoded)
    disp('Decoded message matches the original message.');
else
    disp('Decoded message does NOT match the original message.');
end
