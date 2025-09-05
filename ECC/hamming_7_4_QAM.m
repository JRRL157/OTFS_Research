clc;
clear;

N = 4;
M = 128;
mod_order = 256;

% OTFS Frame Length (Decimal)
L = N * M;

% OTFS Frame Length (Binary)
bL = L * log2(mod_order); 

% Codeword length (n) is 2^m -1, n = 2^m -1
% Message length (k) is n - m, k = n - m
m = floor(log2(bL + 1));

%disp(m);
[H, G, n, k] = hammgen(m);

%disp(G);
%disp(H);

msg = randi([0 1], [k 1]);

disp('Message: ')
disp(num2str(msg.'));

encoded = mod(G.' * msg, 2);

%disp('Encoded: ');
%disp(num2str(encoded));

% Adding one error
i = 100;
encoded(i) = mod(encoded(i) + 1, 2);

syndrome = mod(H * encoded, 2);

%disp('Syndrome: ');
%disp(syndrome);
%disp(H);

% Find the error position by comparing syndrome with columns of H
error_pos = find(ismember(H.', syndrome.', 'rows'));

if ~isempty(error_pos)
    % Correct the error
    encoded(error_pos) = mod(encoded(error_pos) + 1, 2);
end

disp('Decoded: ');
decoded = encoded(end-k+1:end);
disp(num2str(decoded.'));

if isequal(msg, decoded)
    disp('Decoded message matches the original message.');
else
    disp('Decoded message does NOT match the original message.');
end
