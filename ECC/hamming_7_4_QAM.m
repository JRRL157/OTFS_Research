clc;
clear;

N = 4;
M = 128;

m = floor(log2(N*M + 1));

disp(m);
[H, G, n, k] = hammgen(m);

disp(G);
disp(H);

m = randi([0 1], [k 1]);

encoded = mod(G.' * m, 2);

disp('Encoded: ');
disp(num2str(encoded));

% Adding one error
i = 100;
encoded(i) = mod(encoded(i) + 1, 2);

syndrome = mod(H * encoded, 2);

disp('Syndrome: ');
disp(syndrome);
disp(H);

% Find the error position by comparing syndrome with columns of H
error_pos = find(ismember(H.', syndrome.', 'rows'));

if ~isempty(error_pos)
    % Correct the error
    encoded(error_pos) = mod(encoded(error_pos) + 1, 2);
end

disp('Decoded: ');
decoded = encoded(end-k+1:end);
disp(num2str(decoded));

if isequal(m, decoded)
    disp('Decoded message matches the original message.');
else
    disp('Decoded message does NOT match the original message.');
end
