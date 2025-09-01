clc;
clear;

m = [1 0 0 1];
[H, G, n, k] = hammgen(3);

disp(G);
disp(H);

encoded = mod(m * G, 2);

disp('Encoded: ');
disp(num2str(encoded));

% Adding one error
i = 5;
encoded(i) = mod(encoded(i) + 1, 2);
encoded(i+1) = mod(encoded(i+1) + 1, 2);

syndrome = mod(H * encoded.', 2);

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
decoded = encoded(4:7);
disp(num2str(decoded));