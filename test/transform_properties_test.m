clc; clear; close all;pkg load communications;

% Parameters
N = 16; % Matrix size (same as the number of subcarriers)

% Generate Transform Matrices
DZT_matrix = dzt_matrix(N); % Zak Transform Matrix
DFT_matrix = dft_matrix(N,N);
DWHT_matrix = dwht_matrix(N,N);
DCT_matrix = dct_matrix(N,N);
DZT_matrix_opt = zak_transform_matrix(N,N);

% Store matrices in a struct for convenience
transform_matrices = struct('DZT', DZT_matrix,'DZT_OPT', DZT_matrix_opt, 'DFT', DFT_matrix, ...
                            'DCT', DCT_matrix, 'DWHT', DWHT_matrix);

% Loop over all matrices and analyze
for field = fieldnames(transform_matrices)'
    transform_name = field{1};
    transform_matrix = transform_matrices.(transform_name);

    % Eigenvalues
    eigenvalues = eig(transform_matrix);

    % Sparsity
    sparsity_ratio = nnz(transform_matrix) / numel(transform_matrix); % Non-zero elements ratio

    % Display results
    fprintf('Analysis for %s Matrix:\n', transform_name);
    fprintf('Eigenvalues:\n');
    disp(eigenvalues);
    fprintf('Sparsity Ratio: %.4f\n', sparsity_ratio);
    fprintf('-----------------------------\n');
end
