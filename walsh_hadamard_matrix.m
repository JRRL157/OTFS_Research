function H = walsh_hadamard_matrix(N)
    if mod(log2(N), 1) ~= 0
        error('N must be a power of 2');
    end
    H = 1;
    for k = 1:log2(N)
        H = [H, H; H, -H];
    end
end
