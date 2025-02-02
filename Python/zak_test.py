import numpy as np
import random

def zak_transform(x, N):
    """
    Compute the discrete Zak transform (DZT) of a signal x[n].
    
    Parameters:
    x : ndarray
        Input signal (1D array).
    N : int
        Period of the signal for modular mapping.
        
    Returns:
    Z : ndarray
        Zak transform (2D array of shape (N, len(x) // N)).
    """
    L = len(x) // N
    Z = np.zeros((N, L), dtype=complex)
    for k in range(N):
        for m in range(L):
            for n in range(L):
                Z[k, m] += x[n * N + k] * np.exp(-1j * 2 * np.pi * n * m / L)
    return Z

def inverse_zak_transform(Z, N):
    """
    Compute the inverse Zak transform (IDZT) to recover the signal x[n].
    
    Parameters:
    Z : ndarray
        Zak transform (2D array of shape (N, L)).
    N : int
        Period of the signal for modular mapping.
        
    Returns:
    x : ndarray
        Reconstructed signal (1D array).
    """
    L = Z.shape[1]
    x = np.zeros(N * L, dtype=complex)
    for n in range(L):
        for k in range(N):
            for m in range(L):
                x[n * N + k] += Z[k, m] * np.exp(1j * 2 * np.pi * n * m / L)
    return np.real(x) / L

# Function to test with assertions
def test_zak_transform(L, signal_length, period):
    pass_cnt = 0

    for _ in range(L):
        # Generate random signal
        x = np.random.randint(1, 100, size=signal_length)

        # Compute Zak transform and its inverse
        Z = zak_transform(x, period)
        x_reconstructed = inverse_zak_transform(Z, period)

        # Select random indices to compare
        random_indices = random.sample(range(signal_length), min(10, signal_length))

        # Assert values
        for idx in random_indices:
            assert np.isclose(x[idx], x_reconstructed[idx], atol=1e-5), f"Mismatch at index {idx}: {x[idx]} != {x_reconstructed[idx]}"

        #print(f"Test passed for iteration with signal: {x}")
        pass_cnt+=1
    
    return pass_cnt

# Parameters
L = int(10e5)  # Number of iterations
signal_length = 16  # Length of the signal
period = 4  # Period

# Run tests
pass_times = test_zak_transform(L, signal_length, period)

print(f"Pass rate = {(pass_times/L)*100}%")