import numpy as np
import sk_dsp_comm.digitalcom as dc
import matplotlib.pyplot as plt

# Parameters
M = 16  # QAM modulation order (16-QAM)
k = int(np.log2(M))  # Bits per symbol
num_symbols = 16  # Number of symbols to transmit
num_iterations = 1000  # Number of iterations
snr_db_range = np.arange(0, 30, 2)

# OFDM Parameters
num_subcarriers = 64  # Number of OFDM subcarriers
cyclic_prefix_length = 16

# OTFS Parameters
delay_bins = 16
doppler_bins = 64


def add_awgn(signal, snr_db):
    """ Add AWGN noise to a signal """
    snr_linear = 10**(snr_db / 10)
    noise_power = np.mean(np.abs(signal)**2) / snr_linear
    noise = np.sqrt(noise_power / 2) * (np.random.randn(*signal.shape) + 1j * np.random.randn(*signal.shape))
    return signal + noise


def ofdm_modulate(qam_symbols):
    """ OFDM Modulation with QAM symbols """
    # Reshape QAM symbols into OFDM subcarriers
    qam_symbols = np.pad(qam_symbols, (0, num_subcarriers - len(qam_symbols) % num_subcarriers), mode='constant')
    qam_symbols = qam_symbols.reshape(-1, num_subcarriers)
    # Apply IFFT and add cyclic prefix
    ofdm_symbols = np.fft.ifft(qam_symbols, axis=1)
    ofdm_symbols_cp = np.hstack([ofdm_symbols[:, -cyclic_prefix_length:], ofdm_symbols])
    return ofdm_symbols_cp.flatten()


def ofdm_demodulate(received):
    """ OFDM Demodulation with QAM symbols """
    # Remove cyclic prefix
    received = received.reshape(-1, num_subcarriers + cyclic_prefix_length)
    received = received[:, cyclic_prefix_length:]
    # Apply FFT
    return np.fft.fft(received, axis=1).flatten()


def sfft(x):
    """ Symplectic Finite Fourier Transform (SFFT) """
    return np.fft.fft(np.fft.fft(x, axis=0), axis=1)


def isfft(X):
    """ Inverse Symplectic Finite Fourier Transform (ISFFT) """
    return np.fft.ifft(np.fft.ifft(X, axis=1), axis=0)


def otfs_modulate(qam_symbols):
    """ OTFS Modulation with QAM symbols """
    # Reshape QAM symbols into delay-Doppler grid
    qam_symbols = np.pad(qam_symbols, (0, delay_bins * doppler_bins - len(qam_symbols)), mode='constant')
    qam_symbols = qam_symbols.reshape(delay_bins, doppler_bins)
    # Apply SFFT
    return sfft(qam_symbols).flatten()


def otfs_demodulate(received):
    """ OTFS Demodulation with QAM symbols """
    received = received[:delay_bins * doppler_bins]
    received = received.reshape(delay_bins, doppler_bins)
    # Apply ISFFT
    return isfft(received).flatten()


# BER Simulation
schemes = ['OFDM', 'OTFS']
results = {scheme: [] for scheme in schemes}

for snr_db in snr_db_range:
    for scheme in schemes:
        bit_errors = 0
        num_bits = 0

        for i in range(num_iterations):
            # Generate QAM symbols
            _, _, tx_data = dc.qam_gray_encode_bb(n_symb=num_symbols, ns=k, mod=M, pulse='rect', alpha=1.0)

            # Apply multicarrier modulation
            if scheme == 'OFDM':
                modulated = ofdm_modulate(tx_data)
            elif scheme == 'OTFS':
                modulated = otfs_modulate(tx_data)

            # Add AWGN
            received = add_awgn(modulated, snr_db)

            # Demodulation
            if scheme == 'OFDM':
                demodulated = ofdm_demodulate(received)
            elif scheme == 'OTFS':
                demodulated = otfs_demodulate(received)

            # Decode QAM symbols
            rx_data = dc.qam_gray_decode(demodulated, mod=M)

            # Count bit errors
            bit_errors += np.sum(rx_data != tx_data)
            num_bits += len(tx_data)

        # Calculate BER
        ber = bit_errors / num_bits
        results[scheme].append(ber)
        print(f"Scheme: {scheme}, SNR: {snr_db} dB, BER: {ber}")

# Plot BER vs SNR
plt.figure(figsize=(10, 6))
for scheme in schemes:
    plt.semilogy(snr_db_range, results[scheme], marker='o', label=scheme)

plt.title("BER vs SNR for OFDM and OTFS with QAM")
plt.xlabel("SNR (dB)")
plt.ylabel("Bit Error Rate (BER)")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()