import numpy as np
import sk_dsp_comm.digitalcom as dc
import matplotlib.pyplot as plt

# Parameters
M = 16  # Modulation order (16-QAM)
k = int(np.log2(M))  # Bits per symbol
num_symbols = 1000  # Number of symbols to transmit
num_iterations = 100000 # Number of iterations
snr_db_range = np.arange(0, 30, 2)

ber_results = []

for snr_db in snr_db_range:    
    bit_errors = 0
    num_bits = 0
    
    for i in range(num_iterations):
        _, _, tx_data = dc.qam_gray_encode_bb(n_symb=num_symbols, ns=k, mod=M, pulse='rect', alpha=1.0)
        # Add noise for simulation purposes        
        snr_linear = 10**(snr_db / 10)
        noise_power = 1 / (2 * snr_linear)  # Adjust to simulate different SNRs
        noise = np.sqrt(noise_power / 2) * (np.random.randn(len(tx_data)) + 1j * np.random.randn(len(tx_data)))
        received_symbols = tx_data + noise

        # Demodulate the noisy received symbols
        rx_data = dc.qam_gray_decode(x_hat=received_symbols, mod=2)

        # Calculate bit errors
        bit_errors += np.sum(rx_data != tx_data)
        num_bits += len(rx_data)

    ber = bit_errors / num_bits
    ber_results.append(ber)
    print(f"SNR: {snr_db} dB, BER: {ber}")

# Plot BER vs SNR
plt.figure()
plt.semilogy(snr_db_range, ber_results, marker='o')
plt.title("BER vs SNR for 16-QAM")
plt.xlabel("SNR (dB)")
plt.ylabel("Bit Error Rate (BER)")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.show()