#include "LDPC.h"
#include <iostream>
#include <random>
#include <chrono>    // For timing
#include <cstdio>    // For snprintf
#include <vector>
#include <cmath>     // For sqrt
#include <cassert>   // For assert
#include <algorithm> // For std::transform

using namespace std;

using std::chrono::duration_cast;
using std::chrono::steady_clock;
using chrono::milliseconds;

// Function prototypes
vector<double> lin_space(double start, double end, int num);
void example_run_LDPC();

int main() {
	std::cout << "Hello, 5G-NR! (4-QAM Simulation)" << std::endl;

	example_run_LDPC();

	return 0;
}

void example_run_LDPC()
{
	unsigned nMaxIter = 8; // number of decoder iterations
	int M = 1000; // code word length (number of coded bits)
	double rate = 1.0 / 3; // code rate

	unsigned K = unsigned(ceil(M * rate)); // information length

	// instantiates an nrLDPC object
	nrLDPC ldpc = nrLDPC(K, rate);

	// random engines
	default_random_engine random_engine(steady_clock::now().time_since_epoch().count()); // Seed the engine for different results each run
	bernoulli_distribution bern_dist;
	normal_distribution<double> norm_dist(0.0, 1.0);

	// Running parameters
	// NOTE: Adjusted the SNR range for 4-QAM, which typically requires more power than BPSK.
	vector<double> EsN0_dB = lin_space(-2, 2, 9);
	vector<double> N0(EsN0_dB.size(), 0);
	transform(EsN0_dB.begin(), EsN0_dB.end(), N0.begin(), [](const double& x) {return pow(10.0, -x / 10.0); });

	vector<double> ber(N0.size(), 0), bler(N0.size(), 0);
	vector<unsigned> n_bit_errs(N0.size(), 0), n_blk_errs(N0.size(), 0);

	unsigned n_max_blks = 10000;
	unsigned n_min_err_blks = 100; // Stop simulation for an SNR point after 100 block errors
	vector<bool> fillers(ldpc.getFillerLength(), 0);
    
    // Make sure M is a multiple of 2 for pairing bits in 4-QAM
    assert(M % 2 == 0 && "Codeword length M must be even for 4-QAM.");

	// loop each SNR
	for (unsigned i = 0; i < N0.size(); i++) {
		//print progress
		char str[150]; // Increased buffer size
		snprintf(str, sizeof(str), "\nNow running EsN0: %.2f dB [%d of %lu]", EsN0_dB[i], i + 1, N0.size());
		cout << str << endl;
		size_t print_len = 0;

		unsigned n_blks_done = 0;
		auto tStart = steady_clock::now(); // Use steady_clock for more reliable timing

		while ((n_blks_done < n_max_blks) && (n_blk_errs[i] < n_min_err_blks)) {
			// generate random bit stream
			vector<bool> msg;
			msg.reserve(K);
			for (unsigned j = 0; j < K; j++)
				msg.push_back(bern_dist(random_engine));

			// add filler bits
			vector<bool> extMsg = msg;
			extMsg.insert(extMsg.end(), fillers.begin(), fillers.end());

			// LDPC encoding
			vector<bool> enc = ldpc.encode(extMsg);
			assert(ldpc.checkSumCodeWord(enc));

			//rate matching
			vector<bool> rm_enc = ldpc.rateMatch(enc, M);

			// =======================================================================
			// MODIFICATION START: Replaced BPSK with 4-QAM
			// =======================================================================
			
			// This single loop now handles modulation, channel noise, and LLR calculation
			vector<double> llr; 
			llr.reserve(M);
			const double norm_factor = 1.0 / sqrt(2.0); // Normalization factor for Es=1

			// Process bits in pairs for 4-QAM
			for (size_t j = 0; j < M; j += 2) {
				// 1. Get a pair of bits
				bool bit1 = rm_enc[j];
				bool bit2 = rm_enc[j + 1];

				// 2. 4-QAM Modulation (Gray mapping)
				// bit1 -> I-component, bit2 -> Q-component
				// bit 0 maps to +1/sqrt(2), bit 1 maps to -1/sqrt(2)
				double s_I = (1.0 - 2.0 * bit1) * norm_factor;
				double s_Q = (1.0 - 2.0 * bit2) * norm_factor;

				// 3. AWGN Channel
				// Add complex Gaussian noise. Each component (I and Q) has variance N0/2.
				double n_I = sqrt(N0[i] / 2.0) * norm_dist(random_engine);
				double n_Q = sqrt(N0[i] / 2.0) * norm_dist(random_engine);

				// Received complex symbol
				double r_I = s_I + n_I;
				double r_Q = s_Q + n_Q;

				// 4. LLR Calculation (Demodulation)
				// For QPSK with Gray mapping, the LLRs for the two bits can be calculated
				// independently from the I and Q components of the received signal.
				// The LLR for a BPSK-like channel with amplitude A is LLR = 2*A*r / (N0/2).
				// Here, A = 1/sqrt(2).
				// LLR = 2 * (1/sqrt(2)) * r / (N0/2) = 4*r / (sqrt(2)*N0) = 2*sqrt(2)*r / N0
				double llr1 = (2.0 * sqrt(2.0) * r_I) / N0[i];
				double llr2 = (2.0 * sqrt(2.0) * r_Q) / N0[i];

				llr.push_back(llr1);
				llr.push_back(llr2);
			}

			// =======================================================================
			// MODIFICATION END
			// =======================================================================


			// rate recovery
			vector<double> rr_llr = ldpc.rateRecover(llr);

			// ldpc decoding
			vector<bool> msg_cap = ldpc.decode(rr_llr, nMaxIter);

			// count errors
			unsigned n_errs = 0;
			for (unsigned j = 0; j < K; j++) {
				if (msg[j] != msg_cap[j])
					n_errs++;
			}

			if (n_errs) {
				n_bit_errs[i] += n_errs;
				n_blk_errs[i]++;
			}

			n_blks_done += 1;

			ber[i] = n_bit_errs[i] * 1.0 / K / n_blks_done;
			bler[i] = n_blk_errs[i] * 1.0 / n_blks_done;
			
			// print progress for every 10 blocks
			if (n_blks_done % 10 == 0 || n_blks_done == 1) {
				auto tNow = steady_clock::now();
				double elapsed_secs = duration_cast<milliseconds>(tNow - tStart).count() / 1000.0;
				int new_len = snprintf(str, sizeof(str), "Elapsed time: %.1f s, # blocks: %d, # error blocks: %d, ber: %.5e, bler %.5f", elapsed_secs, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
				cout << std::string(print_len, '\b') << str << flush;
				print_len = new_len;
			}
		}

		// print progress when one SNR is finished
		auto tNow = steady_clock::now();
		double elapsed_secs = duration_cast<milliseconds>(tNow - tStart).count() / 1000.0;
		snprintf(str, sizeof(str), "Elapsed time: %.1f s, # blocks: %d, # error blocks: %d, ber: %.5e, bler %.5f", elapsed_secs, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
		cout << std::string(print_len, '\b') << str << flush;
	}

	// print simulation result
	cout << endl << endl;
	cout << "================== Simulation Results ==================" << endl;
	cout << "Modulation: " << "4-QAM" << endl;
	cout << "[M,R] = [ " << M << ", " << rate << "]" << endl;
	cout << "EsN0_dB = [ ";
	for (auto e : EsN0_dB)
		cout << e << " ";
	cout << "]" << endl;

	cout << "BER = [ ";
	for (auto e : ber)
		printf("%.3e ", e);
	cout << "]" << endl;

	cout << "BLER = [ ";
	for (auto e : bler)
		printf("%.3e ", e);
	cout << "]" << endl;
	cout << "========================================================" << endl;
}

vector<double> lin_space(double start, double end, int num) {
	// catch rarely, throw often
	assert(num >= 2 && "The third parameter must be a positive integer >= 2!");

	vector<double> pts;
    pts.reserve(num);
	// length of each segment
	double step = (end - start) / (num - 1);
	for (int i = 0; i < num; i++) {
		pts.push_back(start + i * step);
	}
	return pts;
}