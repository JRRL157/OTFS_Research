/**
 * @file matlab_wrapper.cpp
 * @brief UDP server for MATLAB integration with 5G-NR LDPC encoding/decoding.
 *
 * This file implements a UDP server that listens for data from MATLAB,
 * processes LDPC-related operations based on received op codes, and
 * interacts with the nrLDPC class for encoding and decoding.
 *
 * Supported operations:
 *   - LDPC constructor (op_code 0)
 *   - (Reserved for future: encoding, decoding, etc.)
 *
 * Functions:
 *   - constructorLDPC: Instantiates or reconfigures the LDPC object.
 *   - encodeLDPC: Encodes data using LDPC (to be implemented).
 *   - decodeLDPC: Decodes data using LDPC (to be implemented).
 *
 * Usage:
 *   Compile and run this server. MATLAB can send UDP packets with the
 *   required format to interact with LDPC functionality.
 *
 * @author JRRL157
 * @date October 2025
 */
#include "matlab_wrapper.h"

// instantiates a POLAR object
nrLDPC ldpc = nrLDPC(100, 0.5);

int main() {

    int sockfd;
    struct sockaddr_in server_addr, client_addr;
    socklen_t client_addr_len = sizeof(client_addr);
    
    char buffer[MAX_BUFFER_SIZE];
    uint32_t n_bytes;
    uint8_t response_buffer[MAX_BUFFER_SIZE];
    
    sockfd = socket(AF_INET, SOCK_DGRAM, 0);

    if (sockfd < 0) {
        perror("Socket creation failed");
        exit(EXIT_FAILURE);
    }

    memset(&server_addr, 0, sizeof(server_addr));
    server_addr.sin_family = AF_INET;
    server_addr.sin_addr.s_addr = INADDR_ANY;
    server_addr.sin_port = htons(PORT);

    if (bind(sockfd, (const struct sockaddr *)&server_addr, sizeof(server_addr)) < 0) {
        perror("Bind failed");
        close(sockfd);
        exit(EXIT_FAILURE);
    }

    printf("MATLAB Wrapper listening on port %d for data...\n", PORT);
    while (true) {
        memset(buffer, 0, MAX_BUFFER_SIZE);
        memset(response_buffer, 0, MAX_BUFFER_SIZE);
        
        n_bytes = recvfrom(sockfd, buffer, MAX_BUFFER_SIZE, 0, (struct sockaddr *)&client_addr, &client_addr_len);

        if (n_bytes < 0) {
            perror("Receive failed");
            continue;
        }

        printf("Received %d bytes from %s:%d\n", n_bytes, inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));

        double *received_data = (double *)buffer;
        uint8_t op_code = static_cast<uint8_t>(received_data[0]);
        double* data_ptr = &received_data[1];
        
        printf("#### [%d] Operation code: %d INCOMING!! ####\n", op_code, op_code);

        switch(op_code) {
            case 0: 
                constructorLDPC(data_ptr, response_buffer);
                break;
            case 1:
                encodeLDPC(data_ptr, response_buffer);
                break;
            case 2:
                decodeLDPC(data_ptr, response_buffer);
                break;
            default:
                break;
        }

        if (sendto(sockfd, response_buffer, 512, 0, (const struct sockaddr *)&client_addr, client_addr_len) < 0) {
            perror("Send failed");
        } else {
            printf("Sent %zu bytes back to %s:%d\n", 512, inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));
        }
    }

    close(sockfd);
    return 0;
}

void constructorLDPC(double* data_ptr, uint8_t* response) {
    std::size_t M = static_cast<std::size_t>(data_ptr[0]);
    double CR = data_ptr[1];

    printf("Received constructor parameters: M=%zu, CR=%.2f\n", M, CR);

    unsigned infoLen = unsigned(ceil(M * CR));

    printf("LDPC constructor called with M=%zu, CR=%.2f, infoLen=%u\n", M, CR, infoLen);

    ldpc = nrLDPC(infoLen, CR);

    printf("LDPC object created with K=%zu, N=%zu, R=%.4f, BGn=%u, Zc=%u, F=%zu\n", ldpc.mK, ldpc.mN, ldpc.mR, ldpc.mBGn, ldpc.mZc, ldpc.mF);
    response[0] = static_cast<uint8_t>(1);
}

void encodeLDPC(double* data_ptr, uint8_t* encoded_data) {
    std::size_t K = static_cast<std::size_t>(data_ptr[0]);
    printf("Encoding request for K=%zu\n", K);

    std::vector<bool> msg(K);
    printf("msg = ");
    for (uint32_t i = 0; i < K; i++) {
        msg[i] = static_cast<bool>(data_ptr[1 + i]);
        printf("%.0f ", data_ptr[1 + i]);
    }
    printf("\n");

    printf("Encoding message of length %zu\n", msg.size());

    // Filler bits
    std::vector<bool> fillers(ldpc.getFillerLength(), 0);
    std::vector<bool> extMsg = msg;
    extMsg.insert(extMsg.end(), fillers.begin(), fillers.end());

    printf("Message with fillers length added: %zu\n", extMsg.size());

    // LDPC encoding
    std::vector<bool> enc = ldpc.encode(extMsg);
    assert(ldpc.checkSumCodeWord(enc));

    //rate matching
    std::size_t M = static_cast<std::size_t>(ceil(K / ldpc.mR));
    std::vector<bool> rm_enc = ldpc.rateMatch(enc, M);

    printf("Encoded and rate-matched message length: %zu\n", rm_enc.size());

    printf("encoded = ");
    // Copy encoded data to output buffer
    for(std::size_t i = 0; i < rm_enc.size(); i++) {
        encoded_data[i] = static_cast<uint8_t>(rm_enc[i]);
        printf("%d ", encoded_data[i]);
    }
    printf("\n");

    printf("Encoding completed and data copied to response buffer.\n");
}

void decodeLDPC(double* data_ptr, uint8_t* decoded_data) {
    unsigned int nMaxIter = static_cast<int>(data_ptr[0]);
    uint32_t data_len = static_cast<uint32_t>(data_ptr[1]);

    printf("Decoding request with nMaxIter=%u, data_len=%u\n", nMaxIter, data_len);

    std::vector<double> llr(data_len);
    printf("r = ");
    for (uint32_t i = 0; i < data_len; i++) {
        llr[i] = data_ptr[2 + i];
        printf("%.2f ", llr[i]);
    }
    printf("\n");

    // rate recovery
    std::vector<double> rr_llr = ldpc.rateRecover(llr);
    printf("Rate recovery completed. LLR length: %zu\n", rr_llr.size());

    // scl decoding
    std::vector<bool> msg_cap = ldpc.decode(rr_llr, nMaxIter);
    printf("Decoding completed. Decoded message length: %zu\n", msg_cap.size());

    // Copy decoded data to output buffer
    printf("decoded = ");
    for(std::size_t i = 0; i < msg_cap.size(); i++) {
        decoded_data[i] = static_cast<uint8_t>(msg_cap[i]);
        printf("%d ", decoded_data[i]);
    }
    printf("\n");
}