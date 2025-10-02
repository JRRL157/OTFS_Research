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
        n_bytes = recvfrom(sockfd, buffer, MAX_BUFFER_SIZE, 0, (struct sockaddr *)&client_addr, &client_addr_len);

        if (n_bytes < 0) {
            perror("Receive failed");
            continue;
        }

        printf("Received %d bytes from %s:%d\n", n_bytes, inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));

        double *received_data = (double *)buffer;
        uint16_t op_code = (uint16_t)received_data[0];
        double* data_ptr = &received_data[1];
        
        printf("Operation code: %d INCOMING!!\n", op_code);

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

        if (sendto(sockfd, response_buffer, 64, 0, (const struct sockaddr *)&client_addr, client_addr_len) < 0) {
            perror("Send failed");
        } else {
            printf("Sent %d bytes back to %s:%d\n", 64, inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));
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
    uint32_t data_len = static_cast<uint32_t>(data_ptr[0]);
    std::size_t M = static_cast<std::size_t>(data_ptr[1]);

    std::vector<bool> msg(data_len);
    for (uint32_t i = 0; i < data_len; i++) {
        msg[i] = static_cast<bool>(data_ptr[1 + i]);
    }

    // Filler bits
    std::vector<bool> fillers(ldpc.getFillerLength(), 0);
    std::vector<bool> extMsg = msg;
    extMsg.insert(extMsg.end(), fillers.begin(), fillers.end());

    // LDPC encoding
    std::vector<bool> enc = ldpc.encode(extMsg);
    assert(ldpc.checkSumCodeWord(enc));

    //rate matching
    std::vector<bool> rm_enc = ldpc.rateMatch(enc, M);

    // Copy encoded data to output buffer
    for(std::size_t i = 0; i < rm_enc.size(); i++) {
        encoded_data[i] = static_cast<uint8_t>(rm_enc[i]);
    }
}

void decodeLDPC(double* data_ptr, uint8_t* decoded_data) {
    unsigned int nMaxIter = static_cast<int>(data_ptr[0]);
    uint32_t data_len = static_cast<uint32_t>(data_ptr[1]);
    
    std::vector<double> llr(data_len);
    for (uint32_t i = 0; i < data_len; i++) {
        llr[i] = data_ptr[2 + i];
    }
    
    // rate recovery
    std::vector<double> rr_llr = ldpc.rateRecover(llr);

    // scl decoding
    std::vector<bool> msg_cap = ldpc.decode(rr_llr, nMaxIter);

    // Copy decoded data to output buffer
    for(std::size_t i = 0; i < msg_cap.size(); i++) {
        decoded_data[i] = static_cast<uint8_t>(msg_cap[i]);
    }
}