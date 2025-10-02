/**
 * @file matlab_wrapper.h
 * @brief Header file for MATLAB wrapper functions for LDPC encoding and decoding.
 *
 * This file contains the declarations of functions that serve as a bridge between MATLAB and C++ implementations
 * of LDPC (Low-Density Parity-Check) coding. The functions allow for constructing an LDPC object, encoding data,
 * and decoding soft data using LDPC coding.
 */
#include "LDPC.h"
#include <math.h>
#include <iostream>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <cstring>
#include <unistd.h>

#define PORT 5016
#define MAX_MATRIX_DIMENSION 128
#define MAX_BUFFER_SIZE ((MAX_MATRIX_DIMENSION * MAX_MATRIX_DIMENSION * sizeof(double) * 2) + sizeof(double) * 3)

/**
 * @brief Instantiates or reconfigures the LDPC object.
 *
 * This function takes a pointer to the input data and initializes or reconfigures the LDPC object
 * based on the provided parameters.
 *
 * @param data_ptr Pointer to the input data array containing configuration parameters (as double).
 * @param response Pointer to the output buffer where the response (as uint8_t) will be stored.
 */
void constructorLDPC(double* data_ptr, uint8_t* response);

/**
 * @brief Encodes input data using LDPC (Low-Density Parity-Check) coding.
 *
 * This function takes a pointer to the input data and produces the LDPC encoded output.
 *
 * @param data_ptr Pointer to the input data array to be encoded (as double).
 * @param encoded_data Pointer to the output buffer where the encoded data (as uint8_t) will be stored.
 */
void encodeLDPC(double* data_ptr, uint8_t* encoded_data);

/**
 * @brief Decodes soft data using LDPC (Low-Density Parity-Check) coding.
 *
 * This function takes a pointer to the input data (real number) and produces the LDPC decoded output.
 *
 * @param data_ptr Pointer to the input data array to be decoded (soft decoding).
 * @param encoded_data Pointer to the output buffer where the decoded data (as uint8_t) will be stored.
 */
void decodeLDPC(double* data_ptr, uint8_t* decoded_data);