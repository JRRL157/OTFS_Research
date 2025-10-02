#!/bin/bash
EXECUTABLE_FILE=example_16QAM.bin
g++ -o $EXECUTABLE_FILE example_16QAM.cpp LDPC.cpp nrLDPCTables.cpp