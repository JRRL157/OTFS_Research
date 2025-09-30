#!/bin/bash
EXECUTABLE_FILE=example_4QAM.bin
g++ -o $EXECUTABLE_FILE example_4QAM.cpp LDPC.cpp nrLDPCTables.cpp