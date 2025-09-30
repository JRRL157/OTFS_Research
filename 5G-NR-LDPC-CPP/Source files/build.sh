#!/bin/bash
EXECUTABLE_FILE=example_run.bin
g++ -o $EXECUTABLE_FILE example_run.cpp LDPC.cpp nrLDPCTables.cpp