#ifndef CDHIT_LIB_H
#define CDHIT_LIB_H

#include "cdhit-common.h"
#include "input_sequence.h"
#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>

extern Options options;
extern SequenceDB seq_db;

void cdhit_vector(std::vector<Input_Sequence*>& seq, std::string output, std::string command);

void cdhit_file(std::string input, std::string output, std::string command);

#endif
