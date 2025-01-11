#ifndef INPUT_SEQUENCE_H
#define INPUT_SEQUENCE_H

#include <iostream>
#include <cstring>
#include <string>

struct Input_Sequence
{
	int seq_id;
    char* identifier;
    char* data;
};

struct Sequence_new
{
	int seq_id;
    char* data;
	Sequence_new() {}
//    char* identifier = null;
};
#endif
