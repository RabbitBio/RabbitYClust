#include <zlib.h>
#include "CLI11.hpp"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
	CLI::App app{"read_fai v.0.0.1, read fasta index of protein files"};
	string filename = "";
	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	option_input->required();
	CLI11_PARSE(app, argc, argv);
	
	return 0;
}
