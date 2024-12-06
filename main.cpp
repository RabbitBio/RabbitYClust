#include "include/cdhit-common.h"
#include "include/input_sequence.h"
#include "include/cdhit-lib.h"
#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>

bool read(std::vector<Input_Sequence*>& vec) {
    const char* filename = "/home/user_home/maguiliang/ycluster/RabbitYClust/python/result/bothtestgroups/AWA27912.1.fa";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file!" << std::endl;
        return false;
    }
    std::string line;  // 用于存储每行内容
    while (std::getline(file, line)) {
        Input_Sequence* seq = new Input_Sequence;
        seq->identifier = new char[line.size() + 1];
        std::strcpy(seq->identifier, line.c_str());

        if (!std::getline(file, line)) {
            std::cerr << "Data Error!" << std::endl;
            return false;
        }
        seq->data = new char[line.size() + 1];
        std::strcpy(seq->data, line.c_str());

        vec.push_back(seq);
    }
    file.close();  // 关闭文件
    return true;
}

int main(int argc, char* argv[]) {
    std::vector<Input_Sequence*> seq;
    read(seq);
    // std::cout << "<<<<<<<<<<" << seq.size() << std::endl;
    std::cout << std::endl;
    std::string input = "/home/user_home/maguiliang/codes/c_cpp/cdhit_test/test/input/test.fa";
    std::string output_file = "/home/user_home/maguiliang/codes/c_cpp/cdhit_test/test/output/test.fa";
    std::string output_vector = "/home/user_home/maguiliang/codes/c_cpp/cdhit_test/test/output/test.fa";
    cdhit_file(input, output_file, "-t 0.9 -M 160000 -T 1");
    std::cout << "<<<<<<<" << std::endl;
    // cdhit_vector(seq, output_vector, "-t 0.8 -M 16000 -T 4");
}

// g++ -I./include -L./libs /home/user_home/maguiliang/codes/c_cpp/cdhit_test/libs -fopenmp -lz -o main main.cpp
//  ./cd-hit -i /home/user_home/maguiliang/codes/c_cpp/cdhit_test/test/input/test.fa -o /home/user_home/maguiliang/cdhit/result/test_file.fa -t 0.9 -M 1600000 -T 1