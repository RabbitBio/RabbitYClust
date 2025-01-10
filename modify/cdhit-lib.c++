#include "cdhit-lib.h"

Options options;
SequenceDB seq_db;

char** split_string_to_char_array(const std::string& input, int& count) {
    std::vector<std::string> tokens;
    std::istringstream iss(input);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    count = tokens.size();
    char** result = new char* [count];
    for (int i = 0;i < count;i++) {
        result[i] = new char[tokens[i].length() + 1];
        strcpy(result[i], tokens[i].c_str());
    }
    return result;
}

void free_char_array(char** arr, int count) {
    for (int i = 0;i < count;i++) {
        delete[] arr[i];
    }
    delete[] arr;
}

void cdhit_vector(std::vector<Input_Sequence*>& seq, std::string output, std::string command) {

    options = Options();
    seq_db = SequenceDB();

    int argc = 0;
    char** argv = split_string_to_char_array("cdhit-lib -o " + output + " " + command, argc);

    string db_out;

    if (seq.size() == 0) {
        std::cout<<"Empty input sequence"<<std::endl;
        return;
    }

    float begin_time = current_time();
    float end_time;

    // ***********************************    parse command line and open file
    if (argc < 5) print_usage(argv[0]);
    if (options.SetOptions(argc, argv) == 0)
        print_usage(argv[0]);
    options.SetOptionsInputVector();
    options.Validate();

    db_out = options.output;
    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    seq_db.Readvector(seq, options);

    cout << "total seq: " << seq_db.sequences.size() << endl;
    seq_db.SortDivide(options);

    seq_db.DoClustering(options);

    printf("writing new database\n");

    seq_db.WriteClustersvector(seq, db_out.c_str(), options);

    // write a backup clstr file in case next step crashes
    seq_db.WriteExtra1D(options);
    cout << "program completed !" << endl << endl;
    end_time = current_time();
    printf("Total CPU time %.2f\n", end_time - begin_time);

    free_char_array(argv, argc);

    return;
}

void cdhit_file(std::string input, std::string output, std::string command) {

    options = Options();
    seq_db = SequenceDB();

    int argc = 0;
    char** argv = split_string_to_char_array("cdhit-lib -i " + input + " -o " + output + " " + command, argc);

    string db_in;
    string db_out;

    float begin_time = current_time();
    float end_time;

    // ***********************************    parse command line and open file
    if (argc < 5) print_usage(argv[0]);
    if (options.SetOptions(argc, argv) == 0)
        print_usage(argv[0]);
    options.Validate();

    db_in = options.input;
    db_out = options.output;

    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    //printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );

    seq_db.Read(db_in.c_str(), options);

    cout << "total seq: " << seq_db.sequences.size() << endl;
    seq_db.SortDivide(options);

    seq_db.DoClustering(options);

    printf("writing new database\n");

    seq_db.WriteClusters(db_in.c_str(), db_out.c_str(), options);

    // write a backup clstr file in case next step crashes
    seq_db.WriteExtra1D(options);
    cout << "program completed !" << endl << endl;
    end_time = current_time();
    printf("Total CPU time %.2f\n", end_time - begin_time);

    free_char_array(argv, argc);

    return;
}

// int main(int argc, char* argv[]) {
//     std::vector<Input_Sequence*> seq;
//     read(seq);
//     std::string input = "/home/user_home/maguiliang/ycluster/RabbitYClust/python/result/bothtestgroups/AAF93073.1.fa";
//     std::string output = "/home/user_home/maguiliang/cdhit/result/test_file2.fa";
//     cdhit_file(input, output, "-t 0.8 -M 16000 -T 8");

//     // cdhit_vector(seq, "-o /home/user_home/maguiliang/cdhit/result/test_vector2.fa -t 0.8 -M 16000 -T 8");
// }