using namespace std;

void* consumer(void* arg);

void reorder_swap(vector<uint64_t>& matrix, vector<uint64_t>& indices);

void compactSeqIds(int threads);

template <typename T>
void write_to_file(const std::vector<std::vector<T>>& data, const std::string& filename);

int buildSketches(string filename, string output_filename, int num_threads, int k, int m, bool xxhash_flag, int min_len, bool output_on);

int read_fa(string filename, int min_len);
