#ifndef CONFIG_DATA_H
#define CONFIG_DATA_H
#include <string>
#include <vector>

struct Data {
	std::vector<uint64_t> value;
	int id;
	Data(int id, std::vector<uint64_t> hash) : id(id), value(hash) {}
	Data() {}
};

struct sharedData {
	uint32_t seq_id = 0;
	uint32_t group_id = 0;
	uint64_t value = 0;

	sharedData(uint32_t s_id, uint32_t g_id, uint64_t val) : seq_id(s_id), group_id(g_id), value(val) {}

	sharedData() {}
};

struct GroupNode {
	int id;
	int root;
	GroupNode(int id, int root_id) : id(id), root(root_id) {}
	GroupNode() {}
};

#endif
