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

struct GroupNode {
	int id;
	int root;
	GroupNode(int id, int root_id) : id(id), root(root_id) {}
	GroupNode() {}
};

#endif
