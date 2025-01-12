#ifndef CLUSTER_H
#define CLUSTER_H
#include "cdhit-common.h"
#include "input_sequence.h"

extern Options options;
class cluster{
public:
	vector<int> parent;
	cluster() {
		options = Options();
	}

	void cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent);
};
#endif
