#include "unionfind.h"




class GroupStream {
public:
	UnionFind uf;
	int items;
	int M;
	int R;
	bool slide = true;
	vector<Data> hash_vec;
//	vector<pair<int, int>> hash_vec;

	GroupStream(int n, int m, int r) : uf(n), items(n), M(m), R(r) {
		hash_vec.resize(items);
//		cerr << "hash_vec size is: " << hash_vec.size() << endl;
//		cerr << "groupstream gerneration ends." << endl;
	}

	void setSlideOff(){ slide = false; }
	void setR(int r) { R = r; }
	void setM(int m) { M = m; }
	void Group(vector<vector<uint64_t>>& hashes, unordered_map<int, vector<int>>& group_map);
	// grouping sequences with m hash-functions

	void Unite(vector<Data> dataList, UnionFind& uf);
 	// use unionfind to unite group results by per column
	void Sort(vector<Data>& dataList);
	// sort hash-vec by value 1.used for fast unite; 2.used for constructing GroupResMap
	void GroupByCol(vector<Data>& hash_vec, UnionFind& uf);
	// grouping by column

	void fillHashVec(const vector<vector<uint64_t>>& vec, vector<Data>& hash_vec, int m);
	// construct a sorted struct Data(hash-vec) for a column of hash-funtions(vec)

	void countGroupSize(UnionFind& uf, int m);
	// top 10 largest groups in m column
	
	void getGroupMap(UnionFind& uf,unordered_map<int, vector<int>>& group_map);
	//unordered_map<int, vector<int>> getGroupMap();

};
