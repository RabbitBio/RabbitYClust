#include "unionfind.h"




class GroupStream {
public:
	UnionFind uf;
	int items;
	int M;
	int R = 1;
	int L = 1;
	bool slide = true;
	vector<Data> hash_vec;
//	vector<pair<int, int>> hash_vec;

	GroupStream(int n, int m, int r, int l) : uf(n), items(n), R(r), L(l), M(m) {
		hash_vec.resize(items);
		resize();
	}

	GroupStream(int n, int r, int l) : uf(n), items(n), R(r), L(l) {
		hash_vec.resize(items);
		resize();
//		cerr << "hash_vec size is: " << hash_vec.size() << endl;
//		cerr << "groupstream gerneration ends." << endl;
	}
	
	GroupStream(int n) : uf(n), items(n) { hash_vec.resize(items); }

	void setSlideOff(){ slide = false; }
	void setM(int m) { M = m; }
	void setR(int r) { R = r; }
	void setL(int l) { L = l; }
	void resize() {
		for(auto& data : hash_vec)
			data.value.resize(L * R);
	}
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
