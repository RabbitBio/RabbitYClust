#include "unionfind.h"

void UnionFind::sortSetsByRoot(vector<pair<int, int>>& Id_Root_set) {
	auto construct_start = chrono::high_resolution_clock::now();
	Id_Root_set.reserve(parent.size());
	findRoot(Id_Root_set);
	auto construct_end = chrono::high_resolution_clock::now();
	auto construct_duration = chrono::duration_cast<chrono::seconds>(construct_end - construct_start).count();
	cout << "construct sort_set time needed is " << construct_duration << endl;

	auto sortset_start = chrono::high_resolution_clock::now();
	sort(Id_Root_set.begin(), Id_Root_set.end(), [](const pair<int, int>& a, const pair<int, int>& b){
		return a.second < b.second;
	});	
	auto sortset_end = chrono::high_resolution_clock::now();
	auto sortset_duration = chrono::duration_cast<chrono::seconds>(sortset_end - sortset_start).count();
	cout << "sort for sort_set time needed is " << sortset_duration << endl;
}

unordered_map<int, vector<int>> UnionFind::getGroupMap() {
	vector<pair<int, int>> sort_sets;
	sortSetsByRoot(sort_sets);
	auto map_start = chrono::high_resolution_clock::now();
	unordered_map<int, vector<int>> groups;
	int countSize = 0;
    for (const auto& p : sort_sets) {
    	groups[p.second].push_back(p.first); 
    }

	cout << "size of hashmap is " << groups.size() << endl;
	auto map_end = chrono::high_resolution_clock::now();
	auto map_duration = chrono::duration_cast<chrono::seconds>(map_end - map_start).count();
	cout << "time need for hashmap is " << map_duration << endl;
	return groups;
}


