#pragma once
#include <vector>
#include <unordered_map>
#include <string>

struct ProteinSketchData {
	struct Config {
		int k = 8;
		int sketch_size = 15;
		int min_len = 10;
		bool use_xxhash = true;
	};
	Config config_;

	explicit ProteinSketchData(const Config& cfg) : config_(cfg) {
		hashes.assign(config_.sketch_size, {});
	}
	std::vector<std::vector<uint64_t>> hashes;  // hashes[i] = 第 i 列哈希
	//void reserve_items(int n) {
	//	cfg.items = n;
	//	names.reserve(n);
	//	sequence_map.reserve(n);
	//}

	//void reserve_sketch_size(int sketch_size) {
	//	if(sketch_size > 0) {
	//		cfg.sketch_size = sketch_size;
	//		has
	//	}
	//}
};

struct ProteinData {
	std::vector<std::string>  names;
	std::unordered_map<uint64_t, std::string> sequence_map; //序列id -> 序列内容 
};
//class SketchResult {
//public:
//	void append(uint64_t seq_id,
//			std::string name,
//			std::vector<uint64_t> hash_row);  // hash_row 长度 = m
//
//	void append_batch(const std::vector<uint64_t>& ids,
//			const std::vector<std::string>& names,
//			const std::vector<std::vector<uint64_t>>& hashes);
//
//	// === 只读接口（返回副本或视图）===
//	std::vector<uint64_t>          get_seq_ids() const;
//	std::vector<std::string>       get_names()   const;
//	std::vector<std::vector<uint64_t>> get_hashes() const;
//
//	// C++20: 零拷贝只读视图
//	// std::span<const uint64_t> seq_ids_view() const;
//	// std::span<const std::string> names_view() const;
//
//	size_t size() const;
//	void clear();
//	bool empty() const;
//
//	void reserve(size_t n, size_t m);  // 预分配
//	void set_hash_columns(size_t m);   // 设置 hashes 列数
//
//private:
//	mutable std::shared_mutex mtx_;
//	std::vector<uint64_t>          seq_ids_;
//	std::vector<std::string>       names_;
//	std::vector<std::vector<uint64_t>> hashes_;  // hashes_[i] = 第 i 列哈希
//};

//class SharedData {
//public:
//	SharedData() {}
//
//	std::vector<std::string>& getNames() { return names; }
//	std::unordered_map<uint64_t, std::string>& getAA() { return aa_map; }
//	std::vector<std::vector<uint64_t>>& getSketches() { return all_sketches; }
//
//private:
//	//SharedData(const SharedData&) = delete; // 禁止拷贝
//	//SharedData& operator=(const SharedData&) = delete;
//
//	std::vector<std::vector<uint64_t>> all_sketches;
//	std::vector<std::string> names;
//	std::unordered_map<uint64_t, std::string> aa_map;
//};
