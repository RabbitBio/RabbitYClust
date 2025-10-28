#pragma once
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
//#include <filesystem>
#include <atomic>
#include <shared_mutex>
#include "KHFMinHash.h"
#include "SharedData.h"
class ProteinProcessor {
	public:
		struct Config {
			int k = 8;
			int m = 15;
			int min_len = 10;
			bool use_xxhash = true;
			int num_threads = 1;
			bool output_binary = true;
		};

		// 构造函数
		explicit ProteinProcessor(const Config& cfg);

		// 构建sketch
		int build_sketches(const std::string input_fa,
				const std::string output_sketch,
				ProteinSketchData& protein_sketch_data
				);

		// 读取FA
		int load_sequences(const std::string input_fa,
				ProteinData& proteindata
				);

	private:
		// 线程局部数据
		struct ThreadLocal {
			std::vector<uint64_t> seq_ids;
			std::vector<std::string> names;
			std::vector<std::vector<uint64_t>> hashes;  // [m][...]
		};

		std::mutex file_read_mtx;

		// 线程任务
		void worker_thread(int thread_id,
				gzFile file,
				kseq_t* kseq,
				std::atomic<bool>& done);

		// 汇总线程结果
		void merge_thread_results(
			std::vector<uint64_t>& seq_ids_,
			std::vector<std::vector<uint64_t>>& hashes_
		);

		// 写二进制文件
		void write_binary(const std::string path, std::vector<std::vector<uint64_t>>& hashes) const;

		// 重新排序（多线程乱序 → 按 seq_id 排序）
		void reorder_hashes(
			std::vector<uint64_t>& seq_ids_,
			std::vector<std::vector<uint64_t>>& hashes_
		);

		// 配置
		Config config_;

		// 线程局部结果
		std::vector<ThreadLocal> thread_locals_;

		// 线程安全的结果容器
		mutable std::shared_mutex result_mtx_;
		//std::vector<uint64_t> seq_ids_;
		//std::vector<std::string> names_;
		//std::vector<std::vector<uint64_t>> hashes_;  // hashes_[m][n_seqs]
		//std::unordered_map<uint64_t, std::string> seq_map_;

		// 线程同步
		std::atomic<int> next_seq_id_{0};
};
