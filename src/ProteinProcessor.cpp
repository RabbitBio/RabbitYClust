#include "ProteinProcessor.h"
//#include "KHFMinHash.h"
//#include "kseq.h"
//#include <zlib.h>

#include <iostream>
#include <fstream>
#include <thread>
#include <algorithm>
#include <numeric>

ProteinProcessor::ProteinProcessor(const Config& cfg) : config_(cfg) {
	if (config_.num_threads <= 0) {
		config_.num_threads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
	}
	thread_locals_.resize(config_.num_threads);
	for (auto& tl : thread_locals_) {
		tl.hashes.assign(config_.m, {});
	}
}

void ProteinProcessor::worker_thread(
		int thread_id,
		gzFile file,
		kseq_t* kseq,
		std::atomic<bool>& done) {
	auto& local = thread_locals_[thread_id];

	while (!done) {
		std::string seq, name;
		int len = -1;
		uint64_t seq_id;

		// 线程安全读取序列
		{
			std::lock_guard<std::mutex> lock(file_read_mtx);  // 保留全局 mtx1 给 kseq_read
			len = kseq_read(kseq);
			if (len < 0) {
				done = true;
				break;
			}
			if (len < config_.min_len) continue;

			seq = kseq->seq.s;
			name = kseq->name.s;

			seq_id = next_seq_id_.fetch_add(1);
		}

		// 计算 MinHash
		Sketch::KHFMinHash mh;
		mh.setK(config_.k);
		mh.setM(config_.m);
		if (config_.use_xxhash) {
			mh.buildSketch(seq.c_str());
		} else {
			mh.buildSketchByNoSeedAAHash(seq.c_str());
		}
		auto& sketch = mh.getSektch();

		// 分配全局 seq_id
		//uint64_t seq_id = next_seq_id_.fetch_add(1, std::memory_order_relaxed);

		local.seq_ids.push_back(seq_id);
		local.names.push_back(std::move(name));
		for (int i = 0; i < config_.m; ++i) {
			local.hashes[i].push_back(sketch.hashes[i]);
		}
	}
}

int ProteinProcessor::build_sketches(
		const std::string input_fa,
		const std::string output_sketch,
		ProteinSketchData& protein_sketch_data
		) {
	auto& hashes_ = protein_sketch_data.hashes;
	std::vector<uint64_t> seq_ids_;

	gzFile fp = gzopen(input_fa.c_str(), "r");
	if (!fp) {
		std::cerr << "Failed to open: " << input_fa << '\n';
		return 0;
	}
	kseq_t* ks = kseq_init(fp);

	next_seq_id_ = 0;
	for (auto& tl : thread_locals_) {
		tl.seq_ids.clear();
		tl.names.clear();
		for (auto& row : tl.hashes) row.clear();
	}

	std::atomic<bool> done{false};
	std::vector<std::thread> threads;
	threads.reserve(config_.num_threads);

	for (int i = 0; i < config_.num_threads; ++i) {
		threads.emplace_back(&ProteinProcessor::worker_thread, this, i, fp, ks, std::ref(done));
	}

	for (auto& t : threads) t.join();

	kseq_destroy(ks);
	gzclose(fp);

	std::cout << "Finish building sketch! Total: " << next_seq_id_ << " sequences.\n";

	merge_thread_results(seq_ids_, protein_sketch_data);
	if (config_.num_threads > 1) {
		reorder_hashes(seq_ids_, hashes_);
	}

	if (config_.output_binary && !output_sketch.empty()) {
		write_binary(output_sketch, protein_sketch_data);
	}

	return next_seq_id_;
}

void ProteinProcessor::merge_thread_results(
	std::vector<uint64_t>& seq_ids_,
	ProteinSketchData& protein_sketch_data
	//std::vector<std::vector<uint64_t>>& hashes_
	) {
	std::unique_lock lock(result_mtx_);

	seq_ids_.clear();
	ProteinSketchData::Config sketch_config{
		next_seq_id_,
		config_.m,
		config_.min_len
	};
	protein_sketch_data.reserveSketchesSize(sketch_config);
	auto& hashes_ = protein_sketch_data.hashes;
	//hashes_ = std::vector<std::vector<uint64_t>>(config_.m);

	for (const auto& tl : thread_locals_) {
		seq_ids_.insert(seq_ids_.end(),
				std::make_move_iterator(tl.seq_ids.begin()),
				std::make_move_iterator(tl.seq_ids.end()));
		for (int i = 0; i < config_.m; ++i) {
			hashes_[i].insert(hashes_[i].end(),
					std::make_move_iterator(tl.hashes[i].begin()),
					std::make_move_iterator(tl.hashes[i].end()));
		}
	}
}

void ProteinProcessor::reorder_hashes(
	std::vector<uint64_t>& seq_ids_,
	std::vector<std::vector<uint64_t>>& hashes_
) {
	if (seq_ids_.empty()) return;

	std::vector<size_t> order(seq_ids_.size());
	std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
			return seq_ids_[a] < seq_ids_[b];
			});

	// 重新排列 names 和 hashes
	std::vector<std::vector<uint64_t>> sorted_hashes(config_.m, std::vector<uint64_t>(seq_ids_.size()));

	for (size_t i = 0; i < seq_ids_.size(); ++i) {
		size_t old_idx = order[i];
		for (int j = 0; j < config_.m; ++j) {
			sorted_hashes[j][i] = hashes_[j][old_idx];
		}
	}

	hashes_ = std::move(sorted_hashes);
}

void ProteinProcessor::write_binary(
	const std::string path,
	ProteinSketchData& protein_sketch_data
	//std::vector<std::vector<uint64_t>>& hashes_
	) const {
	auto& hashes_ = protein_sketch_data.hashes;

	std::ofstream ofs(path, std::ios::binary);
	if (!ofs) {
		std::cerr << "Cannot write to: " << path << '\n';
		return;
	}

	std::shared_lock lock(result_mtx_);
	protein_sketch_data.saveConfig(ofs);
	for (int i = 0; i < config_.m; ++i) {
		const auto& row = hashes_[i];
		ofs.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(uint64_t));
		std::cout << "Writing " << row.size() * sizeof(uint64_t) << " bytes for row" << std::endl;
	}
	std::cout << "Sketch written to: " << path << '\n';
}

int ProteinProcessor::load_sequences(const std::string input_fa, int min_len, ProteinData& proteindata) {
	auto& seq_map_ = proteindata.sequence_map;
	auto& names_ = proteindata.names;

	gzFile fp = gzopen(input_fa.c_str(), "r");
	if (!fp) return 0;
	kseq_t* ks = kseq_init(fp);

	seq_map_.clear();
	names_.clear();
	int cnt = 0;

	while (true) {
		int len = kseq_read(ks);
		if (len < 0) break;
		if (len < min_len) continue;

		seq_map_[cnt] = ks->seq.s;
		names_.emplace_back(ks->name.s);
		++cnt;
	}

	kseq_destroy(ks);
	gzclose(fp);
	return cnt;
}

//std::vector<std::string> ProteinProcessor::get_names() const {
//	std::shared_lock lock(result_mtx_);
//	return names_;
//}
//
//std::vector<std::vector<uint64_t>> ProteinProcessor::get_hashes() const {
//	std::shared_lock lock(result_mtx_);
//	return hashes_;
//}
//
//const std::unordered_map<uint64_t, std::string>& ProteinProcessor::get_sequence_map() const {
//	return seq_map_;
//}
