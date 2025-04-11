#include "cluster.h"

Options options;
//SequenceDB seq_db;

void cluster::cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent)
{

//	Options options;
	SequenceDB seq_db;
//    options = Options();
    seq_db = SequenceDB();

    string db_out;

    if (seq.size() == 0) {
        std::cout<<"Empty input sequence"<<std::endl;
        return;
    }

    options.SetOptionsInputVector();
    options.Validate();

    db_out = options.output;
    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    seq_db.Readvector(seq, options);

    cout << "total seq: " << seq_db.sequences.size() << endl;
    seq_db.SortDivide(options);

    seq_db.DoClustering(options);

	seq_db.updateParent(parent);
}

int main()
{
	vector<Sequence_new> sequences;
	string filename = "test-update.fa";
	ifstream file(filename);
	std::string line;
    int currentId = -1; // 当前的 ID，初始值为无效值
    std::string currentData; // 用于临时存储序列内容

	while(getline(file, line)) {
        stringstream ss(line);
        int idCandidate;

		if(line.empty()) continue;

        if (line[0] == '>') {  // 如果以 '>' 开头
            if (currentId != -1 && !currentData.empty()) {
				// 保存之前的记录
				Sequence_new seq;
				seq.seq_id = currentId;
				seq.data = new char[currentData.size() + 1];
				strcpy(seq.data, currentData.c_str());
				sequences.push_back(seq);
            }
				// 更新当前 ID，并清空当前序列内容
			currentId++;
			currentData.clear();
        } else {
            currentData += line;  // 拼接序列行
        }

//        int idCandidate;
//		if (ss >> idCandidate && ss.eof()) {
//				// 如果读取成功并且整行是一个整数
//				if (currentId != -1 && !currentData.empty()) {
//						// 保存之前的记录
//						Sequence_new seq;
//						seq.seq_id = currentId;
//						seq.data = new char[currentData.size() + 1];
//						strcpy(seq.data, currentData.c_str());
//						sequences.push_back(seq);
//				}
//				// 更新当前 ID，并清空当前序列内容
//				currentId = idCandidate;
//				currentData.clear();
//		} else {
//				// 当前行是序列的一部分，追加到 currentData 中
//				currentData += line;
//		}
	}

    if (currentId != -1 && !currentData.empty()) {
        Sequence_new seq;
        seq.seq_id = currentId;
        seq.data = new char[currentData.size() + 1];
        std::strcpy(seq.data, currentData.c_str());
        sequences.push_back(seq);
    }

    file.close();
	cout << sequences.size() << endl;
	vector<int> parent(currentId + 1, 0);
	cluster cls;
	cls.cdhit_cluster(sequences, parent);
//	for(int i = 0 ; i <= currentId; i++){
//		cout << i << " " << parent[i] << endl;
//	}
	vector<int> t_parent1(currentId + 1, 1);
	vector<int> t_parent2(currentId + 1, 0);

	cls.cdhit_cluster(sequences, t_parent1);
	cls.cdhit_cluster(sequences, t_parent2);
	for(int i = 0; i <= currentId; i++){
		cout << parent[i] << " " << t_parent1[i] << t_parent2[i] << endl;
	}
}
