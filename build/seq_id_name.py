import re

def load_sequence_ids(sequence_id_file):
    seq_dict = {}
    with open(sequence_id_file, 'r') as f:
        for line in f:
            parts = line.rsplit(" ", 1)  # 拆分序列和ID
            if len(parts) == 2:
                seq, seq_id = parts[0], parts[1].strip()
                seq_dict[seq] = seq_id
    return seq_dict

def match_sequences(huge_fa_file, seq_dict):
    matches = []
    with open(huge_fa_file, 'r') as f:
        header, sequence = None, ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header and sequence in seq_dict:
                    matches.append((header, seq_dict[sequence]))
                header = re.split(r'\s|\x01', line)[0][1:]  # 提取名称部分
                sequence = ""
            else:
                sequence += line
        if header and sequence in seq_dict:
            matches.append((header, seq_dict[sequence]))
    return matches

def save_results(matches, output_file):
    with open(output_file, 'w') as f:
        for header, seq_id in matches:
            f.write(f"{header}\t{seq_id}\n")

def main():
    sequence_id_file = "sequence-id.txt"
    huge_fa_file = "../../pure-aahash/huge.fa"
    output_file = "huge_results_1.txt"
    
    seq_dict = load_sequence_ids(sequence_id_file)
    matches = match_sequences(huge_fa_file, seq_dict)
    save_results(matches, output_file)
    
    print(f"匹配完成，结果已保存到 {output_file}")

if __name__ == "__main__":
    main()