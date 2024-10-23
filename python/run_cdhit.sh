#!/bin/bash

# 定义输入和输出目录
if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <input_directory> <output_file> [cd-hit]"
  exit 1
fi
INPUT_DIR="$1"     # 输入包含 .fa 文件的目录
OUTPUT_DIR="$2" # 输出目录，存储 CD-HIT 聚类结果
CDHIT="$3" # cd-hit程序路径
THREAD_NUM="$4" # 线程数
IDENTITY_THRESHOLD="0.9"      # CD-HIT 的聚类相似性阈值

# 创建输出目录（如果不存在）
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# 遍历输入目录中的所有 .fa 文件
#for file in "$INPUT_DIR"/*.fa; do
#find "$INPUT_DIR" -name "*.fa" | xargs -P "$THREAD_NUM" -I {} bash -c '
#  file="{}"
#find "$INPUT_DIR" -name "*.fa" | parallel -j 4 --halt soon,fail=1 '
#  file="{}"

counter=0
echo "Time needed for cd-hit clustering with $4 threads..."
time{
for file in "$INPUT_DIR"/*.fa; do
  if [ -f "$file" ]; then
    # 提取文件名（不包括路径和扩展名）
    filename=$(basename -- "$file")
    base="${filename%.*}"

    # 输出文件路径
    output_file="$OUTPUT_DIR/${base}_clustered.fa"

    # 运行 CD-HIT，进行聚类
    "$CDHIT" -i "$file" -o "$output_file" -c "$IDENTITY_THRESHOLD" &
	((counter+=1))
	if [ counter == $THREAD_NUM ]; then 
		wait
		counter=0
	fi

    echo "Finished clustering $file. Output saved to $output_file"
  fi
done
}

echo "Time needed for aggregating multiple clustering results..."
time{
if [ -z "$5" ]; then # 最终输出结果
    FINAL_OUTPUT_FILE="result.clstr"
else
    FINAL_OUTPUT_FILE="$5"
fi
}
for file in "$OUTPUT_DIR"/*.clstr; do
  if [ -f"$file" ]; then
	cat "$file" >> "$FINAL_OUTPUT_FILE"
	echo "Meged: $file"
  fi
done

