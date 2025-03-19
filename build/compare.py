def compare_files(file1_path, file2_path):
    # 读取文件1中的所有行
    with open(file1_path, 'r') as file1:
        file1_lines = set(file1.readlines())
    
    # 读取文件2中的所有行
    with open(file2_path, 'r') as file2:
        file2_lines = set(file2.readlines())
    
    # 打印出文件1中独有的行
    print("Lines present in file1 but not in file2:")
    for line in file1_lines - file2_lines:
        print(line.strip())
    
    # 打印出文件2中独有的行
    print("\nLines present in file2 but not in file1:")
    for line in file2_lines - file1_lines:
        print(line.strip())

# 调用函数，比较两个文件
file1_path = 'log'  # 替换成你的文件路径
file2_path = 'log3'  # 替换成你的文件路径
compare_files(file1_path, file2_path)
