## 在cdhit修改的代码：

cdhit-common.h / cdhit-common.c++:

1. 在 `Options` 中 新增 `bool input_vector` 这一选项用来实现非文件读入时的 `Validate` 时的不一致性。

2. 在 `Options` 中增加函数 `bool Options::SetOptionsInputVector()` 来修改 `input_vector` 的值 `false(default)->true`

3. 在 `Options` 中修改 `Validate` 函数，在`input_vector == true` 时屏蔽对 `input` 的检查

4. 新增 `Readvector` 函数 ，实现传入`vector<Input_Sequence*>` 时与 `SequenceDB` 中的 `vector<Sequence*>` 的转换，并处理相关内容使得其可以完成后续的cdhit代码

5. 新增 `WriteClustervector` 函数，实现当传入`vector<Input_Sequence*>`时，可以与文件输入时保持一致的代表序列输出。
    （后续如果需要修改代表序列的输出，比如质量分数等信息可以直接从此处进行修改，在输出的部分继续输出 `Input_Sequence` 的其他信息）

Makefile:

1. 新增对 `cdhit-lib.o` 的编译

2. 新增对 `libcdhitlib.a` 的编译

## 在cdhit中新增的代码：

1. 新增 `input_sequnece.h`，来实现对 Vector 数据的接收，现在其中只有 `char* identifier` 和 `char* data` 两个属性。如果后期有质量分数相关属性可以直接在此处新增一些东西，输出直接修改 `WriteClustervector` 即可。（如果要实现和原函数一致的，保持源文件一样的输入格式，可以采用存储一个 `char** input_data` 的形式，保持多行的格式，同时修改 `WriteClustervector` ，在输出时输出 `char** input_data`，即可。

2. 新增 `cdhit-lib.h` 和 `cdhit-lib.c++` 文件用于包装静态库，其中实现 `cdhit_vector` 和 `cdhit-file` 两个函数实现对外的接口。

## 静态库需要包含的文件：

cdhit-common.h

cdhit-lib.h

input_sequence.h

libcdhitlib.a

ps. 如果单独编译，记得 加上 `-lz` 和 `-fopenmp`

## 函数用法：

1. `void cdhit_vector(std::vector<Input_Sequence*>& seq, std::string output, std::string command)`, 将所需的数据利用 `vector<Input_Sequence*>` 进行传入，以字符串的格式填入输出文件，剩余的command以字符串的形式直接填入对应的cdhit的相关参数，例如 `std::string command = "-t 0.9 -M 0 -T 4"` 这样来表示 threshold, Memory, thread_num 等值的信息，即除去 `-i` 和 `-o` 以外任何格式。

2. `void cdhit_file(std::string input, std::string output, std::string command)`，以字符串的格式填入输入和输出文件，剩余的command以字符串的形式直接填入对应的cdhit的相关参数，例如 `std::string command = "-t 0.9 -M 0 -T 4"` 这样来表示 threshold, Memory, thread_num 等值的信息，即除去 `-i` 和 `-o` 以外任何格式。

ps. 详细可以去看main.cpp 中的调用，并且给出了在如下格式中可以成功编译的 Makefile 和 编译指令

```
/project
	main.cpp
    /include
    	cdhit-common.h
    	cdhit-lib.h
    	cdhit-utility.h (可有可无)
    	input_sequence.h
    /libs
    	libcdhitlib.a
    Makefile

```

`g++ main.cpp -I./include -L./libs /home/user_home/maguiliang/codes/c_cpp/cdhit_test/libs -fopenmp -lz -o main`