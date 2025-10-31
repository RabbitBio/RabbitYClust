# `RabbitYClust`
## Installtion

### Install from source code

#### Dependancy
cmake version 3.0
c++17

#### Libs Dependancy
git clone git@github.com:RabbitBio/Rabbitlibcdhit.git
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/libcdhit
make install

#### Compile and install
git clone git@github.com:RabbitBio/RabbitYClust.git
cp ~/libcdhit/lib . -r
cp ~/libcdhit/include . -r
mkdir build
cd build
cmake ..
make

## Usage
```bash
# build MinHash sketch for a protein FastA file
Usage: ./yclust sketch <-i> input-FA-filename <-o> output-sketch-result-filename [Options]
Required:
	-i,--input			  set the input file name, fasta or gziped fasta formats
	-o,--skechfile-output set the output file name of the MinHash Sketches
Options:
	-h,--help	look for more informations
	-x,--xxhash-flag 	set the flag on to use xxhash to generate MinHash sketch, default is false
	-t,--threads INT	set the thread number, default 1 thread
	--min-len    INT	set the filter minimum length , protein length less than minLen will be ignore, default 10;
	-k,--kmer-size INT  set the kmer size, default 8
	-m,--sketch-size INT	set the number of hash functions will be used, default 15
```

```bash
# clustering
Usage: ./yclust cluster <-i> input-FA-filename <-o> output-sketch-result-filename [Options]
Required:
	-i,--input			  set the input file name, fasta or gziped fasta formats
	-o,--skechfile-output set the output file name of the MinHash Sketches
	-S,--skech-file-name  set the file name of MinHash Sketch
Options:
	-h,--help
	-t,--threads INT	set the thread number, default 1 thread
	-m,--sketch-size INT	set the number of hash functions will be used, default 15
	-r,--r-size INT set the number of block
	-c,--cluster-on 		enable clustering after grouping
	-f,--final-cluster-off  turn off the final clustering, only work with -c is on;
```
## Example

## Output


