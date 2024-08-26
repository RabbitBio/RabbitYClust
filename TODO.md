# TODO List
The basic idea of RabbitYClust is inspired by the KHF Minhash algorithm. In MinHash algorithm we can easily find the candidate matches by finding the hitting hashes. Thus, we can simply grouping the input sequences by finding the candidate (to be likely belong to a cluster) sequences with bounded errors.

The errors can be calculated by $(1-s^r)^b$, where $s$ is the targeting similarity, $b$ is the number of bands and $r$ is the hashes in each band.

- [x] Python testing scripts of the grouping method.
- [ ] Optimization of the scripts using 'Union-Find' data structure and struct sort.
- [x] Generate hashes based on our OrderMinHash implementation.
- [ ] Testing the accuracy of a subset nr dataset or some simulated data.
- [x] Read the paper of Linclust and make the method of linclust clear.
- [ ] Add support aaHash (kind of rolling hash).
- [ ] FIXIT: The single-threaded hash generation is still slow (82MB fasta file consumes about 35 secs).
- [ ] FIXIT: Return the reversible encoded int value instead of hash value
- [ ] FIXIT: Add option to ignore small peptide sequences

