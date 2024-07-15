# TODO List
The basic idea of RabbitYClust is inspired by the KHF Minhash algorithm. In MinHash algorithm we can easily find the candidate matches by finding the hitting hashes. Thus, we can simply grouping the input sequences by finding the candidate (to be likely belong to a cluster) sequences with bounded errors.

The errors can be calculated by $(1-s^r)^b$, where $s$ is the targeting similarity, $b$ is the number of bands and $r$ is the hashes in each band.

- [ ] Python testing scripts of the grouping method.
- [ ] Optimization of the scripts using 'Union-Find' data structure and struct sort.
- [ ] Generate hashes based on our OrderMinHash implementation.
- [ ] Testing the accuracy of a subset nr dataset or some simulated data.
- [ ] Read the paper of Linclust and make the method of linclust clear