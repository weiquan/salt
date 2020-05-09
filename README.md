## Getting Started
```sh
git clone https://github.com/weiquan/salt.git
cd salt && make
# run test
cd Test/Run_test
./run_test.sh 
# download GRCh38 and common SNPs
./Script/download_GRCh38_snp144.sh
python2 Script/extract_snp.py GRCh38.fa snp144Common.txt GRCh38
# build index for GRCh38 and SNP144
./Bin/salt-idx -k 21 GRCh38.fa GRCh38.snp GRCh38.idx
# align single-end illumina reads
./Bin/salt GRCh38.idx Reads.fq >aln.sam
# align paired-end illumina reads
./Bin/salt -p GRCh38.idx Reads1.fq Reads2.fq >aln.sam
```