#echo 'Download dbsnp135(common snp) from NCBI...'
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/common_all_20150603.vcf.gz
#echo 'Reference from NCBI...'
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

###Human 
REF_RELEASE=hg19
UCSC_HG19=ftp://hgdownload.cse.ucsc.edu/goldenPath/${REF_RELEASE}/chromosomes/chr*
wget UCSC_HG19


DBSNP_RELEASE=142
SNP_FILE=snp${DBSNP_RELEASE}Common.txt
UCSC_COMMON_SNP=http://hgdownload.cse.ucsc.edu/goldenPath/${REF_RELEASE}/database/${SNP_FILE}
wget UCSC_COMMON_SNP


###Mouse
REF_RELEASE=mm10
UCSC_HG19=ftp://hgdownload.cse.ucsc.edu/goldenPath/${REF_RELEASE}/chromosomes/chr*
wget UCSC_HG19

DBSNP_RELEASE=138
SNP_FILE=snp${DBSNP_RELEASE}Common.txt
UCSC_COMMON_SNP=http://hgdownload.cse.ucsc.edu/goldenPath/${REF_RELEASE}/database/${SNP_FILE}
wget UCSC_COMMON_SNP


