"""Extract unmapped reads from samfile."""
import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1], 'r')
for read in samfile:
    if read.is_unmapped:
        print '>' + read.query_name
        print read.query_sequence
        # print '+'
        # print read.query_qualities
samfile.close()
