import optparse
import itertools
import sys
if __name__ == '__main__':
    #parse arguments
    usage = 'usage: %prog [Options] <dbsnp.txt>'
    parser = optparse.OptionParser(usage)
    parser.add_option('-n', '--ncbi', action = 'store_true',  dest='refNcbi', default ='False', help='reference genomic seq from NCBI [UCSC]')
    parser.add_option('-s', '--rs', action = 'store_true',  dest='allelesFromRs', default ='False', help='alleles from rs-fasta files [False]')
    parser.add_option('-f', '--frequency', action = 'store', type = 'float', dest='allelesFromFreq', default = 0.1, help='alleles from frequency data and set min alleles freq [0.1]')
    parser.add_option('-d', '--debug', action = 'store_true', dest='debug', default = False, help = 'print log info to stderr [False]')
    opts, args = parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit(1) 
    #core loop
    fp = open(args[0], 'r')
    for line in fp:
        line = line.strip()
        words = line.split('\t')
        chrom = words[1]
        chromStart = int(words[2])
        chromEnd = int(words[3])
        variant_type = words[11]
        if variant_type != 'single' or chromEnd - chromStart !=1:
            if opts.debug:
                print >>sys.stderr, '[Warning] variant_type is not single' 
                print >>sys.stderr, line  
            continue
        strand = words[6]
        ref = None
        if opts.refNcbi:
            ref = words[7]
        else:#ref from ucsc
            ref = words[8]
        min_freq = opts.allelesFromFreq
        neclotide_types = {'A':False, 'C':False, 'G':False, 'T':False}
        neclotide_complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        #
        neclotide_types[ref] = True
        if opts.allelesFromRs:#alleles from rs fasta
            alleles = words[9]
            for nt in alleles.split('/'):
                if len(nt) != 1:
                    if opts.debug:
                        print >>sys.stderr, '[Warning] allele len >1, skip this line'
                        print >>sys.stderr, line  
                    continue
                if strand == '+':
                    neclotide_types[nt.upper()] = True
                else:
                    comp = neclotide_complement[nt.upper()]
                    neclotide_types[comp] = True
        else:#alleles from freq dat
            alleles = words[23]
            for nt, f in itertools.izip( alleles.split(','), words[25].split(',') ):
                if len(nt) != 1 or f < min_freq:
                    if opts.debug:
                        print >>sys.stderr, '[Warning] allele len >1, skip this line' 
                        print >>sys.stderr, line  
                    continue
                if strand  == '+':
                    neclotide_types[nt.upper()] = True
                else:
                    comp = neclotide_complement[nt.upper()]
                    neclotide_types[comp] = True
        neclotide_string = ''
        for nt in ['A', 'C', 'G', 'T']:
            if neclotide_types[nt]:
                neclotide_string += nt+'/'
        neclotide_string = neclotide_string.rstrip('/')
        print '%s\t%u\t%s\t%s'%(chrom, chromEnd, neclotide_string, ref)
    fp.close()