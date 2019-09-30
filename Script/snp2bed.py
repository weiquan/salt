
import optparse
import sys
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
def extract_variants(fn_snp):
    
    snp = {}
    fp = open(fn_snp, 'r')
    for line in fp:
        if line[0] == '#': continue
        line = line.strip()
        #print line
        try: 
            (_, chrom, chromStart, chromEnd, name, score, strand, 
            refNCBI, refUCSC, observed, molType, variantType) = line.split('\t')[:12]
        except ValueError:
            continue
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
    
        if variantType != 'single' or chromEnd - chromStart !=1:
            continue

        neclotide_types = {'A':False, 'C':False, 'G':False, 'T':False}
        neclotide_complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        #
        alleles = observed.upper()
        for nt in alleles.split('/'):
            if len(nt) != 1:
                if opts.debug:
                    print >>sys.stderr, '[Warning] allele len >1, skip this line'
                    print >>sys.stderr, line  
                continue
            if nt not in ['A', 'C', 'G', 'T']: continue
            if strand == '+':
                neclotide_types[nt] = True
            else:
                comp = neclotide_complement[nt]
                neclotide_types[comp] = True
     
        neclotide_string = ''
        for nt in ['A', 'C', 'G', 'T']:
            if neclotide_types[nt]:
                neclotide_string += nt+'/'
        neclotide_string = neclotide_string.rstrip('/')
        #print '%s\t%u\t%s\t%s'%(chrom, chromEnd, neclotide_string)
        snp.setdefault(chrom, []).append((chromEnd, neclotide_string))
    fp.close()
    return snp
def main(opts, args):
   
    #fp_out_genome = open(args[2]+'.fa', 'w')
    fp_out_snp = open(args[2], 'w')

    print >>sys.stderr, 'Generating snp hash map...'
    snp = extract_variants(args[1])
    fp_genome = open(args[0], 'r') 
    
    for track_name, seq, _ in readfq(fp_genome):
        print >>sys.stderr, 'Generate snp in {}'.format(track_name)
        rname = track_name.split()[0]
        if rname not in snp:
            print >>sys.stderr, 'chrom {} no variants in file {}'.format(rname, args[1])
            continue
        for v in snp[rname]:
            #print >>sys.stderr, rname, seq[v[0]-1], seq[v[0]], seq[v[0]+1], v
            if seq[v[0]-1].upper() not in v[1]: continue 
            print >>fp_out_snp, '{}\t{}\t{}'.format(rname, v[0]-1, v[0])
    fp_genome.close()
    fp_out_snp.close()



if __name__ == '__main__':
    usage = 'usage: %prog [Options] <reference.fa> <reference.snp> <reference.bed>'
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.print_help()
        sys.exit(1)
 
    
    main(opts, args)
