
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

SEGMENT_LEN = 250
def reverse_complement(seq):
    complement_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C', '-':'-', 'N':'N'}
    rev_comp = ''.join([complement_map[c] for c in seq.upper() if c in complement_map])[::-1]
    if len(rev_comp) > 0: return rev_comp
    return None
def extract_indel(fn_snp):
    
    variants_map = {}
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
        
        observed = observed.upper()
        if strand == '-': 
            observed = [reverse_complement(x) for x in observed.split('/') if reverse_complement(x) is not None]
        else:
            observed = observed.split('/')
        variants_map.setdefault(chrom, []).append((chromStart, chromEnd, variantType, refUCSC, observed))
    fp.close()
    return variants_map
def main(opts, args):
   
    fp_out_genome = open(args[2]+'.fa', 'w')
    fp_out_indel = open(args[2]+'.indel.fa', 'w')

    print >>sys.stderr, 'Generating snp hash map...'
    snp = extract_indel(args[1])
    fp_genome = open(args[0], 'r') 
    
    for track_name, seq, _ in readfq(fp_genome):
        l_seq = len(seq)
        print >>sys.stderr, 'Generate snp in {}'.format(track_name)
        rname = track_name.split()[0]
        if rname not in snp:
            print >>sys.stderr, 'chrom {} no variants in file {}'.format(rname, args[1])
            continue
        print >>fp_out_genome, '>'+rname
        for i in xrange(0, len(seq), 60):
            if i + 60 < len(seq):
                print >>fp_out_genome, seq[i:i+60]
            else:
                print >>fp_out_genome, seq[i:]
        for v in snp[rname]:
            print >>sys.stderr, rname, seq[v[0]-1], seq[v[0]], seq[v[0]+1], v
            if seq[v[0]-1].upper() not in v[4]: continue 
            start, end, variantType, ref, observed = v
            alt =[i for i in observed if i != ref]
            if variantType == 'deletion':
                print >>fp_out_indel, '>{}_{}_{}_{}'.format(rname, start, end, variantType)
                print >>fp_out_indel, seq[max(0, start-1-SEGMENT_LEN): start-1]+seq[end:min(end+SEGMENT_LEN, l_seq)] 
            if variantType == 'insertion':
                for ins in alt:
                    print >>fp_out_indel, '>{}_{}_{}_{}'.format(rname, start, end, variantType)
                    print >>fp_out_indel, seq[max(0, start-1-SEGMENT_LEN): start-1]+seq[end:min(end+SEGMENT_LEN, l_seq)] 



    fp_genome.close()
    fp_out_indel.close()
    fp_out_genome.close()


if __name__ == '__main__':
    usage = 'usage: %prog [Options] <reference.fa> <reference.snp> <Prefix>'
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.print_help()
        sys.exit(1)
 
    
    main(opts, args)
