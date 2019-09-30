
import optparse
import sys, re
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
        snp[(chrom, chromEnd)] = name 
          

        #

       #print '%s\t%u\t%s\t%s'%(chrom, chromEnd, neclotide_string)

    fp.close()
    return snp
'''
    5M2I2M2D1M

'''
def t_shift(cigar, q_shift):
    pattern = re.compile(r'(\d+)([SMIDX=])')
    cigars = pattern.findall(cigar)
    t = 0
    q = 0
    for x in cigars:
        n = int(x[0])
        c = x[1]
        if c == 'S': continue
        if  q > q_shift: break
        if c == 'I':
            q += n
        elif c == 'D':
            t += n
        elif c in ['M', 'X', '=']:
            d = min(n, q_shift-q+1)
            #print >>sys.stderr, '123 {} {}'.format(q_shift, q)
            t += d
            q += d
    return t-1
def main(opts, args):
    import time
    print >>sys.stderr, 'build snp db'
    t = time.clock()
    snp = extract_variants(args[1])
    print >>sys.stderr, 'find db'
    print >>sys.stderr, 'time escaped {}'.format(time.clock()-t)
    
    t = time.clock()
    sam = open(args[0])
    pattern = re.compile(r'(?<=XV:i:)\S+')
    for track in sam:
        track = track.strip()
        if track[0] == '@': 
            print track
            continue

        track_fields = track.split()
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = track_fields[:11]
        if len(track_fields) > 11: 
            rs = ''
            snp_shift_list = pattern.search(track)
            if snp_shift_list is not None:
                snp_shift_list = snp_shift_list.group(0)
                #print >>sys.stderr, snp_shift_list
                for snp_shift in snp_shift_list.split(','):
                    #snp_shift =int(snp_shift)+ seq_start(cigar)
                    rs_pos = t_shift(cigar, int(snp_shift))+int(pos)
                    if (rname, rs_pos) not in snp:
                        #print >>sys.stderr, rname, rs_pos,snp_shift 
                        #print >>sys.stderr, t_shift(cigar, int(snp_shift))
                        print >>sys.stderr, track
                        exit(1)
                    rs_id = snp[(rname, rs_pos)]
                    rs += rs_id+','
                if len(rs) > 0:
                    track += '\tRS:Z:'+rs[:-1]
        print track 
    sam.close()

    print time.clock()-t
    



if __name__ == '__main__':
    usage = 'usage: %prog [Options] <SAM> <reference.snp>'
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit(1)
 
    
    main(opts, args)
