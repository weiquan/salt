usage ='VCF2snpaln.py <file.vcf>'
import sys
nt2number  = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
if __name__ == '__main__':
    if len(sys.argv) < 2:
        #print >>sys.stderr, 'CMD ERROR!'
        print >>sys.stderr, usage
        sys.exit(1)
    fp = open(sys.argv[1])
    for line in fp:
        line = line.strip()
        #skip meta-info and header
        if line[0] =='#':
            continue
        words = line.split('\t')
        chrom = words[0] 
        pos = int(words[1])
        ref = words[3]
        alt = ''
        if len(ref) > 1: continue
        alt_map = [0 for i in xrange(4)]
        alt_map[nt2number[ref]] =1
        for c in words[4].split(','):
            if len(c) >1:
                continue
            alt_map[nt2number[c]] = 1
        for i in xrange(4):
            if alt_map[i] != 0:
                alt += 'ACGT'[i]+'/'
        alt = alt[:-1]#trip last '/'
        if len(alt) >1:
            print '{0}\t{1}\t{2}\t{3}'.format(chrom, pos, alt, ref)
    
    fp.close()
