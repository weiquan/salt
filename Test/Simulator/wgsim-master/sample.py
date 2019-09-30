import sys
import random
import time
N = 100000 #2M
usage = 'sample.py [Read1] [Read2] '
def readfq(fp):
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



if __name__ == '__main__':
    if len(sys.argv) < 2:
        print usage
        exit(1)
    #count reads number in in1.fastq
    t = time.time()
    fp = open(sys.argv[1], 'r')
    n = 0
    for name, seq, qual in fp:
        n += 1
    fp.close()
    if n < N:
        print >>sys.stderr, 'seq number smaller than sample number!'
        exit(1)
    print "count reads number use time: %lu sec"%(time.time()-t)
    #sample
    t = time.time()
    population = [i for i in range(n)]
    sample = random.sample(population, N)
    sample.sort()
    print "sample use time: %lu sec"%(time.time()-t)
    #choice reads from in1.fastq
    for k in xrange(1,3) 
        t = time.time()
        fp = open(sys.argv[k], 'r')
        fp_o = open(sys.argv[k]+'.sample', 'w')
        
        i = 0
        j = 0
        for name, seq, qual in fp:
             
            
            if j == len(sample):
                break
            if i == sample[j]:
                print '@'+name
                print seq 
                print '+'
                print qual
                j += 1
            i += 1
        if j < len(sample) and i >sample[j]:
            print >>sys.stderr, 'error'
        if i %1000 == 0:
            print >>sys.stderr, '%u reads passer!'%(i)
        fp.close()
        fp_o.close()

        print "sample %lu lines from %lu lines use time: %lu sec"%(N, n, time.time()-t)
       
