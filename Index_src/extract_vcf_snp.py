
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

neclotide_complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def extract_variants(fn_snp):
	min_var_distance = 20    
    var = {}
    fp = open(fn_snp, 'r')
    for line in fp:
        if line[0] == '#': continue
        line = line.strip()
        #print line
        try: 
            (chrom, pos, rs_id, ref, alt, qual, filter, info ) = line.split('\t')[:8]
        except ValueError:
            continue
        var = []
        for v in alt.split(','):
        	if(len(ref) == len(v) == 1):#SNP
        		var.appped(v.upper())
        	elif(len(ref) == 2 and len(v) == 1):#DEL
        		var.appped(v.upper())
        	elif(len(ref) == 1 and len(v) == 2):#INS
        		var.appped(v.upper())
    
       	alt = '|'.join(var)
        #print '%s\t%u\t%s\t%s'%(chrom, pos, ref, alt)
        snp.setdefault(chrom, []).append((pos, ref, alt))
    fp.close()
    return snp
def gen_local_variants_string(ref, snp_map):
    max_snp_distance = 20
    snp_chain = []
    for snp in snp_map:
        if len(snp_chain) == 0:
            snp_chain.append(snp)
        elif snp[0] - snp_chain[-1][0] <= max_snp_distance:
            snp_chain.append(snp)
        else:
            #[to do] :generate strings from snp_chain
            print snp_chain
            snp_chain = []

def main(opts, args):
    kmer = 25
    fp_out_genome = open(args[2]+'.fa', 'w')
    fp_out_snp = open(args[2]+'.snp', 'w')

    print >>sys.stderr, 'Generating snp hash map...'
    snp_map = extract_variants(args[1])
    fp_genome = open(args[0], 'r') 
''' 
    for track_name, seq, _ in readfq(fp_genome):
        print >>sys.stderr, 'Generate snp in {}'.format(track_name)
        rname = track_name.split()[0]
        if rname not in snp:
            print >>sys.stderr, 'chrom {} no variants in file {}'.format(rname, args[1])
            continue
        print >>fp_out_genome, '>'+rname
        g = []
        last_pos = 0
        for v in snp[rname]:
            pos, ref, alt = v
            if(len(g) == 0):
                g.append( seq[max(pos-25, 0):pos]+"["+ref+'|'+alt+']')
            elif (pos-last_pos > kmer):
                print >> fp_genome, g
                g = []
            elif( pos - last_pos <= kmer):
                g.append( seq[last_pos:pos]+"["+ref+'|'+alt+']')
            last_pos = pos
'''
    for track_name, seq, _ in readfq(fp_genome):
        if(track_name not in snp_map): continue:
        gen_local_variants_string(seq, snp_map[track_name]);
            


    fp_genome.close()

    fp_out_snp.close()
    fp_out_genome.close()


if __name__ == '__main__':
    usage = 'usage: %prog [Options] <reference.fa> <reference.snp> <Prefix>'
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.print_help()
        sys.exit(1)
 
    
    main(opts, args)
