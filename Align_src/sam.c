#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "bntseq.h"
#include "sam.h"
#include "kstring.h"
#include "editdistance.h"
#define SAM_FLAG_PAIRED  0x0001
#define SAM_FLAG_PROPER_PAIRED  0x0002
#define SAM_FLAG_UNMAPED  0x0004
#define SAM_FLAG_MATE_UNMAPPED 0x0008
#define SAM_FLAG_REVERSE 0x0010
#define SAM_FLAG_MATE_REVERSE 0x0020
#define SAM_FLAG_MATE_READ1 0x0040
#define SAM_FLAG_MATE_READ2 0x0080
#define SAM_FLAG_SECANDARY_ALN 0x0100
#define SAM_FLAG_QC 0x0200
#define SAM_FLAG_PCR_DUP 0x0400


#define CIGAR_ERRO 1
#define CIGAR_CORRECT 0

void sam_add_xa(kstring_t *s, index_t *index, query_t *query, int is_cigar);
void sam_add_md_nm(kstring_t *s, index_t *index, query_t *q);
int check_cigar(char *cigar, uint32_t l_cigar, uint32_t l_aln)
{
    //assert(cigar);

    uint32_t tot_l= 0;
    
    long int n; char *pc;
    pc = cigar;
    while(*pc){
        if(isdigit(*pc)){
            n = strtol(pc, &pc, 10);
            if(*pc == 'M' || *pc == 'I') {
                tot_l += n; 
            }
            ++pc; 
        
        }
    }


    if(tot_l != l_aln) return CIGAR_ERRO;
    
    return CIGAR_CORRECT; 


}

void aln_samhead(const opt_t *opt, bntseq_t *bntseq)
{
    int i;
    kstring_t *s;
    s = calloc(1, sizeof(kstring_t));
    //@HD   VN:format version  
    ksprintf(s,"@HD\tVN:ec1fec2\tSO:unsorted\n");
    //@SQ   SN:Ref seq name LN:Ref seq len
    for(i = 0; i < bntseq->n_seqs; ++i){
        ksprintf(s,"@SQ\tSN:%s\tLN:%d\n", bntseq->anns[i].name, bntseq->anns[i].len);
    }
    //@RG   ID:Read group id
    ksprintf(s, "@RG\tID:%s\n", opt->rg_id);
 
    

    


    //@PG   ID:program id   PN:program name CL:cmd line DS:description VN:program version
    time_t t0 = time(NULL);
    struct tm *t = localtime(&t0);
    ksprintf(s, "@PG\tID:snpaln\tPN:snpaln\tCL:\"%s\"\tDS:%d-%d-%d\tVN:0.1beta",opt->cmd->s, t->tm_year+1900, t->tm_mon+1, t->tm_mday);

    printf("%s\n", s->s);

    free(s->s);
    free(s);
}

#define INIT_STR_LEN 256
void aln_samse(index_t *index, query_t *query, const aln_opt_t *opt)
{
    
    bntseq_t *bntseq = index->bntseq; 
    int i;
    const uint8_t *seq, *rseq;
    const char *name;
    
    //const hit_t *hit;
    uint16_t flag = 0;
    kstring_t *s = query->sam; 

    seq = query->seq;
    rseq = query->rseq;
    name = query->name;
    
    //SAM 
    //Qname FLAG Rname Pos Mapq Cigar Rnext Pnext Tlen Seq Qual
    
    //unalned
    if(query->pos == 0xFFFFFFFF){
        ksprintf(s, "%s\t", name);//name
        flag |= SAM_FLAG_UNMAPED;
        ksprintf(s, "%u\t", flag);//flag
        //Rname=*,  Pos=0,  Mapq=0, Cigar=*, Rnext=*, Pnext=0, Tlen=0
        ksprintf(s, "*\t0\t0\t*\t*\t0\t0\t");
        //Seq 
        for(i = 0; i < query->l_seq; ++i){
            //kputc("ACGT"[seq[i]],s);
            ksprintf(s,"%c", "ACGTN"[seq[i]]);

        }
        //Qual 
        if(query->qual) {
            ksprintf(s, "\t%s", query->qual); 
        } else{
            ksprintf(s, "\t*"); 
        }
        return; 
    } 


    //alned 
    ksprintf(s, "%s\t", name);//name      
    if(query->strand) flag |= SAM_FLAG_REVERSE;
    ksprintf(s, "%u\t", flag);//flag
    int rid;
    bns_coor_pac2real(bntseq, query->pos, query->l_seq, &rid);   
    ksprintf(s, "%s\t", bntseq->anns[rid].name);//Rname
    ksprintf(s, "%lu\t", query->pos - bntseq->anns[rid].offset +1);//Pos
    ksprintf(s, "%u\t", query->mapq);//Mapq
    //for(i = 0; i < query->l_cigar; ++i) ksprintf(s, "%u%c",query->cigar[i]>>4,"MID"[query->cigar[i]&15]);//cigar
    //for(i = 0; i < query->l_cigar; ++i) ksprintf(s, "%u%c",__parse_cigar_num(query->cigar[i]),__parse_cigar_op(query->cigar[i]));//cigar
    ksprintf(s, "%s", query->cigar->s);
    ksprintf(s, "\t*\t0\t0\t");//Rnext=*, Pnext=0, Tlen=0

    //print seq
    if(query->strand) {
        for(i = 0; i < query->l_seq; ++i){
            //kputc("ACGT"[rseq[i]],s);
            ksprintf(s, "%c","ACGTN"[rseq[i]]);

        }
        kputc('\t', s);
        //if(strlen((const char*)query->qual)) {
        if(query->qual != NULL ) {
            for(i = query->l_seq -1; i >= 0; --i){
                //kputc(query->qual[i],s);
                ksprintf(s, "%c",query->qual[i]);

            }
        } else{
            //kputc('*',s);
            ksprintf(s, "*");

        }
    } else{
        for(i = 0; i < query->l_seq; ++i){
            //kputc("ACGT"[seq[i]],s);
            ksprintf(s, "%c","ACGTN"[seq[i]]);

        } 
        //kputc('\t',s);            
        ksprintf(s, "\t");

        if(strlen((const char*)query->qual)) {
            ksprintf(s, "%s", query->qual); 
        } else{
            ksprintf(s, "*");
        }
    } 
     
    sam_add_xa(s, index, query, opt->print_xa_cigar);
    if(opt->print_nm_md) sam_add_md_nm(s, index, query); 
    if(opt->rg_id != NULL) ksprintf(s, "\tRG:Z:%s", opt->rg_id);
}

#define STREAM_STD 1
#define STREAM_KSTRING 2
void sam_add_xa(kstring_t *s, index_t *index, query_t *query, int is_cigar)
{
    bntseq_t *bntseq = index->bntseq;
    int rid = 0;
    int i;

    uint32_t last_pos = query->pos, primary_pos = query->pos;
    int is_print = 1;
    int strand;
    for(strand=0; strand<2; ++strand){
        for(i = 0; i < query->hits[strand].n; ++i){
            hit_t *hit = query->hits[strand].a+i;
            uint32_t p, l, n;
            p = hit->pos;
            n = hit->n_diff;
            l = query->l_seq;
            if(p > index->bntseq->l_pac){
                fprintf(stderr, "[%s]: pos > reference len(%u > %u)", __func__, p, index->bntseq->l_pac); 
            } 
            //if(p == last_pos) continue;
            if(p == primary_pos) continue;
            if (is_print) { 
                ksprintf(s, "\tXA:Z:");
                is_print = 0;
            }
            bns_coor_pac2real(bntseq, p, l, &rid); 
            ksprintf(s, "%s,", bntseq->anns[rid].name);//Rname
            ksprintf(s, "%c%lu,", "+-"[strand], p - bntseq->anns[rid].offset +1);//Pos
            //ksprintf(s, "*,");//cigar not compute
            if(is_cigar){
                if(hit->is_gap){
                    char *cigar = calloc(256, sizeof(char)); 
                    int d = ed_diff_withcigar(index->mixRef->seq, hit->pos, query->l_seq+4, strand==0?query->seq:query->rseq, query->l_seq, hit->n_diff, cigar, 256, 1, COMPACT_CIGAR_STRING);
                    if(d != hit->n_diff) {
                        fprintf(stderr, "[%s]: diff error!\n", __func__); 
                        fprintf(stderr, "[%s]: hit->pos, %d != %d\n", query->name, d, hit->n_diff); 
                        exit(1);
                    }
                    ksprintf(s, "%s,", cigar);
                    free(cigar);
                } else{
                    ksprintf(s, "%dM,", query->l_seq);
                } 
            
            }else{
                ksprintf(s, "*,");//cigar not compute 
            }
            ksprintf(s, "%u", n);//NM
            ksprintf(s, ";"); 
            //last_pos = p;
        }
    }
   
    
}

#define MAX_RS 64
//#define __get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __get_ref(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __get_metaref(pac, l) (((pac)[(l)>>3]>>4*((l)%8))&15)
void sam_add_md_nm(kstring_t *s, index_t *index, query_t *q)
{
    if(q->pos == 0xFFFFFFFF) return;
    int i; 
    int nm = 0;
    uint32_t ref_pos = q->pos;
    uint8_t *ref_pac = index->pac;
    int l_seq = q->l_seq;
    const uint8_t *seq0 = q->strand==0?q->seq:q->rseq; 
    const uint8_t *seq = seq0+q->seq_start;
    int l_cigar = q->cigar->l; 
    const char *cigar = q->cigar->s; 
    int n_rs = 0;
    int rs[MAX_RS] = {0};
    ksprintf(s, "\tMD:Z:");
    int n_match = 0;
    while(strlen(cigar) != 0){
        long n = strtol(cigar, &cigar, 10); 
        char op = *cigar;
       
        switch(op){
            case 'M':
        
                for(i = 0; i < n; ++i){
                    //uint8_t bt = __get_pac(ref_pac, ref_pos);
                    assert(ref_pos < index->bntseq->l_pac);
                    uint8_t bt = __get_ref(ref_pac, ref_pos);
                    if(bt == *seq){
                        n_match += 1;
                    }else{//MISMATCH
                        bntseq_t *bntseq = index->bntseq;
                        uint8_t meta_bt = __get_metaref(index->mixRef->seq, ref_pos);
                        int rid;
                        bns_coor_pac2real(bntseq, ref_pos, q->l_seq, &rid); 
                        //fprintf(stderr, "snp_pos:%s, %u, %x, %x\n", bntseq->anns[rid].name, ref_pos-bntseq->anns[rid].offset+1, meta_bt, 1<<(*seq));
                        if((meta_bt & 1<<(*seq)) != 0) {
                            if(n_rs < MAX_RS){
                                rs[n_rs] = (seq-(seq0+q->seq_start))/sizeof(*seq);
                                n_rs += 1; 
                            }
                        }
                        nm += 1; 
                        if (n_match != 0) ksprintf(s, "%d", n_match); 
                        n_match = 0;
                        kputc("ACGTN"[bt], s); 
               
                    } 
                    ref_pos += 1;
                    seq += 1;
                }
                //if(n_match != 0) ksprintf(s, "%d", n_match);
                break;
            case 'I':
                nm+=n; 
                seq += n;
                break;
            case 'D':
                if(n_match != 0) ksprintf(s, "%d", n_match);
                n_match = 0;
                nm += n;
                kputc('^', s);
                for(i=0; i<n; ++i){
                    kputc("ACGTN"[__get_ref(ref_pac, ref_pos)],s); 
                    ref_pos += 1; 
                }  
                break;
            case 'S': 
                break;
        }
        cigar += 1;
    }
    if(n_match != 0) ksprintf(s, "%d", n_match);
    
    ksprintf(s, "\tNM:i:%u", nm);
    if(n_rs >0) {
        ksprintf(s, "\tXV:i:");
        for(i=0; i < n_rs; ++i){    
            if(i!=0) ksprintf(s, ",");
            ksprintf(s, "%d", rs[i]);
        }
    }
    return;
}
#define STRAND_FORWARD 0
#define STRAND_REV 1
void alnpe_sam(index_t *index, query_t *q, const aln_opt_t *opt)
{

    bntseq_t *bntseq = index->bntseq;
    int i, j;
    const uint8_t *seq, *rseq, *qual;
    const char *name;

    int rid[2] = {-1,  -1};
    uint32_t pos[2] = {0, 0};
    int is_map[2] = {0, 0};
    for(i =0; i <2; ++i) { 
        if(q[i].pos != 0xFFFFFFFF) {
            if(q[i].pos > index->bntseq->l_pac){
                fprintf(stderr, "[%s]: pos > reference len\n", __func__ );
                fprintf(stderr, "[%s]: pos = %u, reference len = %u\n", __func__, q[i].pos, index->bntseq->l_pac );
            }
            is_map[i] = 1; 
            bns_coor_pac2real(bntseq, q[i].pos, q[i].l_seq, &rid[i]); 
            pos[i] = q[i].pos-bntseq->anns[rid[i]].offset+1; 
        }
    } 
    int tlen = 0;
    if(is_map[0] &&is_map[1]){
        if(rid[0] != rid[1]) tlen = 0;
        else if(pos[0]<pos[1]) tlen = pos[1]+q[1].seq_end-q[1].seq_start+1 - pos[0];
        else tlen = pos[0]+q[0].seq_end-q[1].seq_start+1 - pos[1];
        if(tlen > opt->max_tlen || tlen < opt->min_tlen) tlen = 0;
    } 
    //const hit_t *hit;
    uint16_t flag;
    
    //kstring_t *s = calloc(1, sizeof(kstring_t));
    //s->m = 1024;
    //s->s = calloc(1024, 1); 
    

    for(i = 0; i < 2; ++i){
        kstring_t *s = q[i].sam;
        seq = q[i].seq;
        rseq = q[i].rseq;
        qual = q[i].qual;
  
        //Qname
        name = q[i].name;
        ksprintf(s, "%s\t", name);        
        //Flag 
        flag =0;
        flag |= SAM_FLAG_PAIRED;
        if(is_map[i] != 1) flag |= SAM_FLAG_UNMAPED;
        if(is_map[1-i] != 1) flag |= SAM_FLAG_MATE_UNMAPPED;
        if(q[i].strand == STRAND_REV) flag |= SAM_FLAG_REVERSE;
        if(q[1-i].strand == STRAND_REV) flag |= SAM_FLAG_MATE_REVERSE;
        if(tlen != 0) flag |= SAM_FLAG_PROPER_PAIRED;
        flag |= (i==0)?SAM_FLAG_MATE_READ1:SAM_FLAG_MATE_READ2;
        ksprintf(s, "%u\t", flag);
        //Rname Pos Mapq Cigar
        if(is_map[i]) {
            ksprintf(s, "%s\t", bntseq->anns[rid[i]].name);//Rname
            ksprintf(s, "%lu\t", pos[i]);//Pos
            ksprintf(s, "%u\t", q[i].mapq);//Mapq
            //Cigar
            if(q[i].seq_start != 0) ksprintf(s,"%dS",q[i].seq_start);
            ksprintf(s, "%s", q[i].cigar->s);
            if(q[i].seq_end != q[i].l_seq-1) ksprintf(s,"%dS",q[i].l_seq -q[i].seq_end-1);
#ifdef DEBUG
            if(check_cigar(q[i].cigar->s, q[i].cigar->l, q[i].seq_end-q[i].seq_start+1) == CIGAR_ERRO) {
                fprintf(stderr, "[Cigar error]: %s\n", q[i].name);
                fprintf(stderr, "[Cigar error]: %s\n", q[i].cigar->s);
                exit(1);
            }
#endif
            kputc('\t',s);
        } else{
            if(is_map[1-i]){
                ksprintf(s, "%s\t", bntseq->anns[rid[1-i]].name);//Rname
                ksprintf(s, "%lu\t", pos[1-i]);//Pos
                ksprintf(s, "255\t*\t");
            } else{
                ksprintf(s, "*\t0\t255\t*\t");                 
            }

       
        }
        //Rnext and Pnext
        if(is_map[1-i]){
            if(rid[i] == rid[1-i] || is_map[i]!=1) { ksprintf(s, "=\t");} 
            else{ ksprintf(s, "%s\t", bntseq->anns[rid[1-i]].name);}
            ksprintf(s, "%lu\t", pos[1-i]);
        } else{
            ksprintf(s,"*\t0\t");
        }

        //Tlen
        if(tlen != 0){
           if(q[i].pos >= q[1-i].pos) ksprintf(s, "-%d\t", tlen); 
           else ksprintf(s, "%d\t", tlen); 
        }  else {ksprintf(s, "0\t");} 
        //Seq and Qual
        if(q[i].strand == STRAND_REV) {
            for(j = 0; j < q[i].l_seq; ++j) kputc("ACGTN"[rseq[j]],s);
            kputc('\t', s);
            if(strlen((const char*)q[i].qual)) {
                for(j = q[i].l_seq -1; j >= 0; --j) kputc(qual[j],s);
            } else{
                kputc('*',s);
            }
        } else{
            for(j = 0; j < q[i].l_seq; ++j) kputc("ACGTN"[seq[j]],s);
            kputc('\t',s);
            if(strlen((const char*)qual)) {
                ksprintf(s, "%s", qual);
            } else{
                kputc('*',s);
            }
        }
        //Optional filed  
        sam_add_xa(s, index, q+i, opt->print_xa_cigar);
        if(opt->print_nm_md) sam_add_md_nm(s, index, q+i); 
        if(opt->rg_id != NULL) ksprintf(s, "\tRG:Z:%s", opt->rg_id);
        kputc('\n',s);
    }

    //if(s->s != NULL)    printf("%s",s->s);    
    //fflush(stdout);
    //free(s->s);
    //free(s);
}
