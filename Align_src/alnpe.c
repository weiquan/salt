/*
 * =====================================================================================
 *
 *       Filename:  alnpe.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/01/2013 02:09:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include "khash.h"
#include "aln.h"
#include "sam.h"
#include "ssw.h"
#include "utils.h"


//MACRO
#define POS_UNMAPPED 0xFFFFFFFF
#define SW_PAIRED 1
#define SW_UNPAIRED 0
#define PAIRED_ALNED 1
#define PAIRED_UNALNED 0

#define LOW_BOUNDARY -1
#define IN_RANGE 0 
#define HIGH_BOUNDARY 1

#define STRAND_CORRECT 0
#define STRAND_WRONG 1
#define __distance(a, b)(a<b?b-a:a-b)
#define __min(a,b) ((a)<(b)?(a):(b))
#define __set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))

const int GAP_OP = 5;
const int GAP_EX = 1;
int8_t score_mat[25] = { 1, -3, -3, -3,-1,
                                -3,  1, -3, -3,-1, 
                                -3, -3, 1, -3,-1, 
                                -3, -3, -3,  1,-1,
                                -1,  -1,  -1,  -1,  -1};
                                    //  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
int8_t score_mat2[256] = {       -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                /* 1 */ -3,  1, -3,  1, -3,  1, -3,  1, -3,  1, -3,  1, -3,  1, -3,  1,
                                /* 2 */ -3, -3,  1,  1, -3, -3,  1,  1, -3, -3,  1,  1, -3, -3,  1,  1,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                /* 4 */ -3, -3, -3, -3,  1,  1,  1,  1, -3, -3, -3, -3,  1,  1,  1,  1, 
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                /* 8 */ -3, -3, -3, -3, -3, -3, -3, -3,  1,  1,  1,  1,  1,  1,  1,  1, 
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
                                        -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3};


static inline int CHECK_IN_RANGE(uint32_t a, uint32_t b, uint32_t small, uint32_t large){
    uint32_t r = __distance(a, b);
    if(a>b || r < small) return LOW_BOUNDARY;
    else if(r>large) return HIGH_BOUNDARY;
    else return IN_RANGE;
}
static inline int CHECK_STRAND(uint32_t pos0, uint32_t pos1, int strand0, int strand1){
    if(pos0 < pos1){
        if(strand0 == STRAND_FORWARD && strand1 == STRAND_BACKWARD){ return STRAND_CORRECT;}
    } else{
        if(strand1 == STRAND_FORWARD && strand0 == STRAND_BACKWARD){ return STRAND_CORRECT;}
    }
    
    return STRAND_WRONG;

}


int pairing2(index_t *index, query_t *q0, query_t *q1, const aln_opt_t *aln_opt)
{
#define __swap_hit(q, hit) do{\
            SWAP((q)->pos, (hit)->pos);\
            SWAP((q)->strand, (hit)->strand);\
            SWAP((q)->n_diff, (hit)->n_diff);\
            SWAP((q)->is_gap, (hit)->is_gap);\
        } while(0)

    uint32_t l2 = q0->l_seq+q1->l_seq; 
    uint32_t min_isize = aln_opt->min_tlen>l2?aln_opt->min_tlen-l2:0;
    uint32_t max_isize = aln_opt->max_tlen>l2?aln_opt->max_tlen-l2:0;
    uint32_t min_erros = (uint32_t)-1;
    unsigned int l_pac = index->bntseq->l_pac; 
    int range;
    if(q0->strand ==STRAND_FORWARD && q1->strand == STRAND_BACKWARD && q0->pos < q1->pos) {
        range = CHECK_IN_RANGE(q0->pos+q0->l_seq, q1->pos, min_isize, max_isize);
        if(range == IN_RANGE){
            min_erros = q0->n_diff + q1->n_diff;
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
            return PAIRED_ALNED;
        }
    
    } else if( q1->strand ==STRAND_FORWARD && q0->strand == STRAND_BACKWARD && q1->pos < q0->pos){
        range = CHECK_IN_RANGE(q1->pos+q1->l_seq, q0->pos, min_isize, max_isize);
        if(range == IN_RANGE){
            min_erros = q0->n_diff+q1->n_diff; 
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
            return PAIRED_ALNED;
        }
    
    }
    //[toFix]: should check q0 primary pos with all q1 alternative pos
    //qo as forward, q1 as backward_hit
    
    hit_t b0, b1;
    uint32_t n_forward = q0->hits[0].n;
    uint32_t n_backward = q1->hits[1].n;
    if(n_forward>0 && n_backward >0){
        hit_t *forward_hits = q0->hits[0].a; 
        hit_t *backward_hits = q1->hits[1].a; 
        uint32_t l0 = q0->l_seq; 
        uint32_t l1 = q1->l_seq; 
        uint32_t i = 0, j = 0, jj;
        for(i=0; i< n_forward; ++i){
            uint32_t pos0 = forward_hits[i].pos;
            for(jj=j; jj < n_backward; ++jj){
                uint32_t pos1 = backward_hits[jj].pos;
                range = CHECK_IN_RANGE(pos0+l0, pos1, min_isize, max_isize);   
                if(range == IN_RANGE){
                    if (forward_hits[i].n_diff+backward_hits[jj].n_diff < min_erros){
                        min_erros = forward_hits[i].n_diff+backward_hits[jj].n_diff;
                        //__swap_hit(q0, &forward_hits[i]);
                        //__swap_hit(q1, &backward_hits[jj]);
                        b0 = forward_hits[i];
                        b1 = backward_hits[jj];
                    }
                } else if(range == LOW_BOUNDARY) {
                     j == jj;
                } else if(range == HIGH_BOUNDARY) {
                    break; 
                }
            }
        }
    } 
    //q1 as forward ,q0 as forward 
    n_forward = q1->hits[0].n;
    n_backward = q0->hits[1].n;
    if(n_forward>0 && n_backward >0){
        hit_t *forward_hits = q1->hits[0].a; 
        hit_t *backward_hits = q0->hits[1].a; 
        uint32_t l0 = q1->l_seq; 
        uint32_t l1 = q0->l_seq; 
        uint32_t i = 0, j = 0, jj;
        for(i=0; i< n_forward; ++i){
            uint32_t pos0 = forward_hits[i].pos;
            for(jj=j; jj < n_backward; ++jj){
                uint32_t pos1 = backward_hits[jj].pos;
                range = CHECK_IN_RANGE(pos0+l0, pos1, min_isize, max_isize);   
                if(range == IN_RANGE){
                    if (forward_hits[i].n_diff+backward_hits[jj].n_diff < min_erros){
                        min_erros = forward_hits[i].n_diff+backward_hits[jj].n_diff;
                        //__swap_hit(q1, &forward_hits[i]);
                        //__swap_hit(q0, &backward_hits[jj]);
                        b1=forward_hits[i];
                        b0=backward_hits[jj];
                    }
                } else if(range == LOW_BOUNDARY) {
                     j == jj;
                } else if(range == HIGH_BOUNDARY) {
                    break; 
                }
            }
        }

    }
/*proper paired */
    if(min_erros != (uint32_t)-1){
        q0->pos = b0.pos; q0->strand = b0.strand; q0->n_diff = b0.n_diff;q0->is_gap=b0.is_gap;
        q1->pos = b1.pos; q1->strand = b1.strand; q1->n_diff = b1.n_diff;q1->is_gap=b1.is_gap;
        
        query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);
        query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
        return PAIRED_ALNED; 
    }
    
    
    
/*Singletons mapped*/
    
    uint32_t start, end;
    uint32_t max_fragment_len = max_isize+q0->l_seq+q1->l_seq;
    
    uint32_t n_hits0 = q0->hits[0].n+q0->hits[1].n;
    uint32_t n_hits1 = q1->hits[0].n+q1->hits[1].n;
    //use q0 as anchor

    if(q0->strand == STRAND_FORWARD){
        start = q0->pos+min_isize+q0->l_seq;
        end = q0->pos+max_isize+q0->l_seq+q1->l_seq;
        end = end >= l_pac?l_pac:end; 

        //if(snpaln_sw(index, &start, &end, q1, STRAND_BACKWARD, aln_opt) == SW_PAIRED){
        if(snpaln_sw_snpaware(index, &start, &end, q1, STRAND_BACKWARD, aln_opt) == SW_PAIRED){
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);               
            return PAIRED_ALNED;
       }
    } else{
        start = q0->pos>max_isize+q1->l_seq?q0->pos-max_isize-q1->l_seq:0;
        end = q0->pos>min_isize?q0->pos-min_isize:0; end = end >= l_pac?l_pac:end; 

        //if(snpaln_sw(index, &start, &end, q1, STRAND_FORWARD, aln_opt) == SW_PAIRED){
        if(snpaln_sw_snpaware(index, &start, &end, q1, STRAND_FORWARD, aln_opt) == SW_PAIRED){
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);                
            return PAIRED_ALNED;
        }
    }

    //use q1 as anchor        
    if(q1->strand == STRAND_FORWARD){
        start = q1->pos+min_isize+q1->l_seq;
        end = q1->pos+max_isize+q1->l_seq+q0->l_seq; end = end >= l_pac?l_pac:end; 
        //if(snpaln_sw(index, &start, &end, q0, STRAND_BACKWARD, aln_opt) == SW_PAIRED){
        if(snpaln_sw_snpaware(index, &start, &end, q0, STRAND_BACKWARD, aln_opt) == SW_PAIRED){
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
            return PAIRED_ALNED;
        }
    } else{
        start = q1->pos>max_isize+q0->l_seq?q1->pos-max_isize-q0->l_seq:0;
        end = q1->pos> min_isize?q1->pos-min_isize:0;        
        end = end >= l_pac?l_pac:end; 
        //if(snpaln_sw(index, &start, &end, q0, STRAND_FORWARD, aln_opt) == SW_PAIRED){
        if(snpaln_sw_snpaware(index, &start, &end, q0, STRAND_FORWARD, aln_opt) == SW_PAIRED){
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
            return PAIRED_ALNED;
        }
    }
    if(q0->pos != POS_UNMAPPED)  query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);
    if(q1->pos != POS_UNMAPPED)  query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
    return PAIRED_UNALNED;


}

#define __get_mixref(pac, l) (((pac)[(l)>>3]>>4*((l)%8))&15)
int snpaln_sw_snpaware(index_t *index, uint32_t *start_pos, uint32_t *end_pos, query_t *q, int strand, const aln_opt_t *opt)
{
    int i;
    int ret = SW_UNPAIRED; 
    
    if((*start_pos)>=index->bntseq->l_pac){ 
        fprintf(stderr, "error %s\n", q->name);
        exit(1);
    }
    const uint32_t start = *start_pos, end = *end_pos;
    
    uint8_t *seq = strand?q->rseq:q->seq;
    int32_t l_seq = q->l_seq; 
   
    
    int32_t l_ref = end -start +1; 
    uint8_t *ref = (uint8_t *)calloc(l_ref, 1);
    uint8_t *read = (uint8_t *)calloc(l_seq, 1);
    
    for(i = 0; i != end - start +1; ++i){
        ref[i] = __get_mixref(index->mixRef->seq, start+i);
    }        
    for(i = 0; i < l_seq; ++i) read[i] = 1<<seq[i]; 
    s_profile *p = ssw_init((int8_t *)read, l_seq, score_mat2, 16, 1); 
    //fprintf(stderr, "Ref:\n");
    //for(i = 0; i != end - start +1; ++i) fprintf(stderr, "%c", "ACGT"[ref[i]]);
    //fprintf(stderr, "\nSeq:\n");
    //for(i = 0; i != l_seq; ++i) fprintf(stderr, "%c", "ACGT"[seq[i]]);
    uint8_t flag = 2;
    //int filter = 80;
    //int filtered = 0;
    int maskLen = q->l_seq/2;
    s_align *result = ssw_align(p, (int8_t *)ref, l_ref, opt->gap_op, opt->gap_ex, flag, opt->filters, opt->filterd, maskLen);//note: filterd no effect

    if(result->score1 >= opt->filters && result->read_end1-result->read_begin1+1 >=opt->filterd){
        q->b0 = result->score1;
        q->b1 = result->score2;
        q->mapq = gen_mapq(q->b0, q->b1);
        q->pos = result->ref_begin1+start;
        q->strand = strand;
        q->seq_start = result->read_begin1;
        q->seq_end = result->read_end1;
        q->cigar->l = 0;//clean cigar data 
        int j;
        for(j = 0; j < result->cigarLen; ++j) {
            uint32_t cigar = result->cigar[j];    
            ksprintf(q->cigar, "%u%c",cigar>>4,"MID"[cigar&15]);//Cigar
        }
        
        *start_pos = result->ref_begin1+start;
        *end_pos = result->ref_end1 +start;

        q->seq_start = result->read_begin1;
        q->seq_end = result->read_end1;
        ret = SW_PAIRED; 

    } else{}
    //fprintf(stderr, "\nRef Pos:[%d, %d]\n", result->ref_begin1, result->ref_end1);
    //fprintf(stderr, "Seq Pos:[%d, %d]\n", result->read_begin1, result->read_end1);


    free(ref);
    free(read);
    align_destroy(result);
    init_destroy(p);
    return ret;
}


int snpaln_sw(index_t *index, uint32_t *start_pos, uint32_t *end_pos, query_t *q, int strand, const aln_opt_t *opt)
{
    int i;
    int ret = SW_UNPAIRED; 
    
    if((*start_pos)>=index->bntseq->l_pac){ 
        fprintf(stderr, "error %s\n", q->name);
        exit(1);
    }
    const uint32_t start = *start_pos, end = *end_pos;

    
    uint8_t *seq = strand?q->rseq:q->seq;
    int32_t l_seq = q->l_seq; 
    s_profile *p = ssw_init((int8_t *)seq, l_seq, score_mat, 5, 1);
    
    int32_t l_ref = end -start +1; 
    uint8_t *ref = calloc(l_ref, 1);
    
    for(i = 0; i != end - start +1; ++i){
        ref[i] = __get_pac(index->pac, start+i);
    }        
    //fprintf(stderr, "Ref:\n");
    //for(i = 0; i != end - start +1; ++i) fprintf(stderr, "%c", "ACGT"[ref[i]]);
    //fprintf(stderr, "\nSeq:\n");
    //for(i = 0; i != l_seq; ++i) fprintf(stderr, "%c", "ACGT"[seq[i]]);
    uint8_t flag = 2;
    //int filter = 80;
    //int filtered = 0;
    int maskLen = q->l_seq/2;
    s_align *result = ssw_align(p, (int8_t *)ref, l_ref, opt->gap_op, opt->gap_ex, flag, opt->filters, opt->filterd, maskLen);//note: filterd no effect

    if(result->score1 >= opt->filters && result->read_end1-result->read_begin1+1 >=opt->filterd){
        q->b0 = result->score1;
        q->b1 = result->score2;
        q->mapq = gen_mapq(q->b0, q->b1);
        q->pos = result->ref_begin1+start;
        q->strand = strand;
        q->seq_start = result->read_begin1;
        q->seq_end = result->read_end1;
        q->cigar->l = 0;//clean cigar data 
        int j;
        for(j = 0; j < result->cigarLen; ++j) {
            uint32_t cigar = result->cigar[j];    
            ksprintf(q->cigar, "%u%c",cigar>>4,"MID"[cigar&15]);//Cigar
        }
        
        *start_pos = result->ref_begin1+start;
        *end_pos = result->ref_end1 +start;

        q->seq_start = result->read_begin1;
        q->seq_end = result->read_end1;
        ret = SW_PAIRED; 

    } else{}
    //fprintf(stderr, "\nRef Pos:[%d, %d]\n", result->ref_begin1, result->ref_end1);
    //fprintf(stderr, "Seq Pos:[%d, %d]\n", result->read_begin1, result->read_end1);


    free(ref);
    align_destroy(result);
    init_destroy(p);
    return ret;
}

int pairing_singleton(index_t *index, query_t *q0, query_t *q1, aln_opt_t *aln_opt)
{
    int n0 = q0->hits[0].n+q0->hits[1].n;
    int n1 = q1->hits[0].n+q1->hits[1].n;
    uint32_t l2 = q0->l_seq+q1->l_seq; 
    uint32_t min_isize = aln_opt->min_tlen>l2?aln_opt->min_tlen-l2:0;
    uint32_t max_isize = aln_opt->max_tlen>l2?aln_opt->max_tlen-l2:0;
  
   
    
    if(q0->pos == 0xFFFFFFFF && q1->pos == 0xFFFFFFFF){
       return PAIRED_UNALNED;
    }
    /*
    if(n0 >0 && n1 >0){
        return pairing0(q0, q1, p, opt);
    }*/
    
    uint32_t start, end;
    if(q0->pos != 0xFFFFFFFF){
        
        if(q0->strand == STRAND_FORWARD){
            start = q0->pos+min_isize+q0->l_seq;
            start = __min(start, index->bntseq->l_pac-1);
            end = q0->pos+max_isize+q0->l_seq+q1->l_seq;
            end = __min(end, index->bntseq->l_pac-1);
            //fprintf(stderr, "[SW]:%s\n", q1->name);
            if (snpaln_sw(index, &start, &end, q1, STRAND_BACKWARD, aln_opt) == SW_PAIRED){
                query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);        
                //p->isize = p->pos[1]==0xFFFFFFFF ?0:end - q0->pos +1;
                return PAIRED_ALNED;
            }
        } else{
            //start = q0->pos-max_isize-q1->l_seq;
            start = q0->pos>max_isize+q1->l_seq?q0->pos-max_isize-q1->l_seq:0;
            start = __min(start, index->bntseq->l_pac-1);
            //end = q0->pos- min_isize;
            end = q0->pos> min_isize?q0->pos-min_isize:0;
            end = __min(end, index->bntseq->l_pac-1);
            //fprintf(stderr, "[SW]:%s\n", q1->name);
            if(snpaln_sw(index, &start, &end, q1, STRAND_FORWARD, aln_opt) == SW_PAIRED){
                query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);        
                //p->isize = p->pos[1]==0xFFFFFFFF ?0:q0->pos+q0->l_seq - start;
                return PAIRED_ALNED;
            }
        }
    } 
    if(q1->pos != 0xFFFFFFFF){

        if(q1->strand == STRAND_FORWARD){
            start = q1->pos+min_isize+q1->l_seq;
            start = __min(start, index->bntseq->l_pac-1);
            end = q1->pos+max_isize+q1->l_seq+q0->l_seq;
            end = __min(end, index->bntseq->l_pac-1);
            //fprintf(stderr, "[SW]:%s\n", q0->name);
            if(snpaln_sw(index, &start, &end, q0, STRAND_BACKWARD, aln_opt) == SW_PAIRED){
                query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
                return PAIRED_ALNED;
            }
        } else{
            //start = q1->pos-max_isize-q0->l_seq;
            start = q1->pos>max_isize+q0->l_seq?q1->pos-max_isize-q0->l_seq:0;
            start = __min(start, index->bntseq->l_pac-1);
            //end = q1->pos- min_isize;
            end = q1->pos> min_isize?q1->pos-min_isize:0;
            end = __min(end, index->bntseq->l_pac-1);
            //fprintf(stderr, "[SW]:%s\n", q0->name);
            if(snpaln_sw(index, &start, &end, q0, STRAND_FORWARD, aln_opt) == SW_PAIRED){
                query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
                return PAIRED_ALNED;
            }
        }
    } 
    if(q0->pos != POS_UNMAPPED)  query_gen_cigar(index->mixRef->l, index->mixRef->seq, q0);
    if(q1->pos != POS_UNMAPPED)  query_gen_cigar(index->mixRef->l, index->mixRef->seq, q1);
    return PAIRED_UNALNED;

    
    //fprintf(stderr, "There must be something wrong with the func %s\n", __func__);
  
   




}
#define MAX_N_PERSEQ 5
void alnpe_core1(int tid, index_t *index, int n_query, query_t *multi_query, aln_opt_t *aln_opt, aux_t* aux[2])
{
    int i;
    for(i=0; i < n_query/2; ++i){ 
#ifdef HAVE_THREAD
        if (i % aln_opt->n_threads != tid) continue;
#endif
        //fprintf(stderr, "[tid:%d]run seq %d", tid, i); 
        //fprintf(stderr, "aln_opt:%d", aln_opt->n_threads); 
        
        int j;
        for(j=0; j< 2; ++j){
            query_t *query = multi_query+i*2+j;
            if(query->n_ambiguous > MAX_N_PERSEQ) continue;
            if(query->l_seq -aln_opt->l_seed+1 > aux[0]->n_sai_range){
                int n_sai_range = query->l_seq-aln_opt->l_seed+1; 
                aux_resize(aux[0], n_sai_range); 
                aux_resize(aux[1], n_sai_range); 
            }
            aux_reset(aux[0]); 
            aux_reset(aux[1]); 
            
            if(aln_opt->l_overlap > 0) alnse_overlap(index, query, aln_opt, aux);
            else alnse_nonoverlap(index, query, aln_opt, aux); 
        }
        query_t *q0 = multi_query+i*2;
        query_t *q1 = multi_query+i*2+1;
        
        //fprintf(stderr, "pos = %u %u\n", q0->pos, q1->pos);
        if(q0->pos != 0xFFFFFFFF  && q1->pos != 0xFFFFFFFF){
            pairing2(index, q0, q1, aln_opt);
        } else if(q0->pos != 0xFFFFFFFF || q1->pos != 0xFFFFFFFF){
            pairing_singleton(index, q0, q1, aln_opt);
        } else{
            
        }
        //fprintf(stderr, "pos = %u %u\n", q0->pos, q1->pos);
        alnpe_sam(index, multi_query+i*2, aln_opt);
    }
}
#ifdef HAVE_THREAD
static void *alnpe_worker(void *data)
{
    thread_aux_t *d = (thread_aux_t*)data;
    alnpe_core1(d->tid, d->index, d->n_query, d->query, d->aln_opt, d->aux); 
    return 0;
}
#endif
int alnpe_core(const opt_t *opt)
{

    int i;
    double t = 0.0;
    index_t *index;
    
    aln_opt_t *aln_opt = aln_opt_init(opt);   
    fprintf(stderr, "[alnpe_core]:  Start paired end alignment!\n");
    fprintf(stderr, "[alnpe_core]:  Reload Index...\n"); t = clock();   
    index = alnpe_index_reload(opt->fn_index);
    fprintf(stderr, "%.2f sec\n", (float)(clock()-t)/CLOCKS_PER_SEC);
    
    queryio_t *qs[2];
    qs[0] = query_open(opt->fn_read1);
    qs[1] = query_open(opt->fn_read2);

    aux_t *aux[2];
#ifdef HAVE_THREAD
    pthread_attr_t attr;
    pthread_t *tid; 
    thread_aux_t *thread_data;
    if(opt->n_threads<=1){
        aux[0] = aux_init(opt->l_read, opt->l_seed);        
        aux[1] = aux_init(opt->l_read, opt->l_seed); 
    }else{
    
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
        thread_data = calloc(opt->n_threads, sizeof(thread_aux_t)); 
        int j;
        for(j=0; j < opt->n_threads;++j){
            thread_data[j].tid = j; thread_data[j].index = index;
            thread_data[j].aln_opt= aln_opt;
            thread_data[j].aux[0] = aux_init(opt->l_read, opt->l_seed);        
            thread_data[j].aux[1] = aux_init(opt->l_read, opt->l_seed);        
        } 
    }       
#else
    aux[0] = aux_init(opt->l_read, opt->l_seed);        
    aux[1] = aux_init(opt->l_read, opt->l_seed); 
#endif



   

      
    aln_samhead(opt, index->bntseq);
    
    int n, tot =0, max_lseq = -1;
    query_t *multi_seqs = calloc(N_SEQS, sizeof(query_t));
    while( (n = query_read_multiPairedSeqs(qs, N_SEQS, multi_seqs)) > 0){
        //fprintf(stderr, "read %d paired seqs!\n", n);

        if(opt->max_tlen == 0){
            fprintf(stderr, "infer isize func haven't been implemented\n");
            break;
        }

        /****************Aligning Reads************************************/ 
        t = realtime();       


#ifdef HAVE_THREAD
        if(opt->n_threads <=1){
            alnpe_core1(0, index, n, multi_seqs, aln_opt, aux); 
        }else{
            int j;
            for (j = 0; j < opt->n_threads; ++j) {
                thread_data[j].query = multi_seqs;
                thread_data[j].n_query = n; 
                pthread_create(&tid[j], &attr, alnpe_worker, thread_data + j);
            }
            for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0); 
        }
#else            
        alnpe_core1(0, index, n, multi_seqs, aln_opt, aux); 
#endif       
 
        /****************Output SAM ******************************************/ 
        tot += n;
        if(tot > 0&&tot % 100 == 0 ) fprintf(stderr, "alned %d reads!\n", tot);
        t = realtime();
	    for(i =0; i < n/2; ++i){
            query_t *q[2];
            q[0] = multi_seqs+i*2;
            q[1] = multi_seqs+i*2+1;
            //alnpe_sam(index, q, opt);
            printf("%s\n", q[0]->sam->s);
            printf("%s\n", q[1]->sam->s);
            query_destroy(q[0]);
            query_destroy(q[1]);
        }
        memset(multi_seqs, '\0', sizeof(query_t)*N_SEQS);



    	//fprintf(stderr, "print SAM escaped %f\n", (realtime()-t)); 
    }

#ifdef HAVE_THREAD
    if(opt->n_threads<=1){
        aux_destroy(aux[0]);
        aux_destroy(aux[1]);
    }else{
        int j;
        for(j=0; j < opt->n_threads;++j){
            aux_destroy(thread_data[j].aux[0]);
            aux_destroy(thread_data[j].aux[1]);
        } 

        free(tid);
        free(thread_data);
    

    }
#else
    aux_destroy(aux[0]);
    aux_destroy(aux[1]);
#endif    
    

    query_close(qs[0]);
    query_close(qs[1]);
    free(multi_seqs);

    alnpe_index_destroy(index);
    aln_opt_destroy(aln_opt);
    return 0;
}
