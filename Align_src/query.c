/*
 * =====================================================================================
 *
 *       Filename:  query.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/28/2012 10:11:20 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kstring.h"
#include "query.h"
#include "editdistance.h"
#include "aln.h"

#define UINT_MAX 0xFFFFFFFF
#define SCORE_MATCH 1
#define SCORE_MISMATCH -2
#define SCORE_GAP -5
#define CIGAR_MATCH 0
#define CIGAR_INS 1
#define CIGAR_DEL 2


extern unsigned char nst_nt4_table[256];

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  query_seq_reverse
 *  Description:  from bwa seqio.c
 *                if comp == 1 return reverse complement
 *                if comp == 0 return mirror seq 
: * =====================================================================================
 */
void query_seq_reverse(int len, uint8_t *seq, int is_comp)
{
	
    int i;
	if (is_comp) {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			if (tmp < 4) tmp = 3 - tmp;
			seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
			seq[i] = tmp;
		}
		if (len&1) seq[i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	} else {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			seq[len-1-i] = seq[i]; seq[i] = tmp;
		}
	}
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  query_init
 *  Description:  query init func
 * =====================================================================================
 */


void query_destroy(query_t *query)
{
    //fprintf(stderr, "%s\n", __func__);

    if(query->name != NULL) free(query->name);
    if(query->comment != NULL) free(query->comment);
    if(query->seq != NULL) free(query->seq);
    if(query->rseq != NULL) free(query->rseq);
    if(query->qual != NULL) free(query->qual);
    if(query->cigar != NULL) {
        free(query->cigar->s);
        free(query->cigar);
    }
    if(query->sam != NULL ) {
        free(query->sam->s);
        free(query->sam);
    }
    kv_destroy(query->hits[0]);
    kv_destroy(query->hits[1]);


    //free(query->hits[0].a);
    //free(query->hits[1].a);
    //kv_destroy(query->ordered_index);
    //kv_destroy(query->unique_index);
}
queryio_t *query_open(const char *fn_fa)
{
    //fprintf(stderr, "%s\n", __func__);

    gzFile fp = 0;
    kseq_t *ks;
    queryio_t *qs = NULL;

    qs = calloc(1, sizeof(queryio_t));
    if(qs == NULL){
    	fprintf(stderr, "[query_init]: allocate mem fail!\n");
    	exit(EXIT_FAILURE);
    }
    fp = gzopen(fn_fa, "r");
    if(fp == Z_NULL){
        fprintf(stderr, "[query_open]: file %s open fail!\n", fn_fa);    
        exit(1);
    }
    ks = kseq_init(fp);
    qs->kseq = ks;
    qs->fp = fp;

    return qs;
}
void query_close(queryio_t *qs)
{    
    //fprintf(stderr, "%s\n", __func__);
    
	kseq_destroy(qs->kseq);
    gzclose(qs->fp);
	free(qs);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  query_read
 *  Description:  read query 
 * =====================================================================================
 */

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}
#define INIT_HIT_NUM 32
int query_read_seq(queryio_t *qs, query_t *query)
{
	int i;
    kseq_t *kseq = qs->kseq;      
    int l_seq = kseq_read(kseq);   
    
    //fprintf(stderr, "[query_read_seq]:%s\n%s\n%s\n",kseq->name.s, kseq->seq.s, kseq->qual.s);
    //fprintf(stderr, "%s\n", __func__);
    
    if(l_seq <= 0) return 0;    
    //read name and comment
    if(kseq->name.s != NULL) {
        trim_readno(&kseq->name);
        query->name = strdup(kseq->name.s);}      
    if(kseq->comment.s != NULL) query->comment = strdup(kseq->comment.s);      

    //read seq and generate reverse seq
    query->l_seq = l_seq;
    query->seq = calloc(l_seq, 1);
    if(  query->seq  == NULL  ){
        fprintf(stderr, "[query_read]: realloc mem fail, mem leak!!\n");
        exit(1);
    }
    query->seq_start = 0;
    query->seq_end = l_seq-1;
    query->rseq = calloc(l_seq, 1);
    if(  query->rseq  == NULL  ){
        fprintf(stderr, "[query_read]: realloc mem fail, mem leak!!\n");
        exit(1);
    }
    query->n_ambiguous = 0;
    for(i = 0; i < l_seq; ++i){
        unsigned nt = nst_nt4_table[(uint8_t)kseq->seq.s[i]];
        if(nt>3) { ++query->n_ambiguous;} 
        query->seq[i] = nt;        
    }
    memcpy(query->rseq, query->seq, l_seq);
    query_seq_reverse(l_seq, query->rseq, 1);
    if(kseq->qual.l > 0){
        /*
        query->qual = calloc(kseq->qual.l, 1);
        if(  query->qual  == NULL  ){
            fprintf(stderr, "[query_read]: realloc mem fail, mem leak!!\n");
            exit(1);
        }
        strncpy(query->qual, kseq->qual.s, kseq->qual.l);
        */
        query->qual = (uint8_t *)strdup(kseq->qual.s);
        //convert to sam base qual
        /*
        uint8_t *q = query->qual;
        for(i = 0; i < kseq->qual.l; ++i){
            q[i] += 33;
        }*/
    }
    query->is_gap = -1;
    query->pos = 0xFFFFFFFF;
    query->n_diff = -1;
    query->strand = 3;
    query->b0 = -1;      
    query->b1 = -1;
    query->cigar = calloc(1, sizeof(kstring_t));
    query->cigar->s = calloc(128, 1);
    query->cigar->m = 128;
    if(query->cigar == NULL) {
        fprintf(stderr, "[query_read_seq]: query->cigar allocate mem fail!\n");
        exit(1);
    }
    query->sam = calloc(1, sizeof(kstring_t));
    query->sam->s = calloc(128, 1);
    query->sam->m = 128;


    //init candidate vec
    //int n_c, n_for_r, n_back_r;
    
    //n_c = (l_seq+24)/25;
    //n_for_r = (l_seq+24)/25*25;
    //n_back_r = (l_seq+24)/25*25;
    
    //kv_resize(sa_interval_t, query->candidate.C,  n_c);
    //kv_resize(sa_interval_t, query->candidate.for_R,  n_for_r);
    //kv_resize(sa_interval_t, query->candidate.back_R,  n_back_r);
    //kv_resize(sa_interval_t, query->rcandidate.C,  n_c);
    //kv_resize(sa_interval_t, query->rcandidate.for_R,  n_for_r);
    //kv_resize(sa_interval_t, query->rcandidate.back_R,  n_back_r);

    //init hit
    kv_resize(hit_t, query->hits[0], INIT_HIT_NUM);
    kv_resize(hit_t, query->hits[1], INIT_HIT_NUM);
    //kv_resize(uint32_t, query->ordered_index, INIT_HIT_NUM);
    //kv_resize(uint32_t, query->unique_index, INIT_HIT_NUM);
    return l_seq;
}		/* -----  end of function query_read  ----- */
int query_read_multiSeqs(queryio_t *qs, int n_seq, query_t *multiSeqs)
{
    int i=0;
    //fprintf(stderr, "%s\n", __func__);

    while(i < n_seq && query_read_seq(qs, multiSeqs+i) >0  ){
        ++i;
    }

    return i;
}

int query_read_multiPairedSeqs(queryio_t *qs[], int n_seq, query_t *multiSeqs)
{
    int i,j;
    i = 0, j = 1; 
    while(i < n_seq  && query_read_seq(qs[0], multiSeqs +i) >0 ){
        i += 2;
    }
    while(j < n_seq && query_read_seq(qs[1], multiSeqs+j) >0  ){
        j += 2;
    }
    if(i != j-1) {
        fprintf(stderr, "[%s]: the number of paired seqs are not same!\n", __func__);
        exit(1);
    }
    return i;

}
//[to do]: change
uint32_t gen_mapq(uint32_t b0, uint32_t b1)
{
    //uint32_t mapq = -4.343 * log(1-(double)abs(b0 - b1)/(double)b0);
    //mapq = (uint32_t)(mapq+4.99);
    if(b0 == 0) return 0;
    double a = 255.0;
    uint32_t mapq = a*((double)abs(b0-b1)/(double)b0);

     
    mapq = mapq <254?mapq:254;
    return mapq; 
}
void query_gen_cigar(uint32_t l_ref, const uint32_t *mixRef, query_t *query)
{
    query->seq_start = 0;
    query->seq_end = query->l_seq-1;
    if(query->pos != 0xFFFFFFFF){
        if(query->is_gap){//have gap
           int ret = ed_diff_withcigar(mixRef, query->pos, query->l_seq+4, query->strand==0?query->seq:query->rseq, query->l_seq, query->n_diff, query->cigar->s, query->cigar->m, 1, COMPACT_CIGAR_STRING);
            if(ret == -1) fprintf(stderr, "[%s]:Erro:%s\n", __func__, query->name);
            //fprintf(stderr, "[cigar]: %s\n", query->cigar->s);           
        } else ksprintf(query->cigar, "%dM", query->l_seq);
    
    
    } //else query->cigar = NULL;
   
}
void query_set_hits(query_t *query, int max_hits, hits_t *hits0, hits_t *hits1)
{
    size_t i, j;
    uint32_t primary_pos = query->pos;
    int tot_hits = 0; 
    query->b0 = query->n_diff;
    query->b1 = 100000;
    for(i = 0; i < 2; ++i){
        hits_t *hits = i==0?hits0:hits1; 
        hit_t *a = hits->a;
        uint32_t last_pos = (uint32_t)-1; 
        for(j = 0; j < hits->n; ++j){

            if(a[j].strand != i){
                fprintf(stderr, "[%s]: %s hits[%lu] strand %d\n", __func__, query->name, i, a[j].strand);
            
            }

            uint32_t pos = a[j].pos;
            if (pos == last_pos || pos == primary_pos) continue;
            if(a->n_diff <= query->n_diff){ 
                if(a->n_diff <= query->b1) query->b1 = a->n_diff;
                kv_push(hit_t, query->hits[i], a[j]); 
                tot_hits += 1;    
            }
            if(tot_hits == max_hits){ 
                goto end;
            }
        }  

    }
end:
    {
        query->mapq = gen_mapq(query->b0, query->b1);   
    
    }
}
