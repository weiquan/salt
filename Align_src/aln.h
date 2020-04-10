#ifndef __ALN_H
#define __ALN_H
/*
 * =====================================================================================
 *
 *       Filename:  alnseed.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/27/2012 11:20:04 AM
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <pthread.h>
#include "kstring.h"
#include "indexio.h"
#include "query.h"
#include "utils.h"
//#include "bwa.h"

//#define N_SEQS 0x40000
#define N_SEQS 100000
typedef struct{
    int se;//1:SE;0:PE
    unsigned int max_tlen, min_tlen;
    int n_mismatch;
    int n_diff;

    unsigned long seed;

    int n_threads;
    char *fn_index;
    char *fn_read1;
    char *fn_read2; 
    char *rg_id;
    int l_read;
    int l_seed;
    int l_overlap;
    int ref;
    int use_sw_extend;
    int print_xa_cigar;
    int print_nm_md;
    kstring_t *cmd;
    uint32_t max_walk; 
    uint32_t max_seed;
    uint32_t max_locate;
    uint32_t max_hits;
    int extend_algo;
    int mismatch_penalty;
    int gapop_penalty;
    int gapext_penalty;

} opt_t;
typedef struct{
    int n_threads;
    int max_diff;
    int max_hits;
    //locationg arguments
    uint32_t max_locate; //max locations per bwt range
    //seeding arguments 
    int seed_only_ref;
    uint32_t max_seed; //max locations per bwt range
    int l_seed;
    int l_overlap;
    //for PE
    uint32_t max_tlen, min_tlen;
    // smith-watemam arguments 
    int gap_op;
    int gap_ex;
    uint8_t flag;
    int filters;
    int filterd;
    int thres_score;
    //sam format arguments;
    int print_xa_cigar;
    int print_nm_md;
    char *rg_id;

} aln_opt_t;

typedef struct{
    uint32_t sp;
    uint32_t ep;
    uint32_t offset;
} sai_t;


typedef kvec_t(sai_t) vec_sai_t;
typedef struct {
    int n_sai_range;
    //auxiliary data for SA index range
    vec_sai_t sai_C;
    vec_sai_t sai_forwardR;
    vec_sai_t sai_backwardR;
    //candidate locations
    vec_uint32_t loci;
    //all hits
    hits_t hits; 
} aux_t;
#ifdef HAVE_THREAD
typedef struct {
    int tid;
    index_t *index;
    int n_query;
    query_t *query;
    aln_opt_t *aln_opt;
    aux_t *aux[2]; 

} thread_aux_t;
#endif
static inline aln_opt_t* aln_opt_init(const opt_t *opt){
    aln_opt_t *aln_opt = calloc(1, sizeof(aln_opt_t));
    aln_opt->n_threads = opt->n_threads;
    aln_opt->l_seed = opt->l_seed;
    aln_opt->max_diff = opt->n_diff; 
    
    aln_opt->min_tlen = opt->min_tlen;
    aln_opt->max_tlen = opt->max_tlen; 
    aln_opt->seed_only_ref = opt->ref;
    aln_opt->l_overlap = opt->l_overlap;
    aln_opt->max_locate = opt->max_locate; 
    aln_opt->max_seed = opt->max_seed; 
    aln_opt->max_hits = 5;
    
    
    aln_opt->gap_op = 3;
    aln_opt->gap_ex = 1;
    aln_opt->flag = 2;
//    opt->filters = 150;
    aln_opt->filters = 0;
    aln_opt->filterd = 20;
//    opt->thres_score = 150;
    aln_opt->thres_score = 50;
    aln_opt->print_nm_md = opt->print_nm_md;
    aln_opt->print_xa_cigar = opt->print_xa_cigar;
    aln_opt->rg_id = opt->rg_id;

    return aln_opt;

}
static inline void aln_opt_destroy(aln_opt_t *aln_opt){
    free(aln_opt); 
}
static inline void aux_resize(aux_t *aux, int n)
{
    aux->n_sai_range = n;
    kv_resize(sai_t, aux->sai_C, n);
    kv_resize(sai_t, aux->sai_forwardR, n);
    kv_resize(sai_t, aux->sai_backwardR, n);
}
static inline aux_t* aux_init(int max_seqLen, int l_seed)
{
    
    aux_t *aux = calloc(1, sizeof(aux_t));

    int n = max_seqLen - l_seed+1;
    aux->n_sai_range = n;
    aux_resize(aux, n);
    kv_resize(uint32_t, aux->loci, 1024);
    kv_resize(hit_t, aux->hits, 1024);
    return aux;
}

static inline void aux_destroy(aux_t *aux)
{
    kv_destroy(aux->sai_C);
    kv_destroy(aux->sai_forwardR);
    kv_destroy(aux->sai_backwardR);
    kv_destroy(aux->loci);
    kv_destroy(aux->hits); 
    free(aux);

}
static inline void aux_reset(aux_t *aux)
{
    //clear sai
    aux->sai_C.n = 0;
    aux->sai_backwardR.n = 0;
    aux->sai_forwardR.n = 0;
    //clear locations 
    aux->loci.n = 0;
    //clear hits
    aux->hits.n = 0;
}
int usage();

void alnse_overlap(index_t *index, query_t *query, aln_opt_t *aln_opt, aux_t *aux_data[2]);
void alnse_nonoverlap(index_t *index, query_t *query, aln_opt_t *aln_opt, aux_t *aux_data[2]); //extend with lv


int alnse_core(const opt_t *opt);

int pairing2(index_t *index, query_t *q0, query_t *q1, const aln_opt_t *aln_opt);
int pairing_singleton(index_t *index, query_t *q0, query_t *q1, aln_opt_t *aln_opt);
int snpaln_sw(index_t *index, uint32_t *start_pos, uint32_t *end_pos, query_t *q, int strand, const aln_opt_t *opt);

int alnpe_core(const opt_t *opt);

int aln_main(int argc, char *argv[]);


#endif
