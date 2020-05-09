/*
 * =====================================================================================
 *
 *       Filename:  alnse.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/01/2013 02:08:59 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */

//cbwt should use forward-search
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#include "ksort.h"
#include "khash.h"
#include "kvec.h"
#include "bwt.h"
#include "aln.h"
#include "editdistance.h"
#include "ssw.h"
#include "sam.h" 

#define __sai_lt(a, b) ((a).ep -(a).sp < (b).ep - (b).sp)
//KHASH_MAP_INIT_INT(32, unsigned char)
KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT(sai, sai_t, __sai_lt)	



#define MAX_LOC_POS 0x40000
#define __MAX(a, b) ((a)>(b)?(a):(b))
#define __MIN(a,b) ((a)<(b)?(a):(b) )





#define LITTLE_SEED_MAX 25
#define PRINT_TIME() \
do{\
    fprintf(stderr, "[aln_main]aln reads ...");\
    fprintf(stderr, "%.2f sec\n", (float)(clock()-tot_aln)/CLOCKS_PER_SEC);\
} while(0)



void alnse_core1(int tid, index_t *index, int n_query, query_t *multi_query, aln_opt_t *aln_opt, aux_t* aux[2]);


static inline int __lt(uint32_t *a, uint32_t *b)
{
    return *a <*b ; 

}
void alnse_seed(index_t *index, uint32_t l_seq, const uint8_t *seq, int l_seed, int split_num, aux_t *aux)
{
    //int m;
    int i, n;       
    int l_longseed = 2*l_seed; 

    cbwt_t *cbwt = index->cbwt;
    rbwt_t *rbwt0 = index->rbwt2->rbwt0; //forward
    rbwt_t *rbwt1 = index->rbwt2->rbwt1; //backward
    lookupTable_t *lkt = index->lkt;
    int l_lkt = lkt->maxLookupLen;
    const uint32_t n_SHARP0 = rbwt0->cumulativeFreq[NT_SHARP];   
    const uint32_t n_SHARP1 = rbwt1->cumulativeFreq[NT_SHARP];   
    
    int n_C, n_for_R, n_back_R;
    sai_t *sa_C, *sa_for_R, *sa_back_R; 
    
    n_C = 0;
    n_for_R = 0;
    n_back_R = 0; 
    sa_C = aux->sai_C.a;
    sa_for_R = aux->sai_forwardR.a;
    sa_back_R = aux->sai_backwardR.a; 
     
   
    int l_win = l_seq/split_num;
    if(l_win < l_seed){
        l_win = l_seed;
        split_num = l_seq/l_win;
    }    
    for(n = 0; n < split_num; ++n )
    {   
        uint32_t k, l, k1, l1;
        uint8_t c;
        //uint32_t primary = cbwt->primary;
        int win_start = l_win * n; 
        int win_end = n < split_num-1 ? win_start +l_win -1 : l_seq-1;

k = 1;//k=0???? 
        l = cbwt->seq_len;
        LKT_lookup_sa(lkt, seq, win_end -l_lkt+1, win_end, &k, &l); //
        
        //Test lookup result
        //uint32_t test_k =1, test_l = cbwt->seq_len;
        //bwt_match_exact_alt(cbwt, l_lkt, seq+win_end-l_lkt+1, &test_k, &test_l);
        //fprintf(stderr, "[Len = l_lkt][LKT -- exact]: [%u, %u] -- [%u, %u]\n",k, l,test_k, test_l ); 
        
        //uint32_t test1_k = 1, test1_l = cbwt->seq_len;
        //uint32_t test2_k = k, test2_l = l;
        //if(bwt_match_exact_alt(cbwt, l_win, seq+win_start, &test1_k, &test1_l)>0 &&      bwt_match_exact_alt(cbwt, win_end-win_start+1-l_lkt, seq+win_start, &test2_k, &test2_l)>0 )
            //fprintf(stderr, "[Len = 25][exact -- lkt]: [%u, %u] -- [%u, %u]\n",test1_k, test1_l,test2_k, test2_l ); 
	    if(k <= l && bwt_match_exact_alt(cbwt, win_end - win_start+1-l_lkt, seq+win_start, &k, &l) >0 ){//lf mapping seq[win_start:win_end-win_start+1-l_lkt]
            if(n_C >= aux->sai_C.m) fprintf(stderr, "n_C:%lu\n", aux->sai_C.m);

            sa_C[n_C].sp = k;
            sa_C[n_C].ep = l;
            sa_C[n_C].offset = win_start;               
            ++n_C;
        }
        //R aln + C aln           
        k = 0;
        l = rbwt0->textLength;   
        if( Rbwt_exact_match_forward(rbwt0, seq+win_start, l_seed, &k, &l) >0 ){//lf mapping seq[win_start::25]
            
            //fprintf(stderr, "[SEQ]: %u\n",n); 
            //fprintf(stderr, "[R range]: [%u, %u]\n",k, l); 
            i = win_start +l_seed;
            while(k <=l && i <= win_end&& i -win_start < l_longseed){//lf mapping extend to seq[win_start+25::25]
                c = seq[i];
                k = rbwt0->cumulativeFreq[c] + Rbwt_BWTOccValue2(rbwt0, k, c, &k1) +1;
                l = rbwt0->cumulativeFreq[c] + Rbwt_BWTOccValue2(rbwt0, l +1, c, &l1);
                if(k1 < l1){//k1 +1<= l1
                    if(n_for_R >= aux->sai_forwardR.m) fprintf(stderr, "n_for_R:%lu\n", aux->sai_forwardR.m);
                    sa_for_R[n_for_R].sp = k1 + n_SHARP0 +1;
                    sa_for_R[n_for_R].ep = l1 + n_SHARP0;
                    sa_for_R[n_for_R].offset = i;
                    ++n_for_R;
                }
                ++i;
            }
            // window 中不包含'#' 
            if(k <= l&&l_win < l_longseed){
                if(n_for_R >= aux->sai_backwardR.m) fprintf(stderr, "n_for_R%lu\n", aux->sai_backwardR.m);

                sa_for_R[n_for_R].sp = k;
                sa_for_R[n_for_R].ep = l;
                sa_for_R[n_for_R].offset = i+1;
                ++n_for_R;
            }
        }
        //R aln + R aln && C aln + R aln 
        k = 0;
        l = rbwt1->textLength;
        if( Rbwt_exact_match_backward(rbwt1, seq+win_end-l_seed+1, l_seed, &k, &l) > 0 ){//lf mapping seq[::l_seed]
            
            //fprintf(stderr, "[SEQ]: %u\n",n); 
            //fprintf(stderr, "[R range]: [%u, %u]\n",k, l); 
            i = win_end - l_seed;
            while(k <=l && i >= win_start && win_end-i <l_longseed){//lf mapping extend to seq[:win_end-l_seed:l_seed]
                c = seq[i];
                k = rbwt1->cumulativeFreq[c] + Rbwt_BWTOccValue2(rbwt1, k, c, &k1) +1; 
                l = rbwt1->cumulativeFreq[c] + Rbwt_BWTOccValue2(rbwt1, l +1, c, &l1);
                if(k1 < l1){
                    if(n_back_R >= aux->sai_backwardR.m) fprintf(stderr, "n_back_R:%lu\n", aux->sai_backwardR.m);

                    sa_back_R[n_back_R].sp = k1 + n_SHARP1 +1;
                    sa_back_R[n_back_R].ep = l1 + n_SHARP1; 
                    sa_back_R[n_back_R].offset = i;
                    ++n_back_R;
                }
                --i;
            }
            if(k <= l && l_win < l_longseed){
                if(n_back_R >= aux->sai_backwardR.m) fprintf(stderr, "n_back_R%lu\n", aux->sai_backwardR.m);

                sa_back_R[n_back_R].sp = k;
                sa_back_R[n_back_R].ep = l;
                sa_back_R[n_back_R].offset = i+1;
                ++n_back_R;
            }
        }
        
    }

    aux->sai_C.n = n_C;
    aux->sai_forwardR.n = n_for_R;
    aux->sai_backwardR.n = n_back_R;

}
#define SKIP_VARIANT 0
//#define OVERLAP_INTV 5
//void alnse_seed_overlap(index_t *index, const uint8_t *seq,uint32_t l_seq, int l_seed, candidate_t *p_can)
void alnse_seed_overlap(index_t *index, uint32_t l_seq, const uint8_t *seq, aln_opt_t *opt, aux_t *aux_data)
{
    //int m;
    int l_seed = opt->l_seed;
    int n;       

    int l_shortseed = l_seed; 

    int OVERLAP_INTV = opt->l_overlap;
    cbwt_t *cbwt = index->cbwt;
    //rbwt_t *rbwt0 = index->rbwt2->rbwt0; //forward
    rbwt_t *rbwt1 = index->rbwt2->rbwt1; //backward

    lookupTable_t *lkt = index->lkt;
    int l_lkt = lkt->maxLookupLen;
   
     
    int n_C, n_back_R;
    sai_t *sa_C, *sa_back_R; 
    
    n_C = 0;
    n_back_R = 0;
    //n_back_R = 0; 
    sa_C = aux_data->sai_C.a;
    //sa_for_R = aux_data->sai_forwardR.a;
    sa_back_R = aux_data->sai_backwardR.a; 
     
    int seed_start, seed_end; 
    for(seed_start= 0; (unsigned int)seed_start < l_seq-l_shortseed+1; ++seed_start )
    {   
        if (seed_start % OVERLAP_INTV != 0) continue;
        seed_end = seed_start+l_shortseed-1;
        uint32_t k, l; uint8_t c;
        
        /* non-snp seeding*/
        k = 1; //k=0???
        l = cbwt->seq_len;
        LKT_lookup_sa(lkt, seq, seed_end -l_lkt+1, seed_end, &k, &l);
        /*
        uint32_t __k=1, __l=cbwt->seq_len;
        if(bwt_match_exact_alt(cbwt, l_lkt, seq+win_end-l_lkt+1, &__k, &__l)>0){
            if(k != __k) fprintf(stderr, "LKT Error !\n");
            if(l != __l) fprintf(stderr, "LKT Error!\n");
        }*/
        if(k <= l && bwt_match_exact_alt(cbwt, l_shortseed-l_lkt, seq+seed_start, &k, &l) >0 ){//lf mapping seq[seed_start:win_end-seed_start+1-l_lkt]
        //if( 0 ){//lf mapping seq[seed_start:win_end-seed_start+1-l_lkt]
            //if(n_C >= p_can->C.m) fprintf(stderr, "n_C:%lu\n", p_can->C.m);
            int l_extend = 0;
            //extend seed util bwt_range < max locate num
            while(l- k > opt->max_seed && l_extend < seed_start){
                bwtint_t ok, ol;
                c = seq[seed_start-l_extend-1];
                if(c > 3) break;
                bwt_2occ(cbwt, k-1, l, c, &ok, &ol);
                if(ok+1 > ol ) { break; } 
                k = cbwt->L2[c]+ok+1, 
                l = cbwt->L2[c]+ol;
                ++l_extend;
                if(l - k <= opt->max_seed){ break;}
            }
#ifdef DEBUG 
        fprintf(stderr, "[Seed C]  [%u - %u], offset= %u, len = %u\n", k, l, seed_start -l_extend, seed_end+1-seed_start+l_extend);

#endif


            sa_C[n_C].sp = k;
            sa_C[n_C].ep = l;
            sa_C[n_C].offset = seed_start-l_extend;               
            ++n_C;
        }
        /* snp-aware seeding */           
        //if(SKIP_VARIANT) continue;
        if(opt->seed_only_ref) continue;
        k = 0;
        l = rbwt1->textLength;   
        if( Rbwt_exact_match_backward(rbwt1, seq+seed_start, l_shortseed, &k, &l) >0 ){//lf mapping seq[win_start::25]
            //fprintf(stderr, "[SEQ]: %u\n",n); 
            //fprintf(stderr, "[R range]: [%u, %u]\n",k, l); 
            //if(n_for_R >= p_can->for_R.m) fprintf(stderr, "n_for_R:%lu\n", p_can->for_R.m);
            int l_extend = 0;
            while( l- k >opt->max_seed && l_extend < seed_start) {
                bwtint_t ok, ol;
                c = seq[seed_start-l_extend-1];
                ok = Rbwt_BWTOccValue(rbwt1, k, c);
                ol = Rbwt_BWTOccValue(rbwt1, l+1, c);
                if(ok+1 > ol ) { break; } 
                k = l = rbwt1->cumulativeFreq[c];
                k += ok+1;
                l += ol;
                ++l_extend;
                if(l-k <= opt->max_seed){break; }
            }
 #ifdef DEBUG 
        fprintf(stderr, "[Seed R]  [%u - %u], offset= %u, len = %u\n", k, l, seed_start-l_extend, seed_end-seed_start+l_extend);

#endif
           
            sa_back_R[n_back_R].sp = k;
            sa_back_R[n_back_R].ep = l;
            sa_back_R[n_back_R].offset = seed_start-l_extend;
            ++n_back_R;
        }
    }

    aux_data->sai_C.n = n_C;
    aux_data->sai_backwardR.n = n_back_R;
    
    ks_introsort_sai(n_C, aux_data->sai_C.a);
    ks_introsort_sai(n_back_R, aux_data->sai_backwardR.a);
    
    //p_can->back_R.n = n_back_R;

}



#define INDEL_MATCH 4 
#define NO_MATCH -1
#define NO_INDEL_MATCH 5
#define FLAG_SUBSTITION 1
#define FLAG_INDEL 2
#define INDEL_FLAG 256
#define MISMATCH_FLAG 1
#define SET_GAP(a, n) (a = INDEL_FLAG&n) 
#define SET_MISMATCH(a, n) (a = MISMATCH_FLAG&n) 
#define IS_INDEL(erro) ( (erro)&INDEL_FLAG)
#define N_DIFF(erro) (erro&255)
#define MAX_LOC_NUM 200

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  alnse_ex_new
 *  Description:  check candidate pos; 
 *                only support mismatch
 * =====================================================================================
 */
#define FORWARD_STRAND 0
#define BACKWARD_STRAND 1
/*
inline void __set_hit(hit_t *hit, uint32_t pos, uint8_t n_err, uint8_t is_gap, uint16_t strand)
{
    
    hit->is_gap = is_gap;
    hit->n_err = n_err;
    hit->pos = pos;
    hit->strand = strand;

}*/
#define code_kmismatch(key) do{\
    if ( (n_mismatch = ed_mismatch(mref->seq, (key), seq, l_seq, max_diff)) >=0 ){\
        if(n_mismatch < max_diff || flag_match == NO_MATCH) {\
            max_diff = n_mismatch;\
            query->is_gap = 0;\
            query->n_diff = n_mismatch;\
            query->strand = strand;\
            query->pos = (key);\
        }\
        flag_match = 1;\
        if(hits->n == hits->m){\
            hits->m = hits->m?hits->m<<1:2; \
            hits->a = realloc(hits->a, sizeof(hit_t) * hits->m);\
        }\
        hit_t *hit = hits->a+hits->n;\
        hit->is_gap = 0;\
        hit->n_diff = n_mismatch;\
        hit->pos = (key);\
        hit->strand = strand;\
        ++hits->n;\
    }\
}while(0) 

#define code_kdiff(key) do{\
    int n_diff;\
    if ( (n_diff = ed_diff(mref->seq, mref->l, (key),l_seq+4,  seq, l_seq, max_diff)) >=0 ){\
        if(n_diff < max_diff ||flag_match == NO_MATCH) {\
            max_diff = n_diff;\
            query->is_gap = 1;\
            query->n_diff = n_diff;\
            query->strand = strand;\
            query->pos = (key);\
        }\
        flag_match = 1;\
        if(hits->n == hits->m){\
            hits->m = hits->m?hits->m<<1:2; \
            hits->a = realloc(hits->a, sizeof(hit_t) * hits->m);\
        }\
        hit_t *hit = hits->a+hits->n;\
        hit->is_gap = 1;\
        hit->n_diff = n_diff;\
        hit->pos = (key);\
        hit->strand = strand;\
        ++hits->n;\
    }\
}while(0)


int alnse_ex_new(index_t *index, query_t *query, int max_diff, uint32_t max_locate, int strand, aux_t *aux[2])
{
    int i;
    uint32_t j, pos;
    cbwt_t *cbwt;
    rbwt_t *rbwt0, *rbwt1;
    mixRef_t *mref;
    
    cbwt = index->cbwt;
    rbwt0 = index->rbwt2->rbwt0;
    rbwt1 = index->rbwt2->rbwt1;
    mref = index->mixRef;
   
    aux_t *aux_data = aux[strand];
    hits_t *hits = &query->hits[strand]; 
    
    int flag_match = NO_MATCH;
    /******************************C Part****************************************************/
    uint8_t *seq = strand == FORWARD_STRAND?query->seq:query->rseq;
    uint32_t l_seq = query->l_seq;
    for(i = 0; i < aux_data->sai_C.n; ++i){
        const sai_t *C = aux_data->sai_C.a+i;
        for(j = C->sp; j <= C->ep && j-C->sp < max_locate; ++j){
             pos = bwt_sa(cbwt, j);
             pos -= C->offset;

             if(pos > mref->l || pos + l_seq > mref->l) {continue;}             
             int n_mismatch = -1;
             code_kmismatch(pos); 
        } 
    }
     
    /******************************R Part forward****************************************************/
    srand(time(0));
    for(i = 0; i < aux_data->sai_forwardR.n; ++i){
        const sai_t *for_R = aux_data->sai_forwardR.a+i;
        uint32_t random_range, random_iter = (uint32_t) -1;
        
        if(for_R->ep - for_R->sp > max_locate){
            random_range = (for_R->ep - for_R->sp)/max_locate;
            for(j = for_R->sp; j <= for_R->ep; ++j){
            	uint32_t iter;
                iter = j - for_R->sp;
                
                if(iter%random_range == 0) {random_iter = rand()%random_range;}

    	    	if(iter%random_range !=random_iter) {continue;} 
            	pos = Rbwt_for_bwt_sa(rbwt0, j);
            	pos -= for_R->offset;

            	if(pos > mref->l || pos + l_seq > mref->l) {continue;}             
             	int n_mismatch = -1;
                code_kmismatch(pos); 
            }
	   }else{
            for(j = for_R->sp; j <= for_R->ep; ++j){
           	    pos = Rbwt_for_bwt_sa(rbwt0, j);
            	pos -= for_R->offset;

            	if(pos > mref->l || pos + l_seq > mref->l) {continue;}             
                int n_mismatch = -1;
                code_kmismatch(pos); 
            }//end for
       }//end else
    }//end for

    /******************************R Part backward****************************************************/
    for(i = 0; i < aux_data->sai_backwardR.n; ++i){
        const sai_t *back_R;
        back_R = aux_data->sai_backwardR.a+i;
        uint32_t random_range, random_iter = (uint32_t)-1;
    	if(back_R->ep - back_R->sp > max_locate){
    	    random_range = (back_R->ep - back_R->sp)/max_locate;
    	
            for(j = back_R->sp; j <= back_R->ep; ++j){
     	    	uint32_t iter;
            	iter = j - back_R->sp;
            	
                if(iter%random_range == 0){random_iter = rand()%random_range;}
            	if(iter%random_range !=random_iter){continue;}
            	pos = Rbwt_back_bwt_sa(rbwt1, j);
            	pos -= back_R->offset;

            	if(pos > mref->l || pos + l_seq > mref->l){continue;}             
                int n_mismatch = -1;
                code_kmismatch(pos);

            	 
            }
    	} else{
       	    for(j = back_R->sp; j <= back_R->ep; ++j){
         	        pos = Rbwt_back_bwt_sa(rbwt1, j);
                	pos -= back_R->offset;

                	if(pos > mref->l || pos + l_seq > mref->l){continue;}             
                    int n_mismatch = -1;
                    code_kmismatch(pos);

            }
        }
    }
    

    return flag_match==NO_MATCH ?NO_MATCH:max_diff;
}
void alnse_locate(index_t *index, uint32_t l_seq, uint32_t max_locate, aux_t *aux_data)
{
    int i;
    uint32_t j, pos = (uint32_t)-1;
    cbwt_t *cbwt; rbwt_t *rbwt0, *rbwt1;
    mixRef_t *mref;
    
    cbwt = index->cbwt;
    rbwt0 = index->rbwt2->rbwt0;
    rbwt1 = index->rbwt2->rbwt1;
    mref = index->mixRef;
    
    //int ret_erro;
    size_t n;
    const sai_t *p;

 
    //int is_gap;
    n = aux_data->sai_C.n;
    p = aux_data->sai_C.a; 
    for(i = 0; i < n; ++i){
        const sai_t *C = p+i;   
        for(j = C->sp; j <= C->ep && j-C->sp <= max_locate; ++j){
            pos = bwt_sa(cbwt, j);
            //fprintf(stderr, "[Seed C] %u - %u\t", pos, C->offset);
            pos -= C->offset;


            if(pos + l_seq > mref->l) {continue;}             
            //k = kh_put(32, h, pos, &ret);
            //kh_value(h, k) += 1;
            kv_push(uint32_t, aux_data->loci, pos);
            if(aux_data->loci.n == MAX_LOC_POS) goto end;
        } 
    }
    
    //fprintf(stderr, "\n");
    srand(time(0));
    n = aux_data->sai_forwardR.n;
    p = aux_data->sai_forwardR.a;
    for(i = 0; i < n; ++i){
        const sai_t *for_R = p+i;
        uint32_t random_range, random_iter = (uint32_t)-1;
        
        if(for_R->ep - for_R->sp > max_locate){
            random_range = (for_R->ep - for_R->sp)/max_locate;
            for(j = for_R->sp; j <= for_R->ep; ++j){
            	uint32_t iter;
                iter = j - for_R->sp;
                
                if(iter%random_range == 0) {random_iter = rand()%random_range;}
    	    	if(iter%random_range !=random_iter) {continue;} 
            	pos = Rbwt_for_bwt_sa(rbwt0, j);
                //fprintf(stderr, "[Seed R for] %u - %u - %u\t", j, pos, for_R->offset);
                pos -= for_R->offset;

            	if(pos > mref->l || pos + l_seq > mref->l) {continue;}             
             	//k = kh_put(32, h, pos, &ret);
                //kh_value(h, k) += 1;

         
                kv_push(uint32_t, aux_data->loci, pos);
                if(aux_data->loci.n == MAX_LOC_POS) goto end;
            }
	}else{
            for(j = for_R->sp; j <= for_R->ep; ++j){
           	    pos = Rbwt_for_bwt_sa(rbwt0, j);
                //fprintf(stderr, "[Seed R for] %u - %u - %u\t", j, pos, for_R->offset);
                pos -= for_R->offset;

            	if(pos + l_seq > mref->l) {continue;}             
                //k = kh_put(32, h, pos, &ret);
                //kh_value(h, k) += 1;

  
                kv_push(uint32_t, aux_data->loci, pos);
                if(aux_data->loci.n == MAX_LOC_POS) goto end;
            }//end for
        }//end else
    }//end for

    n = aux_data->sai_backwardR.n;
    p = aux_data->sai_backwardR.a;
    for(i = 0; i < n; ++i){
        const sai_t *back_R = p+i;
        uint32_t random_range, random_iter = (uint32_t)-1;
    	if(back_R->ep - back_R->sp > max_locate){
    	    random_range = (back_R->ep - back_R->sp)/max_locate;
    	
            for(j = back_R->sp; j <= back_R->ep; ++j){
     	    	uint32_t iter;
            	iter = j - back_R->sp;
            	
                if(iter%random_range == 0){random_iter = rand()%random_range;}
            	if(iter%random_range !=random_iter){continue;}
            	pos = Rbwt_back_bwt_sa(rbwt1, j);
                //fprintf(stderr, "[Seed R back] %u, %u - %u\t", j, pos, back_R->offset);
                pos -= back_R->offset;

            	if(pos > mref->l || pos + l_seq > mref->l){continue;}             
                //k = kh_put(32, h, pos, &ret);
                //kh_value(h, k) += 1;


                kv_push(uint32_t, aux_data->loci, pos);
                if(aux_data->loci.n == MAX_LOC_POS) goto end;
            }
    	} else{
            for(j = back_R->sp; j <= back_R->ep; ++j){
         	        pos = Rbwt_back_bwt_sa(rbwt1, j);
                    pos -= back_R->offset;
                	if(pos > mref->l || pos + l_seq > mref->l){continue;}             

                    
                    kv_push(uint32_t, aux_data->loci, pos);
                    if(aux_data->loci.n == MAX_LOC_POS) goto end;
            }//end for
        }//end if-else
    }//end for

  


end:   
    {
        ks_introsort(uint32_t, aux_data->loci.n, aux_data->loci.a);
        return;       
    }
}



void alnse_locate_alt(index_t *index, uint32_t l_seq, uint32_t max_locate, aux_t *aux_data)
{
    int i;
    uint32_t j, pos = (uint32_t)-1;
    cbwt_t *cbwt;
    rbwt_t *rbwt0, *rbwt1;
    mixRef_t *mref;
    
    cbwt = index->cbwt;
    rbwt0 = index->rbwt2->rbwt0;
    rbwt1 = index->rbwt2->rbwt1;
    mref = index->mixRef;
    
  


    
    //int ret_erro;
    size_t n;
    const sai_t *p;

    /* locate non-snp seed */ 
    //int is_gap;
    n = aux_data->sai_C.n;
    p = aux_data->sai_C.a; 
    for(i = 0; i < n; ++i){
        const sai_t *C = p+i;   
#ifdef DEBUG 
        fprintf(stderr, "[Seed C] [%u - %u], offset= %u\n", C->sp, C->ep, C->offset);

#endif
        for(j = C->sp; j <= C->ep; ++j){
            
            if(j > cbwt->seq_len) {
                fprintf(stderr, "[Seed C] %u - %u\t", pos, C->offset);
                exit(1);
            }
            pos = bwt_sa(cbwt, j);

            pos -= C->offset;
            if(pos + l_seq > mref->l) {continue;}             
            //k = kh_put(32, h, pos, &ret);
            //kh_value(h, k) += 1;
            kv_push(uint32_t, aux_data->loci, pos);
            //if(aux_data->loci.n == MAX_LOC_POS) goto end;
            if(aux_data->loci.n == max_locate) goto end;
        } 
    }
    

    /* locate  forward snp-aware seed*/
    n = aux_data->sai_forwardR.n;
    p = aux_data->sai_forwardR.a;
    for(i = 0; i < n; ++i){
        const sai_t *for_R = p+i;
   
        

        for(j = for_R->sp; j <= for_R->ep; ++j){
            pos = Rbwt_for_bwt_sa(rbwt0, j);
            //fprintf(stderr, "[Seed R] %u - %u\t", pos, for_R->offset);
            pos -= for_R->offset;
            if(pos + l_seq > mref->l) {continue;}             
            kv_push(uint32_t, aux_data->loci, pos);
            if(aux_data->loci.n ==  max_locate) goto end;
        }//end for
    
    }//end for

    /* locate  backward snp-aware seed*/
    n = aux_data->sai_backwardR.n;
    p = aux_data->sai_backwardR.a;
    for(i = 0; i < n; ++i){
        const sai_t *back_R = p+i;
       	int n_skip = (back_R->ep+1- back_R->sp)/MAX_LOC_POS;
        if(n_skip <=0 ) n_skip = 1; 

#ifdef DEBUG 
        fprintf(stderr, "[Seed R] [%u - %u], offset= %u\n", back_R->sp, back_R->ep, back_R->offset);
#endif
        for(j = back_R->sp; j <= back_R->ep; j+= n_skip){
            pos = Rbwt_back_bwt_sa(rbwt1, j);
            pos -= back_R->offset;
            //fprintf(stderr, "[Seed R] %u\n", pos);
            if(pos > mref->l || pos + l_seq > mref->l){continue;}             
            kv_push(uint32_t, aux_data->loci, pos);
            if(aux_data->loci.n == max_locate) goto end;
        }//end for
    }//end for

  


end:   
    {
        ks_introsort(uint32_t, aux_data->loci.n, aux_data->loci.a);
        return;       
    }
}

#define __get_rpac(rpac, l) (((rpac)[(l)>>3]>>4*((l)%8))&15)
int alnse_check_nogap(index_t *index, query_t *query, int max_diff, int strand, aux_t *aux_data)
{
    
    uint32_t l_seq = query->l_seq;
    uint8_t *seq = strand == 0?query->seq:query->rseq;
    //multi_hits_t *hits = &(strand==0?query->multi_hits0:query->multi_hits1);
    hits_t *hits = &aux_data->hits; 
    vec_uint32_t *p_loci = &(aux_data->loci); 
    
    int i;
    uint32_t j, pos = (uint32_t)-1;
   
    
    mixRef_t *mref = index->mixRef;
    


    int flag_match = NO_MATCH;
    { 
        uint32_t n = p_loci->n;
        uint32_t *a = p_loci->a;
        //ks_introsort(uint32_t, n, a);
        //qsort(vec.a, vec.n, sizeof(uint32_t), __lt);
        
        uint32_t pos0, pos1;
        pos0 = (uint32_t) -1;
        for(i = 0; i < n; ++i){
            pos1 = a[i];
            if(pos1 == pos0 || pos1 >= mref->l) continue;
        /*
            fprintf(stderr, "check pos %u\t", pos1);
            int __i;
            for(__i = 0; __i < l_seq; ++__i) fprintf(stderr, "%u\t", __get_rpac(index->mixRef->seq, pos1+__i));
            fprintf(stderr, "\n");
            for(__i = 0; __i < l_seq; ++__i) fprintf(stderr, "%u\t", seq[__i]);
        */
            //fprintf(stderr, "\n");
            int n_mismatch = -1;
            code_kmismatch(pos1);
#ifdef DEBUG 
            int rid = -1;
            bns_coor_pac2real(index->bntseq, pos1, query->l_seq, &rid);   
            fprintf(stderr, "rid = %s, pos = %u, ed = %d\n", index->bntseq->anns[rid].name, pos1 -index->bntseq->anns[rid].offset+1, n_mismatch);
#endif
            pos0 = pos1; 
        }
        return flag_match!=NO_MATCH ?max_diff:NO_MATCH;
    }
}
#define __get_mixref(pac, l) (((pac)[(l)>>3]>>4*((l)%8))&15)
int sw_snp(index_t *index, uint32_t start_pos, uint32_t end_pos, query_t *q, hits_t *hits, int strand, int thres_sc, const aln_opt_t *opt)
{
    extern int8_t score_mat2[256];
    int i;
    int min_score = thres_sc; 
    if((start_pos)>=index->bntseq->l_pac){ 
        fprintf(stderr, "error %s\n", q->name);
        exit(1);
    }
    const uint32_t start = start_pos, end = end_pos;
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
  
    if(result->score1 >= q->b0 && result->read_end1-result->read_begin1+1 >=opt->filterd){
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
        
        //*start_pos = result->ref_begin1+start;
        //*end_pos = result->ref_end1 +start;

        q->seq_start = result->read_begin1;
        q->seq_end = result->read_end1;
    } else{
        /*  
        if(result->score1 >= min_score ||q->b0 == -1) {
            min_score = result->score1;
            q->is_gap = 1;
            q->n_diff = -result->score1;
            q->strand = strand;
            q->pos = result->ref_begin1+start;
        }
        */
        if(hits->n == hits->m){
            hits->m = hits->m?hits->m<<1:2; 
            hits->a = realloc(hits->a, sizeof(hit_t) * hits->m);
        }
        hit_t *hit = hits->a+hits->n;
        hit->is_gap = 1;
        hit->n_diff = -result->score1;
        hit->pos = result->ref_begin1+start;
        hit->strand = strand;
        ++hits->n;
    }
    
    //fprintf(stderr, "\nRef Pos:[%d, %d]\n",result->ref_begin1, result->ref_end1);
    //fprintf(stderr, "Seq Pos:[%d, %d]\n", result->read_begin1, result->read_end1);
    
    int sc = result->score1;  
    free(ref);
    free(read);
    align_destroy(result);
    init_destroy(p);
    return sc;
}


int alnse_check_withgap(index_t *index, query_t *query, int max_diff, int strand, aux_t *aux_data)
{
    uint32_t l_seq = query->l_seq;
    uint8_t *seq = strand == 0?query->seq:query->rseq;
    int i;
    uint32_t j, pos = (uint32_t)-1;


    mixRef_t *mref = index->mixRef;
    
    hit_t hit;
    hits_t *hits = &aux_data->hits;
    int flag_match = NO_MATCH;
    { 
        uint32_t n = aux_data->loci.n;
        uint32_t *a = aux_data->loci.a;
        //ks_introsort(uint32_t, n, a);
        //qsort(vec.a, vec.n, sizeof(uint32_t), __lt);

        uint32_t pos0, pos1;
        pos0 = (uint32_t) -1;
        for(i = 0; i < n; ++i){
            pos1 = a[i];
            if(pos1 == pos0 || pos1+l_seq+4 >= mref->l) continue;
            //fprintf(stderr, "check pos %u\n", pos1);
            code_kdiff(pos1);
            pos0 = pos1; 
        }
        return flag_match!=NO_MATCH ?max_diff:NO_MATCH;
    }
}
int alnse_check_sw(index_t *index, query_t *query, int thres_sc, int strand, const aln_opt_t *opt, aux_t *aux_data)
{
    int i;
    uint32_t l_seq = query->l_seq;
    uint8_t *seq = strand == 0?query->seq:query->rseq;
    uint32_t j, pos = (uint32_t)-1;

    mixRef_t *mref = index->mixRef;
    hit_t hit;
    hits_t *hits = &aux_data->hits;
    int flag_match = NO_MATCH;
    { 
        uint32_t n = aux_data->loci.n;
        uint32_t *a = aux_data->loci.a;
        //ks_introsort(uint32_t, n, a);
        //qsort(vec.a, vec.n, sizeof(uint32_t), __lt);

        uint32_t pos0, pos1;
        pos0 = (uint32_t) -1;
        for(i = 0; i < n; ++i){
            pos1 = a[i];
            if(pos1 == pos0 || pos1+l_seq+4 >= mref->l) continue;
            thres_sc = sw_snp(index, pos1, pos1+l_seq+4, query, hits, strand, thres_sc, opt);
            pos0 = pos1; 
        }
        return thres_sc >= opt->thres_score ?thres_sc:NO_MATCH;
    }
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  alnse_unalned
 *  Description:  print unalned track in fastq format 
 * =====================================================================================
 */
void alnse_unalned(query_t *query, FILE *fp)
{
    int i;
    fprintf(fp,">%s\t%s\n", query->name, query->comment);
    for(i = 0; i < query->l_seq; ++i)
    {
        fputc("ACGTN"[query->seq[i]], fp);
    }
    fputc('\n',fp);
    if(strlen((char *)query->qual)){
        fprintf(fp,"+\n%s\n", query->qual);
    }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  alnse_printCandidate
 *  Description:  
 * =====================================================================================
 */
void alnse_printCandidate(aux_t *aux)
{
    fprintf(stderr, "%lu %lu %lu\n", aux->sai_C.n, aux->sai_forwardR.n, aux->sai_backwardR.n);    
    fprintf(stderr, "%lu %lu %lu\n", aux->sai_C.n, aux->sai_forwardR.n, aux->sai_backwardR.n);    
}

#define __SET_ERROR(e, e0, e1) do{\
    if((e0) == NO_MATCH){\
        if((e1) == NO_MATCH) (e) = NO_MATCH;\
        else (e) = (e1);\
    } else {\
        if((e1) == NO_MATCH) (e) = (e0);\
        else{\
            (e) = (e0) < (e1)?(e0):(e1);\
        }\
    }\
} while(0)
    
/* 
 *[MaA ===  FUNCTION  ======================================================================
 *l         Name:  alnse_overlap
 *  Description:  1. single end aln main func(overlap seeding version)
 *    [MaA            2. not work
 * =====================================================================================
 */
//void alnse_overlap(index_t *index, query_t *query, int max_diff, uint32_t max_locate, aux_t *aux)
#define SKIP_NOGAP 0
void alnse_overlap(index_t *index, query_t *query, aln_opt_t* aln_opt, aux_t *aux_data[2])
{
   

    uint32_t max_locate = aln_opt->max_locate;
    int max_diff = aln_opt->max_diff; 
    
    int l_seq = query->l_seq;
    uint8_t *seq= query->seq;
    uint8_t *rseq = query->rseq;
#ifdef DEBUG
    if(seq == NULL|| rseq == NULL){
        fprintf(stderr, "[%s]: seq[%u] == '' !!!", __func__, (unsigned int )query); 
        exit(1); 
    } 
#endif
    if(max_diff < 0){ max_diff = l_seq/10; }
    //generate genomic loci 
    //fprintf(stderr, "\n[Seeding] %s\n", query->name);
    alnse_seed_overlap(index, l_seq, seq, aln_opt, aux_data[0]);
    alnse_locate(index, l_seq, max_locate, aux_data[0]);
    alnse_seed_overlap(index, l_seq, rseq, aln_opt, aux_data[1]);
    alnse_locate(index, l_seq, max_locate, aux_data[1]);
    /*  
    fprintf(stderr, "\n");
    //check without indel 
    fprintf(stderr, "[Extending] %s\n", query->name);
    fprintf(stderr, "Check no gap...\n");
    */
    int n_mismatch0 = NO_MATCH, n_mismatch1= NO_MATCH;
    if(!SKIP_NOGAP){
max_diff = 3;
        n_mismatch0 = alnse_check_nogap(index, query, max_diff, 0, aux_data[0]);
        if(n_mismatch0 !=NO_MATCH &&n_mismatch0<max_diff) max_diff = n_mismatch0;
        n_mismatch1 = alnse_check_nogap(index, query, max_diff, 1, aux_data[1]);
        if(n_mismatch1 !=NO_MATCH &&n_mismatch1<max_diff) max_diff = n_mismatch1;
    }

    //fprintf(stderr, "pos = %u\n", query->pos);
    //fprintf(stderr, "Check with gap...\n");
    //check with indel 
    if(n_mismatch0 == NO_MATCH && n_mismatch1 == NO_MATCH){
if(max_diff < 0){ max_diff = l_seq/10; }
        int n_diff0 = alnse_check_withgap(index, query, max_diff, 0, aux_data[0]);
        if(n_diff0!=NO_MATCH &&n_diff0<max_diff) max_diff = n_diff0;
        int n_diff1 = alnse_check_withgap(index, query, max_diff, 1, aux_data[1]);
        if(n_diff1!=NO_MATCH &&n_diff1<max_diff) max_diff = n_diff1;
        if(n_diff0 == NO_MATCH && n_diff1 == NO_MATCH) max_diff = NO_MATCH; 
    }

    //fprintf(stderr, "pos = %u\n", query->pos);
    query_set_hits(query, aln_opt->max_hits, &aux_data[0]->hits, &aux_data[1]->hits); 
    //fprintf(stderr, "pos = %u\n", query->pos);
    
    return;    
    



}
void alnse_overlap_alt(index_t *index, query_t *query, aln_opt_t* aln_opt, aux_t *aux_data[2])
{
   

    uint32_t max_locate = aln_opt->max_locate;
    int max_diff = aln_opt->max_diff; 
    
    int l_seq = query->l_seq;
    uint8_t *seq= query->seq;
    uint8_t *rseq = query->rseq;
#ifdef DEBUG
    if(seq == NULL|| rseq == NULL){
        fprintf(stderr, "[%s]: seq[%u] == '' !!!", __func__, (unsigned int )query); 
        exit(1); 
    } 
    fprintf(stderr, "\n[aligning] %s\n", query->name);
#endif
    if(max_diff < 0){ max_diff = l_seq/10; }
    //generate genomic loci 

    alnse_seed_overlap(index, l_seq, seq, aln_opt, aux_data[0]);
    alnse_locate_alt(index, l_seq, aln_opt->max_locate, aux_data[0]);
    alnse_seed_overlap(index, l_seq, rseq, aln_opt, aux_data[1]);
    alnse_locate_alt(index, l_seq, aln_opt->max_locate, aux_data[1]);
    //fprintf(stderr, "\n");
    //check without indel 

    //fprintf(stderr, "[Extending] %s\n", query->name);

#ifdef DEBUG
    fprintf(stderr, "Check no gap...");
#endif
    int n_mismatch0 = NO_MATCH, n_mismatch1= NO_MATCH;
    if(!SKIP_NOGAP){
max_diff = 3;
        n_mismatch0 = alnse_check_nogap(index, query, max_diff, 0, aux_data[0]);
        if(n_mismatch0 !=NO_MATCH &&n_mismatch0<max_diff) max_diff = n_mismatch0;
        n_mismatch1 = alnse_check_nogap(index, query, max_diff, 1, aux_data[1]);
        if(n_mismatch1 !=NO_MATCH &&n_mismatch1<max_diff) max_diff = n_mismatch1;
    }
#ifdef DEBUG
    fprintf(stderr, "Check with gap...");
#endif
    //check with indel 
    if(n_mismatch0 == NO_MATCH && n_mismatch1 == NO_MATCH){
max_diff = l_seq/10;
        int n_diff0 = alnse_check_withgap(index, query, max_diff, 0, aux_data[0]);
        if(n_diff0!=NO_MATCH &&n_diff0<max_diff) max_diff = n_diff0;
        int n_diff1 = alnse_check_withgap(index, query, max_diff, 1, aux_data[1]);
        if(n_diff1!=NO_MATCH &&n_diff1<max_diff) max_diff = n_diff1;
        if(n_diff0 == NO_MATCH && n_diff1 == NO_MATCH) max_diff = NO_MATCH; 
    }
    query_set_hits(query, aln_opt->max_hits, &aux_data[0]->hits, &aux_data[1]->hits); 
    
    return;    
    



}
void alnse_overlap_sw(index_t *index, query_t *query, aln_opt_t* aln_opt, aux_t *aux_data[2])
{
   

    uint32_t max_locate = aln_opt->max_locate;
    int max_diff = aln_opt->max_diff; 
    
    int l_seq = query->l_seq;
    uint8_t *seq= query->seq;
    uint8_t *rseq = query->rseq;
#ifdef DEBUG
    if(seq == NULL|| rseq == NULL){
        fprintf(stderr, "[%s]: seq[%u] == '' !!!", __func__, (unsigned int )query); 
        exit(1); 
    } 
    fprintf(stderr, "\n[aligning] %s\n", query->name);
#endif
    if(max_diff < 0){ max_diff = l_seq/10; }
    //generate genomic loci 

    alnse_seed_overlap(index, l_seq, seq, aln_opt, aux_data[0]);
    alnse_locate_alt(index, l_seq, aln_opt->max_locate, aux_data[0]);
    alnse_seed_overlap(index, l_seq, rseq, aln_opt, aux_data[1]);
    alnse_locate_alt(index, l_seq, aln_opt->max_locate, aux_data[1]);
    //fprintf(stderr, "\n");
    //check without indel 

    //fprintf(stderr, "[Extending] %s\n", query->name);

#ifdef DEBUG
    fprintf(stderr, "Check no gap...");
#endif
    int n_mismatch0 = NO_MATCH, n_mismatch1= NO_MATCH;
    if(!SKIP_NOGAP){
max_diff = 3;
        n_mismatch0 = alnse_check_nogap(index, query, max_diff, 0, aux_data[0]);
        if(n_mismatch0 !=NO_MATCH &&n_mismatch0<max_diff) max_diff = n_mismatch0;
        n_mismatch1 = alnse_check_nogap(index, query, max_diff, 1, aux_data[1]);
        if(n_mismatch1 !=NO_MATCH &&n_mismatch1<max_diff) max_diff = n_mismatch1;
    }
#ifdef DEBUG
    fprintf(stderr, "Check with gap...");
#endif
    //check with indel 
    if(n_mismatch0 == NO_MATCH && n_mismatch1 == NO_MATCH){
        int thres_score = aln_opt->thres_score;
        int sc0 = alnse_check_sw(index, query, thres_score, 0, aln_opt, aux_data[0]);
        if(sc0>thres_score) thres_score = sc0;
        int sc1 = alnse_check_sw(index, query, thres_score, 1, aln_opt, aux_data[1]);
        if(sc1>thres_score) thres_score = sc1;
    }
    //query_set_hits(query, aln_opt->max_hits, &aux_data[0]->hits, &aux_data[1]->hits); 
    
               
    return;    
    



}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  alnse_once_test2
 *  Description:  1. single end aln main func(non-overlap seeding version)
 *                2. using now
 * =====================================================================================
 */
//void alnse_once_test2(index_t *index, query_t *query, int max_diff, uint32_t max_locate, aux_t *aux) //extend with lv

void alnse_nonoverlap(index_t *index, query_t *query, aln_opt_t *aln_opt, aux_t *aux_data[2]) //extend with lv
{
    int max_diff = aln_opt->max_diff;
    int max_locate = aln_opt->max_locate;
    
    uint8_t *seq= query->seq;
    uint8_t *rseq = query->rseq;
    int l_seq = query->l_seq;
#ifdef DEBUG
    if(seq == NULL|| rseq == NULL){
        fprintf(stderr, "[%s]: seq[%u] == '' !!!", __func__, (unsigned int )query); 
        exit(1); 
    } 



#endif 

    int n_mismatch, n_mismatch0=NO_MATCH, n_mismatch1=NO_MATCH;   
    uint32_t l_shortseed = aln_opt->l_seed; 
    uint32_t l_longseed = aln_opt->l_seed*2; 
    if(max_diff < 0 || max_diff > (l_seq+l_shortseed-1)/l_shortseed){
        if(max_diff < 0) max_diff = 30;
        int n_seed = (l_seq+l_longseed-1)/l_longseed;
        /* 
        alnse_seed(index, l_seq, seq, l_longseed, n_mismatch, aux_data[0]);
        n_mismatch0 = alnse_ex_new(index, query, n_mismatch, max_locate, 0, aux_data);
        if(n_mismatch0 != NO_MATCH && n_mismatch0 < n_mismatch) n_mismatch = n_mismatch0;
        alnse_seed(index, l_seq, rseq, l_longseed, n_mismatch, aux_data[1]);
        n_mismatch1 = alnse_ex_new(index, query, n_mismatch, max_locate, 1, aux_data);
        if(n_mismatch1 != NO_MATCH && n_mismatch1 < n_mismatch) n_mismatch = n_mismatch1;
        */

        //generate genomic loci 
        alnse_seed(index, l_seq, seq, l_longseed, n_seed, aux_data[0]);
        alnse_locate(index, l_seq, max_locate, aux_data[0]);
        alnse_seed(index, l_seq, rseq, l_longseed, n_seed, aux_data[1]);
        alnse_locate(index, l_seq, max_locate, aux_data[1]);

 
        //check without indel 
        int max_diff0 = n_seed;
        if(!SKIP_NOGAP){
            n_mismatch0 = alnse_check_nogap(index, query, max_diff0, 0, aux_data[0]);
            if(n_mismatch0 !=NO_MATCH &&n_mismatch0<max_diff0) max_diff0 = n_mismatch0;
            n_mismatch1 = alnse_check_nogap(index, query, max_diff0, 1, aux_data[1]);
            if(n_mismatch1 !=NO_MATCH &&n_mismatch1<max_diff0) max_diff0 = n_mismatch1;
        }
        //check with indel 
        if(n_mismatch0 == NO_MATCH && n_mismatch1 == NO_MATCH){
            int n_diff0 = alnse_check_withgap(index, query, max_diff0, 0, aux_data[0]);
            if(n_diff0!=NO_MATCH &&n_diff0<max_diff0) max_diff0 = n_diff0;
            int n_diff1 = alnse_check_withgap(index, query, max_diff0, 1, aux_data[1]);
            if(n_diff1!=NO_MATCH &&n_diff1<max_diff0) max_diff0 = n_diff1;
            if(n_diff0 == NO_MATCH && n_diff1 == NO_MATCH) max_diff0 = NO_MATCH; 
        }
        if(max_diff0 == NO_MATCH){       
            n_seed = (l_seq+l_shortseed-1)/l_shortseed; 
            aux_reset(aux_data[0]);
            aux_reset(aux_data[1]);
            //generate genomic loci 
            alnse_seed(index, l_seq, seq, l_shortseed, n_seed, aux_data[0]);
            alnse_locate(index, l_seq, max_locate, aux_data[0]);
            alnse_seed(index, l_seq, rseq, l_shortseed, n_seed, aux_data[1]);
            alnse_locate(index, l_seq, max_locate, aux_data[1]);

     
            //check without indel 
            if(!SKIP_NOGAP){
                n_mismatch0 = alnse_check_nogap(index, query, max_diff, 0, aux_data[0]);
                if(n_mismatch0 !=NO_MATCH &&n_mismatch0<max_diff) max_diff = n_mismatch0;
                n_mismatch1 = alnse_check_nogap(index, query, max_diff, 1, aux_data[1]);
                if(n_mismatch1 !=NO_MATCH &&n_mismatch1<max_diff) max_diff = n_mismatch1;
            }
            //check with indel 
            if(n_mismatch0 == NO_MATCH && n_mismatch1 == NO_MATCH){
                int n_diff0 = alnse_check_withgap(index, query, max_diff, 0, aux_data[0]);
                if(n_diff0!=NO_MATCH &&n_diff0<max_diff) max_diff = n_diff0;
                int n_diff1 = alnse_check_withgap(index, query, max_diff, 1, aux_data[1]);
                if(n_diff1!=NO_MATCH &&n_diff1<max_diff) max_diff = n_diff1;
                if(n_diff0 == NO_MATCH && n_diff1 == NO_MATCH) max_diff = NO_MATCH; 
            }
            
         
        }        
        query_set_hits(query, aln_opt->max_hits, &aux_data[0]->hits, &aux_data[1]->hits); 
    }
    //kv_destroy(multi_hits);
    return;    
    
}

int g_n_seqs;
pthread_rwlock_t rwlock;

#ifdef HAVE_THREAD

static void *alnse_worker(void *data)
{
    thread_aux_t *d = (thread_aux_t*)data;
    //alnse_core1(d->tid, d->index, d->n_query, d->query, d->aln_opt, d->aux); 
    alnse_core_thread(d->tid, d->index, d->n_query, d->query, d->aln_opt, d->aux); 
    return 0;
}
#endif
#define MAX_N_PERSEQ 200
void alnse_core_thread(int tid, index_t *index, int n_query, query_t *multi_query, aln_opt_t *aln_opt, aux_t* aux[2])
{
    int i;
    while(1){
	    pthread_rwlock_wrlock(&rwlock);
        i = g_n_seqs++; 
        pthread_rwlock_unlock(&rwlock);
        if(i >= n_query) break;

        //fprintf(stderr, "[tid:%d]run seq %d", tid, i); 
        //fprintf(stderr, "aln_opt:%d", aln_opt->n_threads); 

        query_t *query = multi_query+i;
        
        if(query->n_ambiguous > MAX_N_PERSEQ) continue;
        if(query->l_seq -aln_opt->l_seed+1 > aux[0]->n_sai_range){
            int n_sai_range = query->l_seq-aln_opt->l_seed+1; 
            aux_resize(aux[0], n_sai_range); 
            aux_resize(aux[1], n_sai_range); 
        }
        aux_reset(aux[0]); 
        aux_reset(aux[1]); 
        if(aln_opt->l_overlap > 0) alnse_overlap_alt(index, query, aln_opt, aux);
        else alnse_nonoverlap(index, query, aln_opt, aux); 
        query_gen_cigar(index->mixRef->l, index->mixRef->seq, query);
        aln_samse(index, query, aln_opt);
   }
    
}





void alnse_core1(int tid, index_t *index, int n_query, query_t *multi_query, aln_opt_t *aln_opt, aux_t* aux[2])
{
    int i;
    for(i=0; i < n_query; ++i){ 
#ifdef HAVE_THREAD
        if (i % aln_opt->n_threads != tid) continue;
#endif
        //fprintf(stderr, "[tid:%d]run seq %d", tid, i); 
        //fprintf(stderr, "aln_opt:%d", aln_opt->n_threads); 

        query_t *query = multi_query+i;
        
        if(query->n_ambiguous > MAX_N_PERSEQ) continue;
        if(query->l_seq -aln_opt->l_seed+1 > aux[0]->n_sai_range){
            int n_sai_range = query->l_seq-aln_opt->l_seed+1; 
            aux_resize(aux[0], n_sai_range); 
            aux_resize(aux[1], n_sai_range); 
        }
        aux_reset(aux[0]); 
        aux_reset(aux[1]); 
        //if(aln_opt->l_overlap > 0) alnse_overlap_alt(index, query, aln_opt, aux);
        if(aln_opt->l_overlap > 0){ 
            if(aln_opt->extend_algo == EXTEND_SW) alnse_overlap_sw(index, query, aln_opt, aux);
            else{ 
                alnse_overlap_alt(index, query, aln_opt, aux);
            }
            query_gen_cigar(index->mixRef->l, index->mixRef->seq, query);
        } else{
            fprintf(stderr, "Shouldn't be here!!!\n", aln_opt->n_threads); 
            exit(1);
            alnse_nonoverlap(index, query, aln_opt, aux); 
        }
        //query_gen_cigar(index->mixRef->l, index->mixRef->seq, query);
        aln_samse(index, query, aln_opt);
    }
    
}
int alnse_core(const opt_t *opt)
{
    fprintf(stderr, "[alnse_core]:  Start single end alignment!\n");
    
    
   
    //clock_t t = clock();     
    double t_real = realtime();     
    aln_opt_t *aln_opt = aln_opt_init(opt);
    fprintf(stderr, "[%s]:  Reload index...\n", __func__); 
    index_t *index = alnse_index_reload(opt->fn_index);
    //if(opt->use_sw_extend) index->pac = pac_restore(index->bntseq); 
    fprintf(stderr, "%lf sec escaped.\n", realtime()-t_real);
    t_real = realtime();
/*
 
    bwt_cal_sa_intv_on_ref(index->cbwt, index->cbwt->sa_intv);
    exit(1);
*/  
    aux_t *aux[2];
#ifdef HAVE_THREAD
    
    pthread_attr_t attr;
    pthread_t *tid; 
    thread_aux_t *thread_data;
    pthread_rwlock_init(&rwlock, NULL); 
    
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

    queryio_t *qs = query_open(opt->fn_read1);
    query_t *multiSeqs = calloc(N_SEQS, sizeof(query_t));
    if(multiSeqs == NULL){
        fprintf(stderr, "[%s]:  alloc mem fail!\n", __func__);
        exit(1);
    }
    fprintf(stderr, "[%s]:  Begin to align short reads...\n", __func__); 
    aln_samhead(opt, index->bntseq);
    int l_seed = opt->l_seed;
    int n, i, n_tot=0;
    while( (n = query_read_multiSeqs(qs, N_SEQS, multiSeqs)) >0){
        n_tot += n;
        //if(n_tot < 3000000) continue;
#ifdef HAVE_THREAD
        g_n_seqs = 0;
        if(opt->n_threads <=1){
            alnse_core1(0, index, n, multiSeqs, aln_opt, aux); 
        }else{
            int j;
            for (j = 0; j < opt->n_threads; ++j) {
                thread_data[j].query = multiSeqs;
                thread_data[j].n_query = n; 
                pthread_create(&tid[j], &attr, alnse_worker, thread_data + j);
            }
            for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0); 
        }
#else            
        alnse_core1(0, index, n, multiSeqs, aln_opt, aux); 
#endif       
        for(i = 0; i< n; ++i){
            query_t *query = multiSeqs+i;
            //aln_samse(index, query, opt);
            //printf("%s\n", query->sam->s);
            puts(query->sam->s);
            query_destroy(query);
        } 
        
        memset(multiSeqs, 0, N_SEQS*sizeof(query_t));

        //if(n_tot%10000==0 && n_tot!=0) fprintf(stderr, "%d reads have been aligned!\n", n_tot);
        fprintf(stderr, "%d reads have been aligned!\n", n_tot);
    }
    //fprintf(stderr, "[%s]: total %.2f sec escaped\n", __func__, (float)(clock()-t)/CLOCKS_PER_SEC);
    fprintf(stderr, "[%s]: total %lf sec escaped\n", __func__, realtime()-t_real);

    
#ifdef HAVE_THREAD
    pthread_rwlock_destroy(&rwlock);
    if(opt->n_threads<=1){
        aux_destroy(aux[0]);
        aux_destroy(aux[1]);
    }else{
        int j;
        for(j=0; j < opt->n_threads; ++j){
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
    free(multiSeqs);
    query_close(qs);
    alnse_index_destroy(index);
    aln_opt_destroy(aln_opt);



    return EXIT_SUCCESS;



}
