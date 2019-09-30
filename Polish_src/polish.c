/*
 * =====================================================================================
 *
 *       Filename:  polish.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/26/2014 03:11:46 PM
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
#include <getopt.h>
#include <stdint.h>
#include "ssw.h"
#include "lv.h"
#include "samParser.h"
#include "bntseq.h"
#include "khash.h"
#include "ksort.h"

#define SAM_FLAG_PAIRED  0x0001
#define SAM_FLAG_PROPER_PAIRED  0x0002
#define SAM_FLAG_UNMAPED  0x0004
#define SAM_FLAG_MATE_UNMAPPED 0x0008
#define SAM_FLAG_STRAND 0x0010
#define SAM_FLAG_MATE_STRAND 0x0020
#define SAM_FLAG_MATE_READ1 0x0040
#define SAM_FLAG_MATE_READ2 0x0080
#define SAM_FLAG_SECANDARY_ALN 0x0100


/***********OPT************************************************************************/
#define ALGO_SW 0
#define ALGO_LV 1

#define __set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))
#define ABS(a, b)(a<b?b-a:a-b)
static int8_t score_mat[25] = { 2, -2, -2, -2,0,
                                -2,  2, -2, -2,0, 
                                -2, -2,  2, -2,0, 
                                -2, -2, -2,  2,0,
                                0,  0,  0,  0,  0};

const int GAP_OP =3;
const int GAP_EX =1;
typedef struct{
    int distance_algorithm;
    int is_paired;
} opt_t;
static inline opt_t* opt_init(){
    opt_t *opt = calloc(1, sizeof(opt_t));
    opt->distance_algorithm = ALGO_LV;
    return opt;
}
static inline void opt_destroy(opt_t *opt){
    free(opt);
}

char* const short_options = "shp";
struct option long_options[] = {
    { "sw",     0,   NULL,    's'   },
    { "help",     0,   NULL,    'h' },
    { "pe",     0,   NULL,    'p'   },
    { 0,     0,   0,    0   }
};

/***********OPT END********************************************************************/
static inline void print_help(){
    fprintf(stderr, "\n");
    fprintf(stderr, "polish  [OPT]  <index.prefix>  <SAM>\n\n");
    fprintf(stderr, "OPT:    -s, --sw    use smith-waterman to compute distance   [landau-vishkin]\n");
    fprintf(stderr, "        -h, --help  print help\n");
    fprintf(stderr, "        -p, --pe    paired end mode\n");
    fprintf(stderr, "\n");
}
static inline uint8_t *pac_restore(const bntseq_t *bns)
{
    uint8_t *pac;
    pac = calloc(bns->l_pac/4+2,1);
    if(pac == NULL){
        fprintf(stderr, "[pac_restore]: allocate mem fail!\n");
        exit(1); 
    }
    fread(pac, 1, bns->l_pac/4+2, bns->fp_pac);
    
    return pac;
}
int __get_refseq(unsigned char *refseq, int l_seq, unsigned char *pac, unsigned int l_pac, unsigned int start)
{
    if(start> l_pac){ 
        fprintf(stderr, "[Error]: Out of reference length!\n"); 
        exit(1);
    }
    if(start+l_seq >l_pac){ l_seq = l_pac -start;}
    int i;
    for(i =0; i < l_seq; ++i) refseq[i] = __get_pac(pac, start+i); 
    return l_seq;
}

KHASH_MAP_INIT_STR(str, int)
#define __lt_offset(a, b) ((a).offset < (b).offset)
KSORT_INIT(hits, hit_t, __lt_offset)
static inline void __SWAP_HIT(hit_t *h0, hit_t* h1)
{
#define __SWAP_GENRIC(type, a, b) do{\
    type tmp = (a);\
    (a) = (b);\
    (b) = tmp;} while(0) 
    
    __SWAP_GENRIC(char *, h0->chrom, h1->chrom);
    __SWAP_GENRIC(unsigned int, h0->offset, h1->offset);
    __SWAP_GENRIC(unsigned int, h0->pos, h1->pos);
    __SWAP_GENRIC(int, h0->score, h1->score);
}
void rm_repeat_hits(multi_hits_t *hits)
{



    if (hits->n == 0) return;
    size_t n = 1, i;
    unsigned int last_offset = hits->a[0].offset;
    for(i =1; i < hits->n; ++i){
        unsigned int offset = hits->a[i].offset;
        if(offset != last_offset) {
            __SWAP_HIT(hits->a+i, hits->a+n);
            ++n;            
            last_offset = offset;
        }        
    }
    hits->n = n;
}
#define LOW_BOUNDARY -1
#define IN_RANGE 0 
#define HIGH_BOUNDARY 1

const int MIN_ISIZE = 350;
const int MAX_ISIZE = 650;
const int MAX_DISTANCE = 13;
static inline int CHECK_IN_RANGE(uint32_t a, uint32_t b, uint32_t small, uint32_t large){
    uint32_t r = a>b?a-b:b-a;
    if(a>b || r < small) return LOW_BOUNDARY;
    else if(r>large) return HIGH_BOUNDARY;
    else return IN_RANGE;
}
int __pairing(multi_hits_t *forward_hit, multi_hits_t *backward_hit)
{
    //make sure forward and backward hits are sorted by offset
    int n = 0;
    
    unsigned int n0 = forward_hit->n;
    hit_t *h0 = forward_hit->a;
    unsigned int n1 = backward_hit->n;
    hit_t *h1 = backward_hit->a;
    if(n0==0|| n1==0) return n;
    unsigned int i=0, j=0;

    
    while(i < n0 && j < n1){
        unsigned int offset0 = h0[i].offset;
        unsigned int offset1 = h1[j].offset;
        int r = CHECK_IN_RANGE(offset0, offset1, MIN_ISIZE, MAX_ISIZE);
        if(r == LOW_BOUNDARY){
            j += 1;
        
        } else if(r == HIGH_BOUNDARY ){
            i += 1;
        } else{
            __SWAP_HIT(h0+n, h0+(i++));
            __SWAP_HIT(h1+n, h1+(j++));
            ++n; 
        }
    
    
    }


    return n;
}
void gen_cigar(SAM_t *sam, unsigned int l_pac, char *pac, unsigned int offset, int ALGO)
{
    int l_refseq = sam->l_seq; 
    char *refseq = calloc((l_refseq+15)/8*8, 1);
    l_refseq = __get_refseq(refseq, l_refseq, pac, l_pac, offset);


    /*
    if(i==1) {//reverse
        for(k=0; k < l_refseq0; ++k) refseq0[k] = 4-refseq0[k];
        for(k=0; k < l_refseq0/2; ++k) SWAP(refseq0[k], refseq0[l_refseq0-1-k]); 
    }*/
    //recaluate editdistance
    const int8_t *query = (sam->strand==0)?((const int8_t*)sam->nst_seq):((const int8_t*)sam->nst_rseq);
    unsigned int l_query = sam->l_seq;
    int d = sam->multi_hits[sam->strand]->a[sam->primary].score;
    if(ALGO == ALGO_SW){

        int d = sam->multi_hits[sam->strand]->a[sam->primary].score;
        s_profile *p = ssw_init(query, l_query,score_mat, 5, 1 ); 
        s_align *result = ssw_align(p, (const int8_t *)refseq, l_refseq, GAP_OP, GAP_EX, 2, d, 0, l_query/2); 
        if(result->score1 != d) {
            fprintf(stderr, "push cigar error!\n");
            exit(1); 
        }
        int j;
        if(result->read_begin1 != 0){
            ksprintf(&sam->cigar, "%dS",result->read_begin1);//Cigar
        }
        for(j = 0; j < result->cigarLen; ++j) {
            unsigned int cigar = result->cigar[j];    
            ksprintf(&sam->cigar, "%u%c",cigar>>4,"MID"[cigar&15]);//Cigar
        }
        if(result->read_end1+1 != sam->l_seq){
            ksprintf(&sam->cigar, "%dS",sam->l_seq-result->read_end1-1);//Cigar
        }
        
   
        init_destroy(p);
        align_destroy(result);
    } else if(ALGO == ALGO_LV){

        if(d == -MAX_DISTANCE){
            ksprintf(&sam->cigar, "*");
        } else{
            int ret_d = computeEditDistanceWithCigar((const char*)refseq, l_refseq, (char *)query, l_query, -(d), sam->cigar.s, sam->cigar.m, 1, COMPACT_CIGAR_STRING );  
            if(-d != ret_d){
                fprintf(stderr, "push cigar error!\n");
                exit(1); 

            } 
        } 
    } else{
        fprintf(stderr, "distance_algorithm is not available");
        exit(1);
    }
    free(refseq);
    return;

    

}
void polish_sam_se(SAM_t *sam0)
{
    
    int i;
    unsigned char flag0 = 0;
    int strand0 = sam0->strand; 
    int iter0 = sam0->primary;
    int read0_mapped = strand0 != -1?1:0;  
 


    unsigned int pos0 = read0_mapped?sam0->multi_hits[strand0]->a[iter0].pos:(unsigned int)-1;
    char *chrom0 = read0_mapped?sam0->multi_hits[strand0]->a[iter0].chrom:NULL;

     
        
    {
        printf("%s\t", sam0->seq_name);//QNMAE 

        if(strand0 ==1) flag0 |= SAM_FLAG_STRAND;

        flag0 |= SAM_FLAG_MATE_READ1; 
        if(read0_mapped != 1) flag0 |= SAM_FLAG_UNMAPED;

        printf("%u\t", flag0);//FLAG
        
        if(read0_mapped != 1) printf("*\t0\t");
        else{ 
            printf("%s\t%u\t", chrom0, pos0);//RNAME POS
        }
        //MAPQ (to be fixed)
        if(sam0->b1 == -100000 && sam0->b0 != -100000) printf("60\t");
        else printf("0\t");
        if(read0_mapped) printf("%s\t", sam0->cigar.s); 
        else printf("*\t");
        //MRNAME MPOS isize
        printf("*\t0\t0\t");
        unsigned char *s = (strand0==0)?(sam0->nst_seq):(sam0->nst_rseq);   
        for(i =0; i < sam0->l_seq; ++i) putchar("ACGT"[s[i]]);//SEQ
        putchar('\t');
        //QUAL
        s = (unsigned char *)sam0->qual; 
        if(sam0->flag & SAM_FLAG_STRAND){
            if(strand0 == 0){
                for(i =sam0->l_seq-1; i >=0; --i) putchar(s[i]);    
            } else{ printf("%s\t", s);} 
        } else{
            if(strand0 == 0){
                printf("%s\t", s);
            } else{
                for(i =sam0->l_seq-1; i >=0; --i) putchar(s[i]);
            }
        }
        putchar('\n'); 
    }
    
}
void polish_sam_pe(SAM_t *sam0, SAM_t *sam1, int is_properpaired)
{
    
    int i;
    unsigned char flag0 = 0, flag1= 0;

    
    int strand0, strand1, iter0, iter1;
    strand0 = sam0->strand; iter0 = sam0->primary;
    strand1 = sam1->strand; iter1 = sam1->primary;
    
    int read0_mapped = strand0 != -1?1:0;  
    int read1_mapped = strand1 != -1?1:0;  


    unsigned int pos0 = read0_mapped?sam0->multi_hits[strand0]->a[iter0].pos:(unsigned int)-1;
    unsigned int pos1 = read1_mapped?sam1->multi_hits[strand1]->a[iter1].pos:(unsigned int)-1;
    char *chrom0 = read0_mapped?sam0->multi_hits[strand0]->a[iter0].chrom:NULL;
    char *chrom1 = read1_mapped?sam1->multi_hits[strand1]->a[iter1].chrom:NULL;
     
        
    {
        printf("%s\t", sam0->seq_name);//QNMAE 
        flag0 |= SAM_FLAG_PAIRED; 
        if(is_properpaired) flag0 |= SAM_FLAG_PROPER_PAIRED;
        if(strand0 ==1) flag0 |= SAM_FLAG_STRAND;
        if(strand1 ==1) flag0 |= SAM_FLAG_MATE_STRAND;
        flag0 |= SAM_FLAG_MATE_READ1; 
        if(read0_mapped != 1) flag0 |= SAM_FLAG_UNMAPED;
        if(read1_mapped != 1) flag0 |= SAM_FLAG_MATE_UNMAPPED;
        printf("%u\t", flag0);//FLAG
        
        if(read0_mapped != 1) printf("*\t0\t");
        else{ 
            printf("%s\t%u\t", chrom0, pos0);//RNAME POS
        }
        //MAPQ (to be fixed)
        if(sam0->b1 == -100000 && sam0->b0 != -100000) printf("60\t");
        else printf("0\t");
        if(read0_mapped) printf("%s\t", sam0->cigar.s); 
        else printf("*\t");     
        //MRNAME MPOS
        if(read1_mapped != 1) {
            printf("*\t0\t");
        } else{
            if(read0_mapped != 1 || strcmp(chrom0, chrom1) != 0) {
                printf("%s\t%u\t", chrom1, pos1);
            } else{ 
                printf("=\t%u\t", pos1);
            }
        }
        if(read1_mapped == 1 && read0_mapped == 1){ 
            int isize = strand0==0?ABS(pos0, pos1):-ABS(pos0, pos1);
            printf("%d\t", isize);
        }else {
            printf("0\t"); 
        }
        unsigned char *s = (strand0==0)?(sam0->nst_seq):(sam0->nst_rseq);   
        for(i =0; i < sam0->l_seq; ++i) putchar("ACGT"[s[i]]);//SEQ
        putchar('\t');
        //QUAL
        s = (unsigned char *)sam0->qual; 
        if(sam0->flag & SAM_FLAG_STRAND){
            if(strand0 == 0){
                for(i =sam0->l_seq-1; i >=0; --i) putchar(s[i]);    
            } else{ printf("%s\t", s);} 
        } else{
            if(strand0 == 0){
                printf("%s\t", s);
            } else{
                for(i =sam0->l_seq-1; i >=0; --i) putchar(s[i]);
            }
        }
        putchar('\n'); 
    }
    {   
        printf("%s\t", sam0->seq_name);//QNMAE 
       
        flag1 |= SAM_FLAG_PAIRED;
        if(read1_mapped != 1) flag1 |= SAM_FLAG_UNMAPED; 
        if(read0_mapped != 1) flag1 |= SAM_FLAG_UNMAPED; 
        if(is_properpaired) flag1 |= SAM_FLAG_PROPER_PAIRED;
        if(strand1 ==1) flag1 |= SAM_FLAG_STRAND;


        if(strand0 ==1) flag1 |= SAM_FLAG_MATE_STRAND;
        flag1 |= SAM_FLAG_MATE_READ2; 
       
         
        printf("%u\t", flag1);//FLAG
        
        if(read1_mapped != 1) printf("*\t0\t");
        else{ 
            printf("%s\t%u\t", chrom1, pos1);//RNAME POS
        }
        //MAPQ (to be fixed)
        if(sam1->b1 == -100000 && sam1->b0 != -100000) printf("60\t");
        else printf("0\t");
        if(read1_mapped) printf("%s\t", sam1->cigar.s);
        else printf("*\t");//CIGAR (to be fixed)
        
        //MRNAME MPOS
        if(read0_mapped != 1) {
            printf("*\t0\t");
        } else{
            if(read1_mapped != 1 || strcmp(chrom0, chrom1) != 0) {
                printf("%s\t%u\t", chrom0, pos0);
            } else{ 
                printf("=\t%u\t", pos0);
            }
        }
        if(read1_mapped == 1 && read0_mapped == 1) {
            int isize = strand1 ==0?ABS(pos0, pos1):-ABS(pos0, pos1);
            printf("%d\t", isize);
        } else { 
            printf("0\t"); 
        }
        unsigned char *s = (strand1==0)?(sam1->nst_seq):(sam1->nst_rseq);   
        for(i =0; i < sam1->l_seq; ++i) putchar("ACGT"[s[i]]);//SEQ
        putchar('\t');
        //QUAL
        s = (unsigned char *)sam1->qual; 
        if(sam1->flag & SAM_FLAG_STRAND){
            if(strand1 == 0){
                for(i =sam1->l_seq-1; i >=0; --i) putchar(s[i]);    
            } else{ printf("%s\t", s);} 
        } else{
            if(strand1 == 0){
                printf("%s\t", s);
            } else{
                for(i =sam1->l_seq-1; i >=0; --i) putchar(s[i]);
            }
        }
   
        putchar('\n'); 
    }
}
#define MAX_SHIFT 10
#define UNMAPPED_SCORE -100000
void polish_core(char *bns_prefix, char *fn_sam, const opt_t *opt)
{

    int i, j;
    int distance_algorithm = opt->distance_algorithm;

    FILE *fp_sam = (FILE *)fopen(fn_sam, "r");
    if(fp_sam == NULL){
        fprintf(stderr, "[Error]: Can't open file %s\n", fn_sam); 
        exit(1); 
    }
    bntseq_t *bns = bns_restore(bns_prefix);  

    //map RefName to RefID
    khash_t(str) *hash=kh_init(str);
    int tid;
    for(tid =0; tid < bns->n_seqs; ++tid) {
        int ret;
        khiter_t k;
        k = kh_put(str, hash, bns->anns[tid].name, &ret); 
        kh_value(hash, k) = tid; 
    }

    uint8_t *pac = pac_restore(bns); 
    sam_skipHeader(fp_sam); 
    
    if(opt->is_paired){//PE
        //int best0_id, best0_strand, best1_id, best1_strand; 
        SAM_t *sam0 = sam_readline(fp_sam); 
        SAM_t *sam1 = sam_readline(fp_sam);
        while (sam0 != NULL && sam1 != NULL){
            int l_refseq0 = sam0->l_seq; 
            int l_refseq1 = sam1->l_seq; 
            unsigned char *refseq0= calloc((l_refseq0+15)/8*8, sizeof(char));                      
            unsigned char *refseq1= calloc((l_refseq1+15)/8*8, sizeof(char));                      
            
            for(i =0; i < 2; ++i){
                int n = sam0->multi_hits[i]->n;
                hit_t *h = sam0->multi_hits[i]->a; 
                for(j = 0; j < n; ++j){

                    //fprintf(stderr, "%s ", h[j].chrom); 
		            khiter_t k = kh_get(str, hash, h[j].chrom);

                    //fprintf(stderr, "%u\n", k); 
		            int tid = kh_value(hash, k);
                    h[j].offset = bns->anns[tid].offset+h[j].pos-1; 
                }  
                ks_introsort(hits, n, sam0->multi_hits[i]->a);//sort hits by offset
                rm_repeat_hits(sam0->multi_hits[i]);//remove repeat offset
                n = sam0->multi_hits[i]->n;
                for(j =0; j < n; ++j){
                    l_refseq0 = __get_refseq(refseq0,  l_refseq0, pac, bns->l_pac, h[j].offset);
                    /*
                    if(i==1) {//reverse
                        for(k=0; k < l_refseq0; ++k) refseq0[k] = 4-refseq0[k];
                        for(k=0; k < l_refseq0/2; ++k) SWAP(refseq0[k], refseq0[l_refseq0-1-k]); 
                    }*/
                    //recaluate editdistance
                    const int8_t *query = (i ==0)?((const int8_t*)sam0->nst_seq):((const int8_t*)sam0->nst_rseq);
                    if(distance_algorithm == ALGO_SW){
                        s_profile *p = ssw_init(query, sam0->l_seq,score_mat, 5, 1 ); 
                        s_align *result = ssw_align(p, (const int8_t *)refseq0, l_refseq0, GAP_OP, GAP_EX, 0, 150, 0, sam0->l_seq); 
                        h[j].score = result->score1;
                        init_destroy(p);
                        align_destroy(result);
                    } else if(distance_algorithm == ALGO_LV){
                        int d = computeEditDistance((const char*)refseq0, l_refseq0, (char *)query, sam0->l_seq, MAX_DISTANCE);  
                        h[j].score = (d==-1)?UNMAPPED_SCORE:-d;
      
                    } else{
                        fprintf(stderr, "distance_algorithm is not available");
                        exit(1);
                    }

                }
            
            } 
     
            for(i =0; i < 2; ++i){
                int n = sam1->multi_hits[i]->n;
                hit_t *h = sam1->multi_hits[i]->a;
                for(j = 0; j < n; ++j){

                    //fprintf(stderr, "%s ", h[j].chrom); 
                    khiter_t k = kh_get(str, hash, h[j].chrom);
                    //fprintf(stderr, "%u\n", k); 
		            int tid = kh_value(hash, k);
                    h[j].offset = bns->anns[tid].offset+h[j].pos-1; 
                }
                ks_introsort(hits, n, sam1->multi_hits[i]->a);//sort hits by offset
                rm_repeat_hits(sam1->multi_hits[i]);//remove repeat offset
                n = sam1->multi_hits[i]->n;
                for(j = 0; j <n; ++j){
                    l_refseq1 = __get_refseq(refseq1, l_refseq1, pac, bns->l_pac, h[j].offset);
                    
                    const int8_t *query = (i ==0)?((const int8_t*)sam1->nst_seq):((const int8_t*)sam1->nst_rseq);
                    /*
                    if(i==1) {//reverse complement        if(result->read_begin1 != 0){
            ksprintf(&sam->cigar, "%dS",result->read_begin1);//Cigar
        }
 
                        for(k=0; k < l_refseq1; ++k) refseq1[k] = 4-refseq1[k];
                        for(k=0; k < l_refseq1/2; ++k) SWAP(refseq1[k], refseq1[l_refseq1-1-k]); 
                    }*/
                    
                    //recalculate edit distance
                    if(distance_algorithm == ALGO_SW){
                        s_profile *p = ssw_init((const int8_t*)((i==0)?(sam1->nst_seq):(sam1->nst_rseq)), sam1->l_seq,score_mat, 5, 1 ); 
                        s_align *result = ssw_align(p, (const int8_t *)refseq1, l_refseq1, GAP_OP, GAP_EX, 0, 150, 0, sam1->l_seq); 
                        h[j].score = result->score1;
                        init_destroy(p);
                        align_destroy(result);
                    } else if(distance_algorithm == ALGO_LV){
                        int d = computeEditDistance((const char*)refseq1, l_refseq1, (char *)query, sam1->l_seq, MAX_DISTANCE);  
                        h[j].score = (d==-1)?UNMAPPED_SCORE:-d; 
                    } else{
                        fprintf(stderr, "distance_algorithm is not available");
                        exit(1);
                    }
       
                }
            }
            //Pairing
            unsigned int n_pp0 ;//proper paired number sam0 forward ,sam1 backward
            unsigned int n_pp1;//proper paired number sam0 backward, sam1 forward
            n_pp0 = __pairing(sam0->multi_hits[0], sam1->multi_hits[1]);
            n_pp1 = __pairing(sam1->multi_hits[0], sam0->multi_hits[1]);
            int is_properpaired = (n_pp0+n_pp1!=0)?1:0;
            int strand0 =-1, strand1=-1, iter0=-1, iter1=-1;
            if(is_properpaired != 1){//Not proper paired
                
                int best0 = -100000, best1 = -100000;
                for(i =  0; i < 2; ++i){
                    hit_t *h = sam0->multi_hits[i]->a;
                    int n = sam0->multi_hits[i]->n;
                    for(j = 0; j < n; ++j){
                        int score = h[j].score;
                        if(score == UNMAPPED_SCORE) continue;
                        if(score > best1){
                            best1 = score;
                            if(best1 >best0){
                                SWAP(best0, best1);
                                strand0 = i;
                                iter0 = j;

                            }
                        }
                    
                    }//END for 
                }//END for
                sam0->b0 = best0; sam0->b1 = best1; 
                
                best0 = -100000, best1 = -100000;
                for(i = 0; i < 2; ++i){
                    hit_t *h = sam1->multi_hits[i]->a;
                    int n = sam1->multi_hits[i]->n;
                    for(j = 0; j < n; ++j){
                        int score = h[j].score;
                        if(score == UNMAPPED_SCORE) continue;
                        if(score > best1){
                            best1 = score;
                            if(best1 >best0){
                                SWAP(best0, best1);
                                strand1 = i;
                                iter1 = j;

                            }
                        }
                    
                    }//END for 
                }//END for
                sam1->b0 = best0; sam1->b1 = best1;

            } else{//Pop best pair
                int best0 =-100000, best1 = -100000;
                hit_t *h0 = sam0->multi_hits[0]->a;
                hit_t *h1 = sam1->multi_hits[1]->a;

                for(i=0; i <n_pp0; ++i){
       
                    int score = h0[i].score+h1[i].score;
                    if(score == UNMAPPED_SCORE) continue;
                    if(score > best1){
                        best1 = score;
                        if(best1 >best0){
                            SWAP(best0, best1);
                            strand0 = 0;
                            strand1 = 1; 
                            iter0 = iter1 = i;
                        }
                    } 
                
                }
                h0 = sam0->multi_hits[1]->a;
                h1 = sam1->multi_hits[0]->a;
                
                for(i=0; i <n_pp1; ++i){
                    int score = h0[i].score+h1[i].score;
                    if(score == UNMAPPED_SCORE) continue;
                    if(score > best1){
                        best1 = score;
                        if(best1 > best0){
                            SWAP(best0, best1); 
                            strand0 = 1;
                            strand1 = 0; 
                            iter0 = iter1 = i;
                        }
                    } 
                
                }
                sam0->b0 = sam1->b0 = best0;
                sam0->b1 = sam1->b1 = best1;
            }
            //Output SAM format 
            sam0->strand = strand0;
            sam0->primary = iter0;
            sam1->strand = strand1;
            sam1->primary = iter1;
            if(strand0 != -1 && iter0 != -1) gen_cigar(sam0, bns->l_pac, pac, sam0->multi_hits[strand0]->a[iter0].offset, distance_algorithm);
            if(strand1 != -1 && iter1 != -1) gen_cigar(sam1, bns->l_pac, pac, sam1->multi_hits[strand1]->a[iter1].offset, distance_algorithm);
            polish_sam_pe(sam0, sam1, is_properpaired); 
            
            free(refseq0);
            free(refseq1);
            sam_destroy(sam0); 
            sam_destroy(sam1); 
            sam0 = sam_readline(fp_sam); 
            sam1 = sam_readline(fp_sam);


        } // END while 
        if(sam0) sam_destroy(sam0);
        if(sam1) sam_destroy(sam1);
    } else{//SE
        SAM_t *sam0 = sam_readline(fp_sam); 
        while (sam0 != NULL){
            int l_refseq0 = sam0->l_seq; 
            unsigned char *refseq0= calloc((l_refseq0+15)/8*8, sizeof(char));                      

            
            for(i =0; i < 2; ++i){
                int n = sam0->multi_hits[i]->n;
                hit_t *h = sam0->multi_hits[i]->a; 
                for(j = 0; j < n; ++j){
                    khiter_t k = kh_get(str, hash, h[j].chrom);
                    int tid = kh_value(hash, k);
                    h[j].offset = bns->anns[tid].offset+h[j].pos-1; 
                }  
                ks_introsort(hits, n, sam0->multi_hits[i]->a);//sort hits by offset
                rm_repeat_hits(sam0->multi_hits[i]);//remove repeat offset
                n = sam0->multi_hits[i]->n;
                for(j =0; j < n; ++j){
                    l_refseq0 = __get_refseq(refseq0,  l_refseq0, pac, bns->l_pac, h[j].offset);
                    const int8_t *query = (i ==0)?((const int8_t*)sam0->nst_seq):((const int8_t*)sam0->nst_rseq);
                    if(distance_algorithm == ALGO_SW){
                        s_profile *p = ssw_init(query, sam0->l_seq,score_mat, 5, 1 ); 
                        s_align *result = ssw_align(p, (const int8_t *)refseq0, l_refseq0, GAP_OP, GAP_EX, 0, 150, 0, sam0->l_seq); 
                        h[j].score = result->score1;
                        init_destroy(p);
                        align_destroy(result);
                    } else if(distance_algorithm == ALGO_LV){
                        int d = computeEditDistance((const char*)refseq0, l_refseq0, (char *)query, sam0->l_seq, MAX_DISTANCE);  
                        h[j].score = (d==-1)?UNMAPPED_SCORE:-d;
      
                    } else{
                        fprintf(stderr, "distance_algorithm is not available");
                        exit(1);
                    }

                }
            
            } 
     

            
            //Pairing

            int strand0 =-1, iter0=-1;
            int best0 = -100000, best1 = -100000;
            for(i =  0; i < 2; ++i){
                hit_t *h = sam0->multi_hits[i]->a;
                int n = sam0->multi_hits[i]->n;
                for(j = 0; j < n; ++j){
                    int score = h[j].score;
                    if(score == UNMAPPED_SCORE) continue;
                    if(score > best1){
                        best1 = score;
                        if(best1 >best0){
                            SWAP(best0, best1);
                            strand0 = i;
                            iter0 = j;

                        }
                    }
                
                }//END for 
            }//END for
            sam0->b0 = best0; sam0->b1 = best1; 
            


            //Output SAM format 
            sam0->strand = strand0;
            sam0->primary = iter0;


            if(strand0 != -1 && iter0 != -1) gen_cigar(sam0, bns->l_pac, pac, sam0->multi_hits[strand0]->a[iter0].offset, distance_algorithm);
            polish_sam_se(sam0); 
            free(refseq0);
            sam_destroy(sam0); 

            sam0 = sam_readline(fp_sam); 



        } // END while 
        if(sam0) sam_destroy(sam0);
    
    
    }

    bns_destroy(bns); 
    free(pac);
    fclose(fp_sam);
    kh_destroy(str, hash); 

}
int main ( int argc, char *argv[] )
{
    opt_t *opt = opt_init();
    int c;

    while((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1){
        switch(c)
        {
            case'h':
                print_help();
                return EXIT_SUCCESS;
            case's':
                opt->distance_algorithm = ALGO_SW; 
                break;
            case'p':
                opt->is_paired = 1;
                break;
            default:
                fprintf(stderr, "Unkown argument!\n");
                print_help();
                return EXIT_FAILURE;   
        
        } 
    
    
    }
   
    if(argc-optind != 2){
        print_help();
        opt_destroy(opt);
        return EXIT_SUCCESS;
    }
    char bns_prefix[1024];
    strcat(strcpy(bns_prefix, argv[optind]), ".C");
    polish_core(bns_prefix, argv[optind+1], opt); 
     
    
    opt_destroy(opt);

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
