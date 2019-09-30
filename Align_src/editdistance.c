/*
 * =====================================================================================
 *
 *       Filename:  seedex.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/26/2012 02:35:10 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * ====================================================i=================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "editdistance.h"
#define BITS_PER_CHAR 4
#define CHAR_PER_BYTE 2
#define CHAR_PER_WORD 8
#define MASK_TAIL(n) (0xFFFFFFFF>>(32-((n)<<2)))
#define MASK_HEAD(n) (0xFFFFFFFF<<(32-((n)<<2)))
#define __set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
static inline uint32_t __POPCOUNT(uint32_t v)
{
    v = v-((v>>1)&0x55555555);
    v = (v&0x33333333)+((v>>2)&0x33333333);
    return ((v+(v>>4)&0xF0F0F0F)*0x1010101)>>24;
}



#define CHECK_MATCH_CHAR(word, nt, char_index) ((word) & ((nt) << (BITS_PER_CHAR * (char_index)  ) ) ) 
static const uint8_t nt2bit[5] = {1,2,4,8, 15}; 


#define MATCH 1
#define UN_MATCH 0

static inline int exact_match_num(uint32_t a, uint32_t b)
{
    return __POPCOUNT(a&b);
}
int ed_mismatch_alt(const uint32_t *mixRef, uint32_t ref_st, const uint8_t *seq, uint32_t l_comp, int max_err)
{
    int n_err = 0;
    uint32_t i, j, n_aln = 0;
    uint32_t word_st, word_ed, char_firstword, char_lastword;
    word_st = ref_st/CHAR_PER_WORD;
    char_firstword = ref_st%CHAR_PER_WORD;// char_firstword chars on the right of ref_st 
    word_ed = (ref_st+l_comp-1)/CHAR_PER_WORD;
    char_lastword = (ref_st+l_comp-1)%CHAR_PER_WORD;

    uint32_t wordpacked_seq = 0;/* word packed seq */
    /*count mismatch in first word*/
    for(i =0; i< CHAR_PER_WORD- char_firstword; ++i){ 
        wordpacked_seq |= (uint32_t)nt2bit[seq[n_aln++]]<<(i*BITS_PER_CHAR); 
    }
    n_err += char_lastword - exact_match_num(wordpacked_seq, mixRef[word_st]>>(BITS_PER_CHAR*char_firstword)); 
    if(n_err > max_err) return -1;
    /* core loop */ 
    for(j = word_st+1; j < word_ed; ++j){
        wordpacked_seq = 0;
        for(i = 0; i < CHAR_PER_WORD; ++i) {
            wordpacked_seq |= (uint32_t)nt2bit[seq[n_aln++]]<<(i*BITS_PER_CHAR);
        } 
        n_err += CHAR_PER_WORD -exact_match_num(wordpacked_seq, mixRef[j]); 
        if(n_err > max_err) return -1;
    }
    /* count mismatch in last word */
    wordpacked_seq = 0;
    for(i =0; i< char_lastword; ++i){ 
        wordpacked_seq |= (uint32_t)nt2bit[seq[n_aln++]]<<(i*BITS_PER_CHAR); 
    }
    n_err += char_lastword- exact_match_num(wordpacked_seq, mixRef[word_ed]&MASK_TAIL(char_lastword)); 

    if(n_err > max_err) return -1;


    return n_err;
}
int ed_mismatch(const uint32_t *mixRef, uint32_t ref_st, const uint8_t *seq, uint32_t l_comp, int max_err)
{
	int n_err = 0;
	int i, seq_iter= 0;
	uint32_t word_st, word_ed;
	uint32_t char_in_first_word, char_in_last_word;

	word_st = ref_st/CHAR_PER_WORD;
	word_ed = (ref_st+l_comp-1)/CHAR_PER_WORD;
	char_in_first_word = ref_st % CHAR_PER_WORD;
	char_in_last_word = (ref_st+l_comp-1)%CHAR_PER_WORD;

	uint32_t word; 
#ifdef DEBUG
    //fprintf(stderr, "seedex!\n");	
#endif
	word = mixRef[word_st];
	word >>= char_in_first_word*BITS_PER_CHAR;
    //fprintf(stderr, "\ncheck pos%u max_diff%d\n", ref_st, max_err);
    for(i = char_in_first_word; i < CHAR_PER_WORD; ++i){
/*  
        fprintf(stderr, "[ref:read]: %u\t%u\n", word&15, nt2bit[seq[seq_iter]]); 
        if(( word&15 & nt2bit[seq[seq_iter]] ) == 0 ){
            fprintf(stderr, "!="); 
        }
*/
        if( ( word&15 & nt2bit[seq[seq_iter]] ) == 0){
            if(++n_err > max_err){
                return -1;
            }

		}         
        ++seq_iter;
        word>>=BITS_PER_CHAR;
	}
	++word_st;
	
	for(; word_st != word_ed; ++word_st){
		word = mixRef[word_st];
		for(i = 0; i < CHAR_PER_WORD; ++i){
/*  
            fprintf(stderr, "[ref:read]: %u\t%u\n", word&15, nt2bit[seq[seq_iter]]); 
            if(( word&15 & nt2bit[seq[seq_iter]] ) == 0 ){
                fprintf(stderr, "!="); 
            }
*/
            if(( word&15 & nt2bit[seq[seq_iter]] ) == 0 ){
                if(++n_err > max_err){
                    return -1;
                }
		    }             
            ++seq_iter;
    	    word >>= BITS_PER_CHAR;
		}
	}

	word = mixRef[word_ed];
	for(i = 0; i <= char_in_last_word; ++i){
#ifdef DEBUG
        fprintf(stderr, "[ref:read]: %u\t%u\n", word&15, nt2bit[seq[seq_iter]]); 
        if(( word&15 & nt2bit[seq[seq_iter]] ) == 0 ){
            fprintf(stderr, "!="); 
        }
#endif
 	    if( ( word&15 & nt2bit[seq[seq_iter]] ) == 0){
            if(++n_err >max_err){
                return -1;
            }
		}         
        ++seq_iter;
    
		word>>=BITS_PER_CHAR;
	}
//	++word_st;
	return n_err;
}
int ed_mismatch_2bitref(const uint32_t *pac, uint32_t ref_st, const uint8_t *seq, int l_comp, int max_err)
{
    int i, n_err = 0;
    for(i = 0; i < l_comp; ++i){
        if(__get_pac(pac, i+ref_st) != seq[i] && ++n_err>max_err) return -1;  
    }	
    return n_err;
}


int ed_diff(const uint32_t *mixRef, uint32_t l_mref, uint32_t ref_st, const uint32_t l_ref, const uint8_t *seq, uint32_t l_seq, int max_k_diff)
{
    //fprintf(stderr , "\n [ed_diff] : %u %u %d\n", l_ref, l_seq, max_k_diff );
    //if(ref_st>l_ref || ref_st+l_ref > l_mref) return -1;
    if(ref_st>l_mref || ref_st+l_ref > l_mref) return -1;
    #define __is_even(n) (n&1 == 0)
    #define __is_odd(n) (n&1 == 1)
    #define __get_char_in_word(word, pos) (((word)>>(BITS_PER_CHAR*(pos))) & 0xF) 
    uint32_t i;
    uint8_t *t = (uint8_t *)calloc((l_ref+15)/8*8, sizeof(uint8_t));
    uint8_t *q = (uint8_t *)calloc((l_seq+15)/8*8, sizeof(uint8_t));
    uint32_t ti = 0, qi = 0;
    
    /* Init text */
    //uint8_t c;
    uint32_t word, wi = ref_st/CHAR_PER_WORD;
    //first word
    if(ref_st%CHAR_PER_WORD != 0){
        word = mixRef[wi];
        for(i = ref_st%CHAR_PER_WORD; i < CHAR_PER_WORD; ++i) {
            t[ti++] = __get_char_in_word(word, i);            
        }
        ++wi;
    }
    //mid words(
    uint32_t n_words = (l_ref+ref_st)/CHAR_PER_WORD;
    for(; wi < n_words; ++wi ) {
        word = mixRef[wi];
        for(i = 0; i < CHAR_PER_WORD; ++i){
            t[ti++] = __get_char_in_word(word, i);            
        } 
    } 
        //last word
    if(ti != l_ref){
        word = mixRef[wi]; 
        for(i = 0; ti <l_ref ; ++i) {
            t[ti++] = __get_char_in_word(word, i);            
        }
    }       
     
    /* Init query */
    for(i = 0; i < l_seq;++i) {
        if(seq[i] >3 ) q[qi++] = 15;   
        else q[qi++] = 1<<seq[i];             
    }
   
/*
    fprintf(stderr, "check %u\n", ref_st);
    uint32_t __i=0, __j=0;
    for(__i=0; __i< 20; ++__i) fprintf(stderr, "%u\t",t[__i]);  
    for(__j=0; __j< 20; ++__j) fprintf(stderr, "%u\t", q[__j]);  
*/
    int r =  computeEditDistance((char *)t, ti, (char *)q, qi, max_k_diff); 
            
    free(t);
    free(q);
    return r;

}

int ed_diff_withcigar(const uint32_t *mixRef, uint32_t ref_st, const uint32_t l_ref, const uint8_t *seq, uint32_t l_seq, int max_k_diff, char *cigarBuf, int cigarLen, int useM, CigarFormat cigarFormat)
{
    #define __is_even(n) (n&1 == 0)
    #define __is_odd(n) (n&1 == 1)
    #define __get_char_in_word(word, pos) ((word>>(BITS_PER_CHAR*pos)) & 0xF) 
    uint32_t i;
    uint8_t *t = (uint8_t *)calloc((l_ref+15)/8*8, sizeof(uint8_t));
    uint8_t *q = (uint8_t *)calloc((l_seq+15)/8*8, sizeof(uint8_t));
    uint32_t ti = 0, qi = 0;
    
    /* Init text */
    //uint8_t c;
    uint32_t word, wi = ref_st/CHAR_PER_WORD;
    //first word
    if(ref_st%CHAR_PER_WORD != 0){
        word = mixRef[wi];
        for(i = ref_st%CHAR_PER_WORD; i < CHAR_PER_WORD; ++i) {
            t[ti++] = __get_char_in_word(word, i);            
        }
        ++wi;
    }
        //mid words(
    uint32_t n_words = (l_ref+ref_st)/CHAR_PER_WORD;
    for(; wi < n_words; ++wi ) {
        word = mixRef[wi];
        for(i = 0; i < CHAR_PER_WORD; ++i){
            t[ti++] = __get_char_in_word(word, i);            
        } 
    } 
        //last word
    if(ti != l_ref){
        word = mixRef[wi]; 
        for(i = 0; ti <l_ref ; ++i) {
            t[ti++] = __get_char_in_word(word, i);            
        }
    }       
     
    /* Init query */
    for(i = 0; i < l_seq;++i) { 
        if(seq[i] >3 ) q[qi++] = 15;   
        else q[qi++] = 1<<seq[i];
    }
   
    /* */
    int ret =  computeEditDistanceWithCigar((char *)t, ti, (char *)q, qi, max_k_diff, cigarBuf, cigarLen, useM, cigarFormat); 
                
    free(t);
    free(q);
    return ret;

}
