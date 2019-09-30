/*
 * =====================================================================================
 *
 *       Filename:  test_ssw_snp.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/10/2018 10:30:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdint.h>
#include "ssw.h"
                                    //  0   A   C   G   T   N   6   7   8   9   10  11  12  13  14  15
static int8_t score_mat2[256] = {       -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* A */  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* C */ -3,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* AC*/  1,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* G */ -3, -3,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /*GA*/   1, -3,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /*GC*/  -3,  1,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /*GAC*/  1,  1,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* T */ -3, -3, -3,  1, -3, -3, -3,  -3, -3, -3, -3, -3, -3, -3, -3,-3, 
                                /* TA */ 1, -3, -3,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* TC */-3,  1, -3,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* TAC*/ 1,  1, -3,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /* TG */-3, -3,  1,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /*TGA*/  1, -3,  1,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /*TGC*/ -3,  1,  1,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3,
                                /*TGAC*/ 1,  1,  1,  1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,-3
                                        


};


int snpaln_sw_snpaware(int l_ref, uint8_t *ref, int l_seq, uint8_t *seq)
{
    int i;
    s_profile *p = ssw_init((int8_t *)seq, l_seq, score_mat2, 16, 1); 
    //fprintf(stderr, "Ref:\n");
    //for(i = 0; i != end - start +1; ++i) fprintf(stderr, "%c", "ACGT"[ref[i]]);
    //fprintf(stderr, "\nSeq:\n");
    //for(i = 0; i != l_seq; ++i) fprintf(stderr, "%c", "ACGT"[seq[i]]);
    uint8_t flag = 2;
    //int filter = 80;
    //int filtered = 0;
    int maskLen = l_seq/2;
    s_align *result = ssw_align(p, (int8_t *)ref, l_ref, 5, 2, flag, 0, 100, maskLen);//note: filterd no effect
    printf("Ref: ");
    for(i = 0; i < l_ref; ++i) printf("%2x ", ref[i]);

    printf("\nSef: ");
    for(i = 0; i < l_seq; ++i) printf("%2u ", seq[i]);

    printf("\nScore = %u, ref:%u, %u-read:%u,%u\n", result->score1, result->ref_begin1, result->ref_end1, result->read_begin1, result->read_end1);
    int j;
    for(j = 0; j < result->cigarLen; ++j) {
        uint32_t cigar = result->cigar[j];    
        printf("%u%c",cigar>>4,"MID"[cigar&15]);//Cigar
    }
    printf("\n"); 
    return 0;
}  

#include	<stdlib.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int main ( int argc, char *argv[] )
{
    uint8_t Ref[19] = {1,3,5,7, 2, 4,8,9,10,11,12,13,14,15,1,2, 4,6 ,1}; 
    uint8_t Seq[20] = {0,0,0,0, 1, 2,2,3,3, 3 , 3,3,3,   3,3, 0,1, 2, 1, 0}; 
    snpaln_sw_snpaware(19, Ref, 20, Seq);
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
