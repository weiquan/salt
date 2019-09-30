/*
 * =====================================================================================
 *
 *       Filename:  samParser.h
 *
 *    
 *
 *        Version:  1.0
 *        Created:  06/27/2014 04:18:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#ifndef __SAMPARSER_H
#define __SAMPARSER_H
#include <stdlib.h>
#include <stdio.h>
#include "kvec.h"
#include "kstring.h"
/********************************************SAM struct***************************************/
typedef struct{
    char *chrom;
    unsigned int pos;
    unsigned int offset;
    int score;
}hit_t;
typedef kvec_t(hit_t) multi_hits_t;
typedef struct{
    int l_seq_name; char *seq_name;
    int flag;
    int strand;
    int primary;    //multi_hits->a[primary].tid, multi_hits->a[primary].pos
    int mapq;
//    int n_cigar; char *cigar;
    kstring_t cigar;
    char *mate_seqname;
    unsigned int mate_pos;
    int isize;
    int l_seq; char *seq; //neclotide sequence 
    unsigned char *nst_seq, *nst_rseq;//number string text
    char *qual;
    multi_hits_t *multi_hits[2]; 
    
    int b0,b1; //--LV only
        //--SW only
    int l_md; char *MD;
    char *buf;    
} SAM_t;

/********************************************SAM struct end***************************************/


SAM_t *sam_readline(FILE *fp);
void sam_destroy(SAM_t *sam);
void sam_skipHeader(FILE *fp);




#endif
