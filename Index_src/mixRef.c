/*
 * =====================================================================================
 *
 *       Filename:  mixRef.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/09/2012 01:25:07 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#include "kseq.h"
#include "hapmap.h"
#include "mixRef.h"
#define CHAR_IN_WORD 8
/*
#if UNIT_TESTING
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>
#endif 
*/
#include<stdint.h>
KSEQ_INIT(gzFile, gzread)
static uint8_t nt5_4bit_table[256] = {
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0 /*'-'*/, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};// {A,C,G,T,#, N} = {0, 1, 2, 3, 4, 5}

#define CHAR_PER_BYTE 2
#define BYTE_PER_WORD 4
#define CHAR_PER_WORD 8
mixRef_t *mixRef_restore(const char *fn)
{
   
    FILE *fp;
    fp = fopen(fn,"rb");
    if(fp == NULL){
        fprintf(stderr, "[mixRef_restore]:file open %s fail!\n", fn);
        exit(EXIT_FAILURE);
    }


    mixRef_t *mixRef = (mixRef_t *)calloc(1, sizeof(mixRef_t));
    if(mixRef == NULL){
        fprintf(stderr, "[mixRef]allocate mem fail!\n");
        exit(EXIT_FAILURE);
    } 

    fread(&mixRef->l, sizeof(uint32_t), 1, fp);   
    mixRef->seq = (uint32_t *)calloc((mixRef->l+CHAR_IN_WORD-1)/CHAR_IN_WORD, sizeof(uint32_t)); 
    if(mixRef->seq == NULL){
        fprintf(stderr, "[mixRef_restore]: mixRef->seq allocate mem fail!\n");
        exit(EXIT_FAILURE);
    }
    fread(mixRef->seq, sizeof(uint32_t),(mixRef->l + CHAR_IN_WORD -1)/CHAR_IN_WORD , fp); 
    fclose(fp);
    return mixRef;
}
void mixRef_destroy(mixRef_t *mixRef)
{
    free(mixRef->seq);
    free(mixRef);
}
#define __set_pac(pac, l, v) (pac)[(l)>>3]|=((v)<<4*((l)%8))
#define __clear_pac(pac, l) (pac)[(l)>>3]&=(unsigned int)(~(15<<(4*((l)%8))))
#define __get_pac(pac, l) (((pac)[(l)>>3]>>4*((l)%8))&15)
int build_mixRef(const char *fn_fa, const char *fn_hapmap, const char *fn_mixRef)
{
    int i;    
    FILE *fp_hapmap =NULL, *fp_mixRef = NULL;
    gzFile fp_fa = Z_NULL;
    kseq_t *seq; hapmap_t *hm;
    mixRef_t *mixRef;
    {//open stream
        fp_fa	= gzopen( fn_fa, "r" );
        if ( fp_fa == NULL ) {
            fprintf ( stderr, "couldn't open file '%s'; %s\n",
                fn_fa, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        fp_hapmap	= fopen( fn_hapmap, "r" );
        if ( fp_hapmap == NULL ) {
            fprintf ( stderr, "couldn't open file '%s'; %s\n",
                    fn_hapmap, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        fp_mixRef	= fopen( fn_mixRef, "w" );
        if ( fp_hapmap == NULL ) {
            fprintf ( stderr, "couldn't open file '%s'; %s\n",
                    fn_mixRef, strerror(errno) );
            exit (EXIT_FAILURE);
        }
    }
    
    mixRef = (mixRef_t *)calloc(1, sizeof(mixRef_t));
    if(mixRef == NULL){
        fprintf(stderr, "[main_gen_mixRef]allocate mem fail!\n");
        exit(EXIT_FAILURE);
    }
    seq = kseq_init(fp_fa);
    hm = hapmap_init(fp_hapmap); 
  
    uint32_t tot_l = 0;int l=0;
    mixRef->seq = NULL;
    while( (l = kseq_read(seq)) >0)
    {
        uint32_t __i;
        uint32_t n_word = (l+tot_l-1+CHAR_IN_WORD)/CHAR_IN_WORD;
        mixRef->seq =   (uint32_t *)realloc(mixRef->seq, n_word*sizeof(uint32_t));     
        if(mixRef->seq  == NULL){
            fprintf(stderr, "[main_gen_mixRef]:realloc mem fail!\n");
            exit(EXIT_FAILURE);
        }
        uint32_t *rpac = mixRef->seq;
        char *q = seq->seq.s;
        for(__i = tot_l; __i < tot_l+l; ++__i){
            __clear_pac(rpac, __i); 
            __set_pac(rpac, __i, nt5_4bit_table[(uint8_t)q[__i-tot_l]]); 
        }

        hapmap_readhm(hm);
        //fprintf(stderr, "snppos:%u\n", hm->snp_pos[0]);
        for(__i=0; __i < hm->snp_num; ++__i){
             __set_pac(mixRef->seq, tot_l+hm->snp_pos[__i], hm->snp_type[__i] &15); 
        }        
/*
        for(__i =tot_l; __i < tot_l+l;++__i){
            uint8_t c = __get_pac(mixRef->seq, __i);
            putchar('\t');
            putchar('[');

            for(i=0; i < 4; ++i){
                if((c&1)!= 0) putchar("ACGT"[i]);
                c >>=1;
            }
            putchar(']'); 
            putchar('\t'); 
        }
*/
        tot_l += l;    
    } 
    mixRef->l = tot_l;
    fwrite(&mixRef->l, sizeof(uint32_t), 1, fp_mixRef);   
    fwrite(mixRef->seq, sizeof(uint32_t), (mixRef->l-1+CHAR_IN_WORD)/CHAR_IN_WORD, fp_mixRef);
    
   
    free(mixRef->seq);
    free(mixRef);
    kseq_destroy(seq);
    hapmap_destroy(hm); 
    
    
    {//close stream
        if( gzclose(fp_fa) == EOF ) {			/* close output file   */
            fprintf ( stderr, "couldn't close file '%s'; %s\n",
                fn_fa, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        if( fclose(fp_hapmap) == EOF ) {			/* close input file   */
            fprintf ( stderr, "couldn't close file '%s'; %s\n",
                    fn_hapmap, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        if( fclose(fp_mixRef) == EOF ) {			/* close input file   */
            fprintf ( stderr, "couldn't close file '%s'; %s\n",
                    fn_hapmap, strerror(errno) );
            exit (EXIT_FAILURE);
        }
    }
    return EXIT_SUCCESS;
}
#ifdef MAIN_GEN_MIXREF
static int usage()
{
    fprintf(stderr, "********************\n");
    fprintf(stderr, "gm <ref.fa.in> <hapmap.in><mixRef.out>\n");
     
    fprintf(stderr, "********************\n");
    return 1;
}

int main(int argc, char *argv[])
{
    if(argc < 4){
        return usage();
    }
    
    char *fn_fa = argv[1];
    char *fn_hapmap = argv[2];
    char *fn_mixRef = argv[3];
    return build_mixRef(fn_fa, fn_hapmap, fn_mixRef); 
}
#endif 
