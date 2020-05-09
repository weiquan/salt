/*
 * =====================================================================================
 *
 *       Filename:  index1.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/11/2013 02:29:22 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <time.h>
#include "bntseq.h"
#include "LookUpTable.h"
#include "bwt.h"
#include "rbwt.h"
#include "localPattern.h"
#include "mixRef.h"
static int index_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "index  [OPT]  <ref.fa>  <variants.txt>  <prefix>\n\n"); 
    fprintf(stderr, "       -k <int>         seedLen [25]\n"); 
    fprintf(stderr, "       -h\n"); 
    fprintf(stderr, "\n");
    return 1;
}
//file name
//  index
//  index.C.pac index.C.rpac index.C.bwt index.C.occ index.C.sa
//  index.R.rpac index.R.for.bwt index.R.for.occ index.R.for.sa 
//  index.R.pac index.R.back.bwt index.R.back.occ index.R.back.sa 
//  index.lp index.ref
int C_sa_intv = 8;
int maxLookupLen = 12;
int index_main(int argc, char *argv[])
{
    clock_t t;
    int l_seed = 25;
    int c;
    while((c = getopt(argc, argv, "k:h"))>=0){
        switch(c){
            case 'k':
                l_seed = atoi(optarg); 
                break;
            case 'h':
                return index_usage();
            default:
                return 1;
              
        }
    
    }
    if (argc - optind  != 3 ){
        return index_usage();
    }
    
    char str0[1024], str1[1024], str2[1024];
    char *fn_fa, *fn_hm, *prefix;
    fn_fa = argv[optind]; fn_hm = argv[optind+1]; prefix = argv[optind+2];    

    /*************index C part*********************************************/ 

                            //fa2pac
    //str0 = "prefix.C"
    strncpy(str0, prefix, 1024);strcat(str0, ".C");
    
    t = clock();
    gzFile fp_fa;
    fp_fa = gzopen(fn_fa, "r");
    fprintf(stderr, "[salt_index]: convert %s file to packedfile!\n",str0); 
    bns_fasta2bntseq(fp_fa, str0);
    gzclose(fp_fa);
    fprintf(stderr, "[salt_index]: C part fa2pac %.2f seconds elapse.\n", 
                                     (float)(clock() - t) / CLOCKS_PER_SEC);

                            //LookUpTable
    //str0 = "prefix.C.pac" str1 = "prefix.C.lkt"
    strncpy(str0, prefix, 1024);strcat(str0, ".C.pac");
    strncpy(str1, prefix, 1024);strcat(str1, ".C.lkt");
    
    lookupTable_t *lkt;
    lkt = LKT_init(maxLookupLen);    
    fprintf(stderr, "[salt_index]: build lookup table from file %s\n", str0);
    LKT_build_lookuptable(str0, lkt);
    fprintf(stderr, "[salt_index]: build lookup table finish!\n");
    LKT_output(str1, lkt);
    fprintf(stderr, "[salt_index]: output lookup table to %s\n", str1);
    LKT_destroy(lkt);

                            //pac2bwt
    //str0 = "prefix.C.pac"; str1 = "prefix.C.bwt"; str2 = "prefix.C.sa"
    strncpy(str0, prefix, 1024);strcat(str0, ".C.pac");
    strncpy(str1, prefix, 1024);strcat(str1, ".C.bwt");
    strncpy(str2, prefix, 1024);strcat(str2, ".C.sa");
    
    fprintf(stderr, "[salt_index]: convert C pac%s to bwt%s\n", str0, str1);
    t = clock(); 
    bwt_bwtgen(str0, str1); 
    fprintf(stderr, "[salt_index]: C part pac2bwt %.2f seconds elapse.\n", 
                                      (float)(clock() - t) / CLOCKS_PER_SEC);
                            //update bwt and call sa
    {
        bwt_t *bwt;
        t = clock();
        fprintf(stderr, "[salt_index]: C part Update BWT... ");
        bwt = bwt_restore_bwt(str1);
        bwt_bwtupdate_core(bwt);
        bwt_dump_bwt(str1, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        bwt_t *bwt;
        t = clock();
        fprintf(stderr, "[salt_index]: C part Construct SA from BWT and Occ... ");
        bwt = bwt_restore_bwt(str1);
        bwt_cal_sa(bwt, C_sa_intv);
        bwt_dump_sa(str2, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    }


    /*************index R part****************************/

    
    strncpy(str0, prefix, 1024);strcat(str0, ".R.seedLen");
    FILE *fp = fopen(str0, "w");
    fwrite(&l_seed,sizeof(int), 1,  fp);
    fclose(fp); 

    
    //str0 = "prefix.lp"; str1 = "prefix.R"
     
    strncpy(str0, prefix, 1024);strcat(str0, ".lp");
    strncpy(str1, prefix, 1024);strcat(str1, ".R");
                            //gen local pattern

    fprintf(stderr, "[salt_index]: Build local pattern %s\n", str0);
    //ss_core(fn_fa, fn_hm, l_seed, str0);
    ss_core_alt(fn_fa, fn_hm, l_seed, str0);

                            //lp to pac
    gzFile fp_lp;
    fp_lp = gzopen(str0, "r");
    fprintf(stderr, "[salt_index]: Convert local pattern to pac with prefix %s\n", str1);
    R_bns_fasta2bntseq(fp_lp, str1);
    gzclose(fp_lp);
    fprintf(stderr, "[salt_index]: R part fa2pac %.2f seconds elapse.\n", 
                                     (float)(clock() - t) / CLOCKS_PER_SEC);
                            //pac to bwt 
    //str1 = "prefix.R.rpac"; str2="prefix.R.forward"
    strncpy(str1, prefix, 1024);strcat(str1, ".R.rpac");
    strncpy(str2, prefix, 1024);strcat(str2, ".R.forward");
    fprintf(stderr, "[salt_index]: convert R pac %s to R bwt %s\n", str1, str2); 
    Rbwt_bwt_bwtgen(str1, str2);
    //str1 = "prefix.R.pac"; str2="prefix.R.backward"
    strncpy(str1, prefix, 1024);strcat(str1, ".R.pac");
    strncpy(str2, prefix, 1024);strcat(str2, ".R.backward");
    fprintf(stderr, "[salt_index]: convert R pac %s to R bwt %s\n", str1, str2); 
    Rbwt_bwt_bwtgen(str1, str2);
    //str1 = "prefix.R" str0 = "prefix.lp" 

    strncpy(str1, prefix, 1024);strcat(str1, ".R");
    fprintf(stderr, "[salt_index]: build R sa from local pattern %s\n", str0); 
    Rbwt2_save_sa(str1, str0);
    /*************index mixref****************************/
    //str0 = "prefix.ref"
    strncpy(str0, prefix, 1024);strcat(str0, ".ref");
    fprintf(stderr, "[salt_index]: build mixRef file %s\n", str0); 
    build_mixRef(fn_fa, fn_hm, str0); 

    return 0;
}
#ifdef MAIN_INDEX

int main ( int argc, char *argv[] )
{
    index_main(argc, argv);
    return 0;
}				/* ----------  end of function main  ---------- */

#endif
