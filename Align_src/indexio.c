/*
 * =====================================================================================
 *
 *       Filename:  indexio.c
 *
 *    Description:  :
 *
 *        Version:  1.0
 *        Created:  09/11/2014 05:37:28 PM
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
#include "indexio.h"

index_t* alnse_index_reload(  const char *prefix)
{
    char cbwt_prefix[128], rbwt_prefix[128], fn_local_pattern[128], fn_mixRef[128], fn_lkt[128];  
    index_t *index = NULL;
    
    index = calloc(1, sizeof(index_t));
    if(index == NULL){
        fprintf(stderr, "[aln_index_reload]: allocate mem fail!\n");
        exit(1);
    }
    
    memcpy(cbwt_prefix, prefix, strlen(prefix) +1); strcat(cbwt_prefix, ".C");
    memcpy(rbwt_prefix, prefix, strlen(prefix) +1); strcat(rbwt_prefix, ".R");
    memcpy(fn_local_pattern, prefix, strlen(prefix) +1); strcat(fn_local_pattern, ".lp");
    memcpy(fn_lkt, prefix, strlen(prefix) +1); strcat(fn_lkt, ".C.lkt");
    memcpy(fn_mixRef, prefix, strlen(prefix) +1); strcat(fn_mixRef, ".ref");
    
    fprintf(stderr, "[aln_index_reload]: reload C part index!\n");
    index->cbwt = cbwt_init(cbwt_prefix);
    fprintf(stderr, "[aln_index_reload]: reload R part index!\n");
    index->rbwt2 = Rbwt2_init(rbwt_prefix);
    fprintf(stderr, "[aln_index_reload]: reload mixRef!\n");
    index->mixRef = mixRef_restore(fn_mixRef); 
    index->bntseq = bns_restore(cbwt_prefix);  
    index->pac = pac_restore(index->bntseq);
    index->lkt = LKT_restore(fn_lkt);
    return index;
}
void alnse_index_destroy(index_t *index)
{
    bwt_destroy(index->cbwt); 
    Rbwt2_destroy(index->rbwt2);
    mixRef_destroy(index->mixRef);
    bns_destroy(index->bntseq);
    free(index->pac);
    LKT_destroy(index->lkt);
    free(index);
}

index_t* alnpe_index_reload(const char *prefix)
{
    char cbwt_prefix[128], rbwt_prefix[128], fn_local_pattern[128], fn_mixRef[128], fn_lkt[128];  
    index_t *index = NULL;
    
    index = calloc(1, sizeof(index_t));
    if(index == NULL){
        fprintf(stderr, "[aln_index_reload]:    allocate mem fail!\n");
        exit(1);
    }
    
    memcpy(cbwt_prefix, prefix, strlen(prefix) +1); strcat(cbwt_prefix, ".C");
    memcpy(rbwt_prefix, prefix, strlen(prefix) +1); strcat(rbwt_prefix, ".R");
    memcpy(fn_local_pattern, prefix, strlen(prefix) +1); strcat(fn_local_pattern, ".lp");
    memcpy(fn_lkt, prefix, strlen(prefix) +1); strcat(fn_lkt, ".C.lkt");
    memcpy(fn_mixRef, prefix, strlen(prefix) +1); strcat(fn_mixRef, ".ref");
  
    fprintf(stderr, "[aln_index_reload]:    reload C part index!\n");
    index->cbwt = cbwt_init(cbwt_prefix);
    fprintf(stderr, "[aln_index_reload]:    reload R part index!\n");
    index->rbwt2 = Rbwt2_init(rbwt_prefix);
    fprintf(stderr, "[aln_index_reload]:    reload mixRef!\n");
    index->mixRef = mixRef_restore(fn_mixRef); 
    index->bntseq = bns_restore(cbwt_prefix);  
    index->pac = pac_restore(index->bntseq);
    index->lkt = LKT_restore(fn_lkt);
    return index;
}
void alnpe_index_destroy(index_t *index)
{
    bwt_destroy(index->cbwt); 
    Rbwt2_destroy(index->rbwt2);
    mixRef_destroy(index->mixRef);
    bns_destroy(index->bntseq);
    free(index->pac);
    LKT_destroy(index->lkt);
    free(index);
}



