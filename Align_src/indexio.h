/*
 * =====================================================================================
 *
 *       Filename:  indexio.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/11/2014 05:36:04 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */

#include <stdio.h>
#include "rbwt.h"
#include "bwt.h"
#include "metaref.h"
#include "bntseq.h"
#include "lookup.h"
#define cbwt_t bwt_t//from bwa
typedef struct{
    rbwt2_t *rbwt2;
    cbwt_t *cbwt;
    mixRef_t *mixRef;
    bntseq_t *bntseq; 
    uint8_t *pac;
    lookupTable_t *lkt;
} index_t;
static inline cbwt_t *cbwt_init(const char *prefix)
{
    char str[1024];
    cbwt_t *bwt; 
    
    strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
    strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
    
    return bwt;    
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
index_t* alnse_index_reload(  const char *prefix);
void alnse_index_destroy(index_t *index);

index_t* alnpe_index_reload(const char *prefix);
void alnpe_index_destroy(index_t *index);

