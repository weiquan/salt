#ifndef __LOOLUPTABLE_H
#define __LOOLUPTABLE_H
/*
 * =====================================================================================
 *
 *       Filename:  LookUpTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/07/2013 02:59:17 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdint.h>
typedef struct{
    uint32_t maxLookupLen;
    uint32_t n_item;
    uint32_t *item;
} lookupTable_t;

lookupTable_t *LKT_init(const int maxLookupLen);
void LKT_destroy(lookupTable_t *lkt);
lookupTable_t *LKT_restore(const char *fn_lkt);
uint32_t LKT_seq2LktItem(const uint8_t *seq, int from, int to);
void LKT_build_lookuptable(const char *fn_pac, lookupTable_t *lkt);
static inline void LKT_lookup_sa(lookupTable_t *lkt, const uint8_t *seq, int from, int to, uint32_t *l, uint32_t *k)
{
    uint32_t i;
    
    i = LKT_seq2LktItem(seq, from, to);
    *l = lkt->item[i];
    *k = lkt->item[i+1]-1;
}

#endif
