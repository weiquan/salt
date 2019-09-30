/*
 * =====================================================================================
 *
 *       Filename:  LookUpTable.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/07/2013 02:59:12 PM
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
#include <assert.h>
#include "LookUpTable.h"

//from bwa 0.6.1
#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#define POWER_OF_2(a) (1<<a*2)

lookupTable_t *LKT_init(const int maxLookupLen)
{
    uint32_t n_item;
    n_item = POWER_OF_2(maxLookupLen) +1;
    
    lookupTable_t *lkt;
    lkt = (lookupTable_t *)malloc(sizeof(lookupTable_t)+n_item*sizeof(uint32_t));
    lkt->maxLookupLen = maxLookupLen;
    lkt->n_item = n_item;
    lkt->item = (uint32_t *)(lkt+1); 
    memset(lkt->item, 0, lkt->n_item*sizeof(uint32_t));
    return lkt;
}
void LKT_destroy(lookupTable_t *lkt)
{
    free(lkt);
//    free(lkt->item);
}
lookupTable_t *LKT_restore(const char *fn_lkt)
{
    //lookupTable_t lkt;
    int maxLookupLen;    
    FILE * fp_lkt;

    fp_lkt = fopen(fn_lkt, "rb");
    if(fp_lkt == NULL){
        fprintf(stderr, "[build_lookuptable]:file %s open fail!\n", fn_lkt);
        exit(1);
    }
    fread(&maxLookupLen, sizeof(int), 1, fp_lkt);

    lookupTable_t *lkt;
    lkt = LKT_init(maxLookupLen);
    fread(lkt->item, sizeof(uint32_t), lkt->n_item, fp_lkt);
    
    return lkt;
}
void LKT_build_lookuptable(const char *fn_pac, lookupTable_t *lkt)
{
    
    uint8_t *pac;
    uint32_t l_ref, l_pac;
    uint8_t char_in_last_byte;

    //load pac 
    FILE *fp_pac = fopen(fn_pac, "rb");
    if(fp_pac == NULL){
        fprintf(stderr, "[build_lookuptable]:file %s open fail!\n", fn_pac);
        exit(1);
    }
    fseek(fp_pac, -1, SEEK_END);
    fread(&char_in_last_byte,1, 1, fp_pac); 
    l_pac = ftell(fp_pac) -1;
    l_ref = (l_pac -1)*4 + char_in_last_byte;
    fseek(fp_pac, 0, SEEK_SET);
    pac = calloc(l_pac, sizeof(uint8_t));
    if(pac ==NULL){
        fprintf(stderr, "[build_lookuptable]:mem allocate fail!\n");
        exit(1);
    }
    fread(pac, l_pac, sizeof(uint8_t), fp_pac);
    fclose(fp_pac);

    //gen lookup table
    uint32_t i = 0, j;
    uint32_t iter_item;
    uint32_t mask = lkt->n_item -2;
    uint32_t *item = lkt->item;
    int LookUpSeqLen = lkt->maxLookupLen; 
    if(l_ref < LookUpSeqLen || LookUpSeqLen > 32){
        fprintf(stderr, "[build_lookuptable]:length of reference is shorter than lookup max length!\n");
        exit(1);
    }
    iter_item = 0;
    for(j = 0; j < LookUpSeqLen; ++j)
    {
        

	iter_item <<= 2;
        iter_item &= mask; 
        iter_item |=  _get_pac(pac, i + j);
    }
    assert(iter_item < lkt->n_item);
    ++item[iter_item+1];
    ++i;
    for(; i <= l_ref - LookUpSeqLen; ++i)
    {
        if(i%1000000 == 0)
	fprintf(stderr,"iter :%u\n", i );
        for(j = 0; j < LookUpSeqLen; ++j)
        {

            iter_item <<= 2;
            iter_item &= mask; 
            iter_item |=  _get_pac(pac, i + j);
        }
	assert(iter_item < lkt->n_item);
        ++item[iter_item+1];

    }
    for(j = 0; j < LookUpSeqLen; ++j)
    {

        iter_item <<= 2;
        iter_item &= mask; 
	++item[iter_item+1];
    }
    /*
    for(i = 0; i < lkt->n_item; ++i)
    { 
        fprintf(stderr, "%u\n", item[i]);
    }*/
    //accumulate number
    for(i = 1; i < lkt->n_item; ++i)
    {
        item[i] += item[i -1];
    }
    
    free(pac);
}
void LKT_output(const char *fn_lkt, lookupTable_t *lkt)
{
    FILE * fp_lkt;
    fp_lkt = fopen(fn_lkt, "wb");
    if(fp_lkt == NULL){
        fprintf(stderr, "[build_lookuptable]:file %s open fail!\n", fn_lkt);
        exit(1);
    }
    //fwrite(lkt, 1, sizeof(lookupTable_t)+lkt->n_item * sizeof(uint32_t), fp_lkt);
    fwrite(&lkt->maxLookupLen, sizeof(int), 1,fp_lkt);
    fwrite(lkt->item, sizeof(uint32_t), lkt->n_item, fp_lkt);
    fclose(fp_lkt);

}
uint32_t LKT_seq2LktItem(const uint8_t *seq, int from, int to)
{
    int i;
    uint32_t iter_item = 0;

    for(i = from; i < to; ++i){
        iter_item |= seq[i];
        iter_item <<= 2;
    }

    iter_item |= seq[i];
    return iter_item;
}
#ifdef MAIN_LKT
int main(int argc, char *argv[])
{
    if(argc < 3){
        fprintf(stderr, "[main]: lkt argc <3!\n");
        exit(1);
    }
    lookupTable_t *lkt;
    lkt = LKT_init(12);    
    LKT_build_lookuptable(argv[1], lkt);
    fprintf(stderr, "build lookuptable finish!\n");
    LKT_output(argv[2], lkt);

    fprintf(stderr, "output lookup table!\n");

    LKT_destroy(lkt);
}
#endif 
/*
uint32_t LKT_forItem2backItem(uint32_t seq, int n)
{
    uint8_t i;
    uint8_t a, b;
    uint32_t mask0, mask1, mask3;
    for(i = 0; i < n/2; ++i)
    {
        shift0 = i <<1;
        shift1 = (n-i-1)<<1;
        mask0 = 3<<shift0;
        mask1 = 3<<shift1;
        mask3 = ~(mask0|mask1);
        a = (seq & mask0) >>shift0;
        b = (seq & mask1) >>shift1;
        seq &= mask3;
        seq |= a<<shift1;
        seq |= b<<shift0; 
    }
    return seq;
}
void build_lookuptable_backward(const char* fn_lkt_forward, const char *fn_lkt_backward)
{

    lookupTable_t *lkt_for, *lkt_back;
    lkt_for = LKT_restore(fn_lkt_forward);
    lkt_back = LKT_init(lkt_for->maxLookupLen); 
    uint32_t i, j;
    for(i = 0; i < lkt_for->n_item; ++i)
    {
        j = LKT_forItem2backItem(i, lkt_for->maxLookupLen);
        lkt_back->item[j] = lkt_for->item[i];
    } 
    
    FILE * fp_lkt_back;
    fp_lkt_back = fopen(fn_lkt_backward, "wb");
    if(fp_lkt_back == NULL){
        fprintf(stderr, "[build_lookuptable]:file %s open fail!\n", fn_lkt_backward);
        exit(1);
    }
    fwrite(lkt_back, 1, sizeof(lookupTable_t)+lkt_back->n_item * sizeof(uint32_t), fp_lkt_back);
    LKT_destroy(lkt);
    fclose(fp_lkt_back);

    LKT_destroy(lkt_for);
    LKT_destroy(lkt_back);
}
*/
