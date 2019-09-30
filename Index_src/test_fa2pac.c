/*
 * =====================================================================================
 *
 *       Filename:  test_fa2pac.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/14/2012 07:47:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <stdint.h>
#include "bntseq.h"

#define _get_pac(pac, l) ((pac)[(l)>>1]>>((~(l)&1)<<2)&15)
int main(int argc, char *argv[])
{
    bntseq_t *bns;
    bns = bns_restore(argv[1]);
    uint8_t *pac;
    pac = calloc(bns->l_pac/2+2, 1);
    fread(pac, 1, bns->l_pac/2+2, bns->fp_pac);
    int i;
    for(i = 0; i < bns->l_pac; ++i){
        putchar( "ACGT#"[_get_pac(pac, i)]);         
    }
    bns_destroy(bns);


}

