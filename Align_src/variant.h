/*
 * =====================================================================================
 *
 *       Filename:  variant.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/14/2012 05:23:38 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#ifndef _HAPMAP_H_
#define _HAPMAP_H_

#include <stdint.h>
#include <stdio.h>

typedef struct {
//  char *file_ann;
    char *rname;
    uint32_t snp_num;
    uint32_t *snp_pos;
    uint8_t *snp_type;
//snp    ref        snp_type ref_tag(4bit) snp_tag[TGCA]
//A\G     T      ->                 0011      0101
    FILE *fp_variant;
} variant_t;

variant_t *variant_init(FILE *fp);
void variant_destroy(variant_t *v);
int variant_get_snpnum(variant_t *v);
int variant_readv(variant_t *v);

#endif
