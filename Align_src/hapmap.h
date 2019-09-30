/*
 * =====================================================================================
 *
 *       Filename:  hapmap.h
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
    char *chrID;
    uint32_t snp_num;
    uint32_t *snp_pos;
    uint8_t *snp_type;
//snp    ref        snp_type ref_tag(4bit) snp_tag[TGCA]
//A\G     T      ->                 0011      0101
    FILE *fp_hapmap;
} hapmap_t;

hapmap_t *hapmap_init(FILE *fp);
void hapmap_destroy(hapmap_t *hm);
int hapmap_get_snpnum(hapmap_t *hm);
int hapmap_readhm(hapmap_t *hm);

#endif
