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
#include<stdint.h>
#include<stdio.h>
typedef struct {
//  char *file_ann;
    char *chrID;
    uint32_t snp_num;
    uint32_t *snp_pos;
    uint8_t *snp_type;
//snp    ref        snp_type ref_tag(4bit) snp_tag[TGCA]
//A\G     G      ->                 0010      0101
    FILE *fp_hapmap;
} hapmap_t;


static int occ1[16] = {  0, 1, 1, 2,
                    1, 2, 2, 3,
                    1, 2, 2, 3,
                    2, 3, 3, 4 };
static uint8_t snptype_map0[16] ={  4, 0, 1, 0, 
                                    2, 0, 1, 0, 
                                    3, 0, 1, 0, 
                                    2, 0, 1, 0 };
static uint8_t snptype_map1[16] ={  4, 4, 4, 1, 
                                    4, 2, 2, 1, 
                                    4, 3, 3, 1, 
                                    3, 2, 2, 1 };
static uint8_t snptype_map2[16] ={  4, 4, 4, 4,
                                    4, 4, 4, 2,
                                    4, 4, 4, 3,
                                    4, 3, 3, 2 };
static uint8_t snptype_map3[16] ={  4, 4, 4, 4,
                                    4, 4, 4, 4,
                                    4, 4, 4, 4,
                                    4, 4, 4, 3 };

hapmap_t *hapmap_init(FILE *fp);
void hapmap_destroy(hapmap_t *hm);
int hapmap_get_snpnum(hapmap_t *hm);
int hapmap_readhm(hapmap_t *hm);
static inline int hapmap_get_snptypenum(uint8_t snp_type)
{
   return occ1[snp_type&15]; 
}
static inline uint8_t hapmap_get_snptype(const uint8_t snptype, const uint8_t snp_type_iter)
{
    uint8_t ret_snp_type;
    switch (snp_type_iter){
      case 0: ret_snp_type = snptype_map0[snptype&15]; break;
      case 1: ret_snp_type = snptype_map1[snptype&15]; break;
      case 2: ret_snp_type = snptype_map2[snptype&15]; break;
      case 3: ret_snp_type = snptype_map3[snptype&15]; break;
      default :fprintf(stderr, "snp_type_iter > 3!\n");exit(1);    
    }
    return ret_snp_type;
}
uint8_t hapmap_nt2snptypei(uint8_t snptype, const uint8_t ref_nt);
#endif
