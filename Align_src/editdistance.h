/*
 * =====================================================================================
 *
 *       Filename:  editdistance.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/29/2012 03:06:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include "LandauVishkin.h"

int ed_mismatch(const uint32_t *mixRef, uint32_t ref_st, const uint8_t *seq, uint32_t l_comp, int max_err);
int ed_diff(const uint32_t *mixRef, uint32_t l_mref, uint32_t ref_st, const uint32_t l_ref, const uint8_t *seq, uint32_t l_seq, int max_k_diff);
int ed_diff_withcigar(const uint32_t *mixRef, uint32_t ref_st, uint32_t l_ref, const uint8_t *seq, uint32_t l_seq, int max_k_diff, char *cigarBuf, int cigarLen, int useM, CigarFormat cigarFormat);
