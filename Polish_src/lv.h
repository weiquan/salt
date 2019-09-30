/*
 * =====================================================================================
 *
 *       Filename:  LandauVishkin.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/11/2014 08:20:41 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
// Computes the edit distance between two strings and returns a CIGAR string for the edits.
#include <stdlib.h>
#define bool int
typedef enum {
    COMPACT_CIGAR_STRING = 0,
    EXPANDED_CIGAR_STRING = 1,
    COMPACT_CIGAR_BINARY = 2,
} CigarFormat;
int computeEditDistance(
        const char* text, int textLen, 
        const char* pattern, int patternLen, 
        int k);

int computeEditDistanceWithCigar(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM, 
    CigarFormat format);

