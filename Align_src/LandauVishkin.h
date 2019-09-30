/*
Copyright 2012, Regents of the University of California.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

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
/*
 *1. text string use 4bit to encode SNP
 *
 */
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

