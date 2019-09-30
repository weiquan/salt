#ifndef RBWT_H
#define RBWT_H
/*
 * =====================================================================================
 *
 *       Filename:  bwt.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/19/2012 05:10:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <stdint.h>
#define ALPHABET_SIZE               5
#define BIT_PER_CHAR                4
#define CHAR_PER_WORD               8
#define CHAR_PER_BYTE               2

#define BITS_IN_WORD 32
#define BITS_IN_BYTE 8
#define BYTES_IN_WORD 4

#define ALL_ONE_MASK 0xFFFFFFFF
#define DNA_OCC_CNT_TABLE_SIZE_IN_WORD  65536

#define BITS_PER_OCC_VALUE          16
#define OCC_VALUE_PER_WORD          2
#define OCC_INTERVAL                256
#define OCC_INTERVAL_MAJOR          65536



#define NT_A 0
#define NT_C 1
#define NT_G 2
#define NT_T 3
#define NT_SHARP 4
#define NT_DOLLOR 7

#define MASK_CHAR 0xF

#define MATCH 1
#define UNMATCH 0

#define COUNT_BITS_PER_CHAR 8
#define MASK_CHAR_IN_WORD 0x000F


#define DECODE_TABLE_T uint64_t
typedef struct SaIndexRange {
	unsigned int startSaIndex;
	unsigned int endSaIndex;
} SaIndexRange;
typedef struct rbwt_t {
	unsigned int textLength;			// length of the text
	unsigned int saInterval;			// interval between two SA values stored explicitly
	unsigned int inverseSa0;			// SA-1[0]
	unsigned int *cumulativeFreq;		// cumulative frequency
	unsigned int *bwtCode;				// bwt_t code
	unsigned int *occValue;				// Occurrence values stored explicitly
	unsigned int *occValueMajor;		// Occurrence values stored explicitly
	//unsigned int *saValue;    			// SA values stored explicitly
    unsigned int *saValueSharp;
	
	DECODE_TABLE_T  *decodeTable;			// For decoding bwt_t by table lookup
    unsigned int decodeTableGenerated;	// == TRUE if decode table is generated on load and will be freed
	unsigned int bwtSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occMajorSizeInWord;	// Temporary variable to hold the memory allocated
	//unsigned int saValueSize;			// Temporary variable to hold the memory allocated
	unsigned int saValueSizeSharp;
} rbwt_t;
typedef struct rbwt2_t{
    rbwt_t *rbwt0, *rbwt1;//forward rbwt && backward rbwt
} rbwt2_t;
typedef struct sharp2Ri_t{
    int n;
    int *sharp2Ri_array;
} sharp2Ri_t;


void Rbwt_bwt_bwtgen(const char *fn_pac, const char *fn);
void Rbwt_gen_sa(rbwt_t *bwt, const sharp2Ri_t *sharp2Ri,  int direction);

void Rbwt2_save_sa(const char *prefix,const char *fn_localpattern);
void Rbwt_restore_sa(rbwt_t *rbwt, const char *fn_sa);
void Rbwt_restore_bwt(rbwt_t *bwt, const char *bwt_file_name);
void Rbwt_restore_occ(rbwt_t *bwt, const char *occ_file_name);

int Rbwt_exact_match_forward(const rbwt_t *rbwt, const uint8_t *query, const int qlen, uint32_t *k, uint32_t *l);
int Rbwt_exact_match_backward(const rbwt_t *bwt, const uint8_t *query, const int qlen, uint32_t *k, uint32_t *l);
unsigned int Rbwt_BWTOccValue(const rbwt_t *bwt, unsigned int index, const unsigned int character); 
uint32_t Rbwt_BWTOccValue2(const rbwt_t *bwt, unsigned int index, const unsigned int character, unsigned int *occSharp);
uint32_t Rbwt_for_bwt_sa(rbwt_t *bwt,  uint32_t sa_index);
uint32_t Rbwt_back_bwt_sa(rbwt_t *bwt,  uint32_t sa_index);

rbwt_t *Rbwt_init(const char *bwt_file_name, const char* occ_file_name);
void Rbwt_destroy(rbwt_t* bwt);
rbwt2_t *Rbwt2_init(const char *prefix);
void Rbwt2_destroy(rbwt2_t *rbwt2);

static inline uint8_t Rbwt_bwt2nt(rbwt_t* bwt, uint32_t pos)
{
    uint32_t pos_in_word, pos_in_last_word;
    uint8_t nt;
    uint8_t shift;

    if(pos == bwt->inverseSa0){
        nt = NT_SHARP;  
    } else {
        if(pos > bwt->inverseSa0){//$ not encoded 
            --pos;
        } 
        pos_in_word = pos/CHAR_PER_WORD;       
        pos_in_last_word = pos%CHAR_PER_WORD;    
        shift = (CHAR_PER_WORD - pos_in_last_word-1)*BIT_PER_CHAR;
        assert(bwt->bwtSizeInWord >= pos_in_word);
        nt =  bwt->bwtCode[pos_in_word] >> shift & MASK_CHAR; 
    }
    return nt;
}

sharp2Ri_t *Rbwt_sharp2Ri_init(const char *fn_localpattern);
void Rbwt_sharp2Ri_destroy(sharp2Ri_t *sharp2Ri);

















#endif
