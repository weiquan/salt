/*
 * =====================================================================================
 *
 *       Filename:  bwt.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/15/2012 10:40:16 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <assert.h>
#include "rbwt.h"
#include "kseq.h"
#include "hapmap.h"

#ifdef DEBUG 
    #include <assert.h>
#endif


KSEQ_INIT(gzFile, gzread)

static void initializeVAL(unsigned int *startAddr, const unsigned int length, const unsigned int initValue)
{
    unsigned int i;
    for (i=0; i<length; ++i) startAddr[i] = initValue;
}
//copy from bwt_gen.c
static inline unsigned int BWTOccValueExplicit(const rbwt_t *bwt, const unsigned int occExplicitIndex,
                                               const unsigned int character)
{
    unsigned int occIndexMajor;

    occIndexMajor = occExplicitIndex * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

    if (occExplicitIndex % OCC_VALUE_PER_WORD == 0) {
#ifdef DEBUG
        if(occIndexMajor*ALPHABET_SIZE+character >= bwt->occMajorSizeInWord){
            fprintf(stderr, "index%u occMajorSizeInWord%u\n", occIndexMajor*ALPHABET_SIZE+character, bwt->occMajorSizeInWord); 
        }
        if(occExplicitIndex/OCC_VALUE_PER_WORD*ALPHABET_SIZE+character >= bwt->occSizeInWord){
            fprintf(stderr, "index%u occSizeInWord%u\n", occExplicitIndex/OCC_VALUE_PER_WORD*ALPHABET_SIZE+character, bwt->occSizeInWord); 
        }
        if(occExplicitIndex/OCC_VALUE_PER_WORD*ALPHABET_SIZE+character ==22){
            fprintf(stderr, "occExplicitIndex%u character%d\n", occExplicitIndex, character); 
        }
#endif
        return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
               (bwt->occValue[occExplicitIndex / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] >> 16);

    } else {
#ifdef DEBUG
        if(occIndexMajor*ALPHABET_SIZE+character >= bwt->occMajorSizeInWord){
            fprintf(stderr, "index%u occMajorSizeInWord%u\n", occIndexMajor*ALPHABET_SIZE+character, bwt->occMajorSizeInWord); 
        }
        if(occExplicitIndex/OCC_VALUE_PER_WORD*ALPHABET_SIZE+character >= bwt->occSizeInWord){
            fprintf(stderr, "index%u occSizeInWord%u\n", occExplicitIndex/OCC_VALUE_PER_WORD*ALPHABET_SIZE+character, bwt->occSizeInWord); 
        }
        if(occExplicitIndex/OCC_VALUE_PER_WORD*ALPHABET_SIZE+character ==22){
            fprintf(stderr, "occExplicitIndex%u character%d\n", occExplicitIndex, character); 
        }
#endif



       return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
               (bwt->occValue[occExplicitIndex / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] & 0x0000FFFF);
    }
}
static DECODE_TABLE_T __ForwardDNAOccCount(const unsigned int*  dna, const unsigned int index,
                                       const DECODE_TABLE_T*  dnaDecodeTable)
{   
    static const unsigned int truncateRightMask[8] = { 0x00000000, 0xF0000000,0xFF000000, 0xFFF00000, 
                                                       0xFFFF0000, 0xFFFFF000,0xFFFFFF00, 0xFFFFFFF0};
    unsigned int wordToCount, charToCount;
    unsigned int i, c;
    DECODE_TABLE_T sum;
    
    wordToCount = index / CHAR_PER_WORD;
    charToCount = index - wordToCount * CHAR_PER_WORD;
    
    sum = 0;    
    for (i=0; i<wordToCount; ++i) {
        sum += dnaDecodeTable[dna[i] >> 16];
        sum += dnaDecodeTable[dna[i] & 0x0000FFFF];
    }

    if (charToCount > 0) {
        c = dna[i] & truncateRightMask[charToCount];    // increase count of 'a' by 16 - c;
        sum += dnaDecodeTable[c >> 16];
        sum += dnaDecodeTable[c & 0xFFFF];
        sum -= (DECODE_TABLE_T)(CHAR_PER_WORD - charToCount); // decrease count of 'a' by 16 - positionToProcess
    }
    
    return sum;
}
static inline unsigned int ForwardDNAOccCount(const unsigned int *dna, const unsigned int index, const unsigned int character,
                                              const DECODE_TABLE_T* dnaDecodeTable)
{
    //fprintf(stderr, "for!!");
    DECODE_TABLE_T sum;
    sum = __ForwardDNAOccCount(dna, index, dnaDecodeTable);
    sum >>= character * COUNT_BITS_PER_CHAR;
    sum &= 0xFF;
    return sum;
}
//copy from bwt_gen.c
static DECODE_TABLE_T __BackwardDNAOccCount(const unsigned int*  dna, const unsigned int index, 
                                        const DECODE_TABLE_T*  dnaDecodeTable)
{
    static const unsigned int truncateLeftMask[8] = {   0x00000000,0x0000000F, 0x000000FF, 0x00000FFF,
                                                        0x0000FFFF,0x000FFFFF, 0x00FFFFFF, 0x0FFFFFFF};
    unsigned int wordToCount, charToCount;
    unsigned int i, c;
    DECODE_TABLE_T sum;
    
    wordToCount = index / CHAR_PER_WORD;
    charToCount = index - wordToCount * CHAR_PER_WORD;

    dna -= wordToCount + 1;
    sum = 0;
    if (charToCount > 0) {
        c = *dna & truncateLeftMask[charToCount];   // increase count of 'a' by 16 - c;
        sum += dnaDecodeTable[c >> 16];
        sum += dnaDecodeTable[c & 0xFFFF];
        sum -= (DECODE_TABLE_T)( CHAR_PER_WORD - charToCount); // decrease count of 'a' by 16 - positionToProcess
    }

    for (i=0; i<wordToCount; ++i) {
        ++dna;
        sum += dnaDecodeTable[*dna >> 16];
        sum += dnaDecodeTable[*dna & 0x0000FFFF];
    }
    
    return sum;
}
static inline unsigned int BackwardDNAOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character,
                                        const DECODE_TABLE_T*  dnaDecodeTable)
{
    DECODE_TABLE_T sum;
    //fprintf(stderr, "back!!");
    sum = __BackwardDNAOccCount(dna, index, dnaDecodeTable);
    sum >>= character * COUNT_BITS_PER_CHAR;
    sum &= 0xFF;
    
    return sum;
}
unsigned int Rbwt_BWTOccValue(const rbwt_t *bwt, unsigned int index, const unsigned int character) 
{
    unsigned int occValue;
    unsigned int occExplicitIndex, occIndex;
    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    if (index > bwt->inverseSa0) {
        index--;
    }
    
    occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;   // Bidirectional encoding
    occIndex = occExplicitIndex * OCC_INTERVAL;
    occValue = BWTOccValueExplicit(bwt, occExplicitIndex, character);
    
    if(index > bwt->textLength){
        fprintf(stderr, "occIndex > textLength!\n");
        exit(1);
    }
/*
    if(index == bwt->textLength){
        return bwt->cumulativeFreq[character+1] - bwt->cumulativeFreq[character]; 
    }
*/
    if (occIndex == index) {
        return occValue;
    }
    if (occIndex < index) {
        return occValue + ForwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, index - occIndex, character, bwt->decodeTable);
    } else {
        //assert(occIndex/CHAR_PER_WORD*CHAR_PER_WORD + occIndex -index <= bwt->textLength);
        return occValue - BackwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, occIndex - index, character, bwt->decodeTable);
    }
}


void GenerateDNAOccCountTable(DECODE_TABLE_T *dnaDecodeTable)
{
	unsigned int i, j, c, t;
    for (i=0; i<DNA_OCC_CNT_TABLE_SIZE_IN_WORD; ++i) {
		dnaDecodeTable[i] = 0;
		c = i;
		for (j=0; j<4; ++j) {
			t = c & MASK_CHAR_IN_WORD;
			dnaDecodeTable[i] += (DECODE_TABLE_T)1 << (t * COUNT_BITS_PER_CHAR);
			c >>= BIT_PER_CHAR;
		}
	}
}
uint32_t Rbwt_BWTOccValue2(const rbwt_t *bwt, unsigned int index, const unsigned int character, unsigned int *occSharp)
{
    unsigned int occValue, occValueSharp;
    unsigned int occExplicitIndex, occIndex;
    unsigned int occIndexMajor;
    DECODE_TABLE_T DNAOccCount;
    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    if (index > bwt->inverseSa0) {
        index--;
    }

    occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;   // Bidirectional encoding
    occIndex = occExplicitIndex * OCC_INTERVAL;
    occIndexMajor = occIndex / OCC_INTERVAL_MAJOR;

    unsigned int i1, i2;
    i1 = occIndexMajor*ALPHABET_SIZE;
    i2 = occExplicitIndex/ OCC_VALUE_PER_WORD *ALPHABET_SIZE;
    if (occExplicitIndex % OCC_VALUE_PER_WORD == 0) {
        occValue = bwt->occValueMajor[i1 + character] +
               (bwt->occValue[i2 + character] >> 16);
        occValueSharp = bwt->occValueMajor[i1 + NT_SHARP] +
               (bwt->occValue[i2 + NT_SHARP] >> 16);
    } else {
        occValue = bwt->occValueMajor[i1 + character] +
               (bwt->occValue[i2 + character] & 0x0000FFFF);
        occValueSharp = bwt->occValueMajor[i1 + NT_SHARP] +
               (bwt->occValue[i2 + NT_SHARP] & 0x0000FFFF);
    }
    
    if (occIndex == index) {
        *occSharp = occValueSharp;
        return occValue;
    }    
    if(occIndex <index){
        DNAOccCount = __ForwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, index - occIndex, bwt->decodeTable);
        occValue += (DNAOccCount >> (character*COUNT_BITS_PER_CHAR)) & 0xFF; 
        occValueSharp += (DNAOccCount >> (NT_SHARP*COUNT_BITS_PER_CHAR)) & 0xFF;
    } else{
        DNAOccCount = __BackwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, occIndex - index, bwt->decodeTable);
        occValue -= (DNAOccCount >> (character*COUNT_BITS_PER_CHAR)) & 0xFF; 
        occValueSharp -= (DNAOccCount >> (NT_SHARP *COUNT_BITS_PER_CHAR))& 0xFF;
    }
    *occSharp = occValueSharp;
    return occValue;
}

//将获得bwt pos位置的字符。

//bwtio
void Rbwt_restore_bwt(rbwt_t *bwt, const char *fn_bwt)
{
    FILE *bwt_file = NULL;
//    uint32_t bwtCodeSizeInWord;   
    uint32_t bwtCodeSizeInWord; 
    bwt_file = (FILE*)fopen(fn_bwt, "rb");
	if (bwt_file == NULL) {
		fprintf(stderr, "restore_bwt: Cannot open BWT code file!\n");
		exit(1);
	}
    
    fread(&bwt->textLength, sizeof(uint32_t),1, bwt_file);
   	fread(&bwt->inverseSa0, sizeof(uint32_t), 1, bwt_file);
	bwt->cumulativeFreq[0] = 0;
    fread(bwt->cumulativeFreq + 1, sizeof(uint32_t), ALPHABET_SIZE, bwt_file);
    fread(&bwt->bwtSizeInWord, sizeof(uint32_t), 1, bwt_file);
    bwtCodeSizeInWord = (bwt->bwtSizeInWord *CHAR_PER_WORD +OCC_INTERVAL) / OCC_INTERVAL * OCC_INTERVAL/ CHAR_PER_WORD +1;
    bwt->bwtCode = (uint32_t *)calloc(bwtCodeSizeInWord, sizeof(uint32_t));
    fread(bwt->bwtCode, sizeof(uint32_t), bwt->bwtSizeInWord, bwt_file);

    fclose(bwt_file);
}
void Rbwt_restore_occ( rbwt_t *bwt, const char *fn_occ)
{
    FILE *occ_file = NULL;
    
    occ_file = (FILE*)fopen(fn_occ, "rb");
	if (occ_file == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open occ value file!\n");
		exit(1);
	}
	
    fread(&bwt->occSizeInWord, sizeof(uint32_t), 1, occ_file);
    bwt->occValue = (uint32_t *)calloc(bwt->occSizeInWord, sizeof(uint32_t));
    fread(bwt->occValue, sizeof(uint32_t), bwt->occSizeInWord, occ_file);
	fread(&bwt->occMajorSizeInWord, sizeof(uint32_t), 1, occ_file);
    bwt->occValueMajor = (uint32_t *)calloc(bwt->occMajorSizeInWord, sizeof(uint32_t));	
    fread(bwt->occValueMajor, sizeof(uint32_t), bwt->occMajorSizeInWord, occ_file);
	fclose(occ_file);
}
uint32_t Rbwt_for_bwt_sa(rbwt_t *bwt,  uint32_t sa_index)
{
    uint32_t step = 0;
    uint8_t c;
    
    while(sa_index <= bwt->cumulativeFreq[NT_SHARP] ){
        c = Rbwt_bwt2nt(bwt, sa_index);
        sa_index = bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, sa_index, c) +1;
        ++step;
    }
    if(sa_index > bwt->cumulativeFreq[NT_SHARP]){
        //plus 1 for  SA store the pos of the character which before '#'
        return  bwt->saValueSharp[sa_index - bwt->cumulativeFreq[NT_SHARP] -1] - step+1;
    }  else {
        fprintf(stderr, "bwt_sa error!\n");
        exit(1);
    }
}
uint32_t Rbwt_back_bwt_sa(rbwt_t *bwt,  uint32_t sa_index)
{
    uint32_t step = 0;
    uint8_t c;
    
    while(sa_index <= bwt->cumulativeFreq[NT_SHARP] ){
        c = Rbwt_bwt2nt(bwt, sa_index);
        sa_index = bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, sa_index, c) +1;
        ++step;
    }
    if(sa_index > bwt->cumulativeFreq[NT_SHARP]){
        //minus 1 for  SA store the pos of the character which before '#'
        return  bwt->saValueSharp[sa_index - bwt->cumulativeFreq[NT_SHARP] -1] + step -1;
    }  else {
        fprintf(stderr, "bwt_sa error!\n");
        exit(1);
    }
}
sharp2Ri_t *Rbwt_sharp2Ri_init(const char *fn_localpattern)
{
#define L_LP 50    
    int i;

    sharp2Ri_t *sharp2Ri = NULL;
    sharp2Ri = calloc(1, sizeof(sharp2Ri_t));
    if(sharp2Ri == NULL) {
        fprintf(stderr, "[Rbwt_sharp2Ri_init]: mem allocate fail!\n");
        exit(EXIT_FAILURE); 
    }
    
    /*****************************************************
     * compute number of local pattern
     */
    kseq_t *seq;
    int n_lp = 0;
    gzFile fp_LP = Z_NULL;
    fp_LP = gzopen(fn_localpattern, "r");
    seq = kseq_init(fp_LP);
    //compute LP number in first iter
    kseq_read(seq);
    for(i = 0; i < seq->seq.l; ++i){
        if(seq->seq.s[i] == '#'){
           ++n_lp; 
        } 
    }
    //compute LP number loop
    while(kseq_read(seq) > 0){
        for(i = 0; i < seq->seq.l; ++i){
            if(seq->seq.s[i] == '#'){
                ++n_lp; 
            } 
        }
    }
    
    kseq_destroy(seq);   
    gzclose(fp_LP);

    /*******************************************************
     * set sharp2Ri array
     */
 
    int j;
    unsigned int pos;
    int *sharp2Ri_array = (int *)malloc(sizeof(int) * n_lp);
    
    fp_LP = gzopen(fn_localpattern, "r");
    seq = kseq_init(fp_LP);
    
    j = 0;
    kseq_read(seq);
    pos =(unsigned int) atol(seq->comment.s);
    for(i = 0; i < seq->seq.l; ++i){
        if(seq->seq.s[i] == '#'){
            sharp2Ri_array[j++] = pos;
        }
    }


    while(kseq_read(seq) > 0){
        pos = (unsigned int)atol(seq->comment.s);
        for(i = 0; i < seq->seq.l; ++i){
             if(seq->seq.s[i] == '#'){
                sharp2Ri_array[j++] = pos;
             }
        
        }


    }
    
    kseq_destroy(seq);
    gzclose(fp_LP);
    
    sharp2Ri->n = n_lp;
    sharp2Ri->sharp2Ri_array = sharp2Ri_array;
    return sharp2Ri;
}




void Rbwt_sharp2Ri_destroy(sharp2Ri_t *sharp2Ri)
{
    free(sharp2Ri->sharp2Ri_array);
    free(sharp2Ri);
}
//forward direction =1;
//backward direction =-1;
void Rbwt_gen_sa(rbwt_t *bwt, const sharp2Ri_t *sharp2Ri,  const int direction)
{
#define SA_INTERVAL 32

    const uint32_t n_SHARP = bwt->cumulativeFreq[NT_SHARP];
    int i;
    uint32_t sa_index, l ;
    uint8_t c;
     
    bwt->saValueSizeSharp = bwt->textLength - bwt->cumulativeFreq[NT_SHARP] +1;
    bwt->saValueSharp     = (uint32_t *)calloc(bwt->saValueSizeSharp, sizeof(uint32_t));   
    //gen SA core loop     
    l = 0;sa_index = 0; i = direction >0?0:sharp2Ri->n;
    if(direction <0){
        c = Rbwt_bwt2nt(bwt, sa_index);
//        putchar("ACGT#"[c]);
//        if(c == 4){
//            printf("0\n");
//        }
        sa_index = bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, sa_index, c) +1;
        i += direction;
    }
    int l_alt_seq = 0;
    while( l < bwt->textLength ){
        c = Rbwt_bwt2nt(bwt, sa_index);

//        putchar("ACGT#56$"[c]);
        if(c ==7 ){
            sa_index = 0;
        } else{ 
            sa_index = bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, sa_index, c) +1;
            ++l_alt_seq;
        }
        if(sa_index > n_SHARP){
            if(direction >0 ){
                bwt->saValueSharp[sa_index - n_SHARP -1] = sharp2Ri->sharp2Ri_array[i];
            } else{
                bwt->saValueSharp[sa_index - n_SHARP -1] = sharp2Ri->sharp2Ri_array[i+1] -l_alt_seq;
            }
            i += direction;
            l_alt_seq = 0;
        } else{
        
        }
/*
        if(c == 4){
            printf("%d\n", bwt->saValueSharp[sa_index - n_SHARP -1]);
        }
*/
        ++l;
    }
}
//haven't tested
rbwt2_t *Rbwt2_init(const char *prefix)
{
    rbwt2_t *rbwt2 = NULL;
    char bwt0[256], bwt1[256], occ0[256], occ1[256], sa0[256], sa1[256]; 
    
    strncpy(bwt0, prefix, 256);strcat(bwt0, ".bwt_0");
    strncpy(bwt1, prefix, 256);strcat(bwt1, ".bwt_1");
    strncpy(occ0, prefix, 256);strcat(occ0, ".occ_0");
    strncpy(occ1, prefix, 256);strcat(occ1, ".occ_1");
    strncpy(sa0, prefix, 256);strcat(sa0, ".sa_0");
    strncpy(sa1, prefix, 256);strcat(sa1, ".sa_1");
      
    rbwt2 = calloc(1, sizeof(rbwt2_t));
    if(rbwt2 == NULL) {
        fprintf(stderr, "[Rbwt2_init]:rbwt2 allocate memory fail!\n");
        exit(1);
    } 
    //restore bwt ,occ and gen SA
      
    rbwt2->rbwt0 = Rbwt_init(bwt0, occ0);
    Rbwt_restore_sa(rbwt2->rbwt0, sa0);
    rbwt2->rbwt1 = Rbwt_init(bwt1, occ1);
    Rbwt_restore_sa(rbwt2->rbwt1, sa1);
    // Generate decode tables
    rbwt2->rbwt1->decodeTable = rbwt2->rbwt0->decodeTable = (DECODE_TABLE_T *)calloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD, sizeof(DECODE_TABLE_T));
	GenerateDNAOccCountTable(rbwt2->rbwt0->decodeTable);
    
    return rbwt2;
}
void Rbwt2_save_sa(const char *prefix,const char *fn_localpattern)
{
    rbwt2_t *rbwt2 = NULL;
    char bwt0[256], bwt1[256], occ0[256], occ1[256]; 
    
    strncpy(bwt0, prefix, 256);strcat(bwt0, ".forward.bwt");
    strncpy(bwt1, prefix, 256);strcat(bwt1, ".backward.bwt");
    strncpy(occ0, prefix, 256);strcat(occ0, ".forward.occ");
    strncpy(occ1, prefix, 256);strcat(occ1, ".backward.occ");
   
    rbwt2 = calloc(1, sizeof(rbwt2_t));
    if(rbwt2 == NULL) {
        fprintf(stderr, "[Rbwt2_init]:rbwt2 allocate memory fail!\n");
        exit(1);
    } 
    //restore bwt ,occ and gen SA
      
    rbwt2->rbwt0 = Rbwt_init(bwt0, occ0);
    rbwt2->rbwt1 = Rbwt_init(bwt1, occ1);
    
    rbwt_t *rbwt0, *rbwt1;
    rbwt0 = rbwt2->rbwt0;
    rbwt1 = rbwt2->rbwt1; 
    // Generate decode tables
    rbwt1->decodeTable = rbwt0->decodeTable = (DECODE_TABLE_T *)calloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD, sizeof(DECODE_TABLE_T));
	GenerateDNAOccCountTable(rbwt0->decodeTable);
    
    sharp2Ri_t *sharp2Ri;
    sharp2Ri = Rbwt_sharp2Ri_init(fn_localpattern);
    Rbwt_gen_sa(rbwt2->rbwt0, sharp2Ri, 1);
    Rbwt_gen_sa(rbwt2->rbwt1, sharp2Ri, -1);
    Rbwt_sharp2Ri_destroy(sharp2Ri);
    
    char sa0[256], sa1[256];
    FILE *fp_sa0, *fp_sa1;
    
    strncpy(sa0, prefix, 256);strcat(sa0, ".forward.sa");
    strncpy(sa1, prefix, 256);strcat(sa1, ".backward.sa");
    fp_sa0 = fopen(sa0, "wb");
    fwrite(&rbwt0->saValueSizeSharp, 1, sizeof(uint32_t), fp_sa0);
    fwrite(rbwt0->saValueSharp, rbwt0->saValueSizeSharp, sizeof(uint32_t), fp_sa0);
    fp_sa1 = fopen(sa1, "wb");
    fwrite(&rbwt1->saValueSizeSharp, 1, sizeof(uint32_t), fp_sa1);
    fwrite(rbwt1->saValueSharp, rbwt0->saValueSizeSharp, sizeof(uint32_t), fp_sa1);
    fclose(fp_sa0);
    fclose(fp_sa1);
    
    Rbwt2_destroy(rbwt2);
}
void Rbwt_restore_sa(rbwt_t *rbwt, const char *fn_sa)
{
    FILE *fp_sa;
    fp_sa = fopen(fn_sa, "rb");
    if(fp_sa == NULL){
        fprintf(stderr, "[Rbwt_restore_sa]:file %s open fail!\n", fn_sa);
        exit(EXIT_FAILURE); 
    } 
    fread(&rbwt->saValueSizeSharp, 1, sizeof(uint32_t), fp_sa);
    rbwt->saValueSharp = calloc(rbwt->saValueSizeSharp, sizeof(uint32_t));
    if(rbwt->saValueSharp == NULL){
        fprintf(stderr, "[Rbwt_restore_sa]:allocate mem fail!\n");
        exit(EXIT_FAILURE); 
    }
    fread(rbwt->saValueSharp, rbwt->saValueSizeSharp, sizeof(uint32_t), fp_sa);
    fclose(fp_sa);
}
void Rbwt2_destroy(rbwt2_t *rbwt2)
{
    free(rbwt2->rbwt0->decodeTable);
    Rbwt_destroy(rbwt2->rbwt0);
    Rbwt_destroy(rbwt2->rbwt1);
    free(rbwt2);
}
rbwt_t *Rbwt_init(const char *fn_bwt, const char* fn_occ)
{
    rbwt_t *bwt;
	bwt = (rbwt_t*)calloc(1, sizeof(rbwt_t));

    bwt->cumulativeFreq = (uint32_t*)calloc((ALPHABET_SIZE + 1), sizeof(uint32_t));
	initializeVAL(bwt->cumulativeFreq, ALPHABET_SIZE + 1, 0);

	bwt->bwtSizeInWord = 0;
	//bwt->saValueOnBoundary = NULL;
    
    Rbwt_restore_bwt(bwt, fn_bwt);  
    Rbwt_restore_occ(bwt, fn_occ);
    //no need
    //bwt->inverseSaInterval = ALL_ONE_MASK;
    //bwt->inverseSaSize = 0;
	//bwt->inverseSa = NULL;
    return bwt;
}
//haven't tested
static inline int Rbwt_bw_search(rbwt_t *bwt, const SaIndexRange *in_saIndex, uint8_t c, SaIndexRange *out_saIndex)
{
#ifdef DEBUG
    if(in_saIndex->endSaIndex > bwt->textLength || in_saIndex->startSaIndex > bwt->textLength){
        fprintf(stderr, "bwt->textLength = %u\n", bwt->textLength);
        fprintf(stderr, "[sp, ep] = [%u, %u]\n", in_saIndex->startSaIndex, in_saIndex->endSaIndex);
        fprintf(stderr, "LF mapping out of range !!!\n");
        exit(1);
    }
#endif
    out_saIndex->startSaIndex = bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, in_saIndex->startSaIndex, c) +1;
    out_saIndex->endSaIndex =  bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, in_saIndex->endSaIndex +1, c);
    return out_saIndex->endSaIndex - out_saIndex->startSaIndex;
}

int Rbwt_exact_match_backward(const rbwt_t *bwt, const uint8_t *query, const int qlen, uint32_t *k, uint32_t *l)
{
    //unsigned int *in_sp, *in_ep, *out_sp, *out_ep; 
    uint32_t k0, l0;
#ifdef DEBUG
    assert(bwt && query && k && l);
#endif
    if(qlen <=0){
        return 0;
    }
   
    k0 = *k;
    l0 = *l;
    
    int step = 0;
    while(k0 <= l0 && step < qlen)
    {
#ifdef DEBUG        
        fprintf(stderr, "[iter %d][sp, ep] = [%u, %u]\n", step, k0, l0);
#endif
        //fprintf(stderr, "[sp, ep] = [%u, %u]\n", out_saIndex->startSaIndex, out_saIndex->endSaIndex);
        uint8_t c = query[qlen-step -1];
        k0 = bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, k0, c) +1;
        l0 =  bwt->cumulativeFreq[c] + Rbwt_BWTOccValue(bwt, l0 +1, c);
        ++step;
    }
    *k = k0;
    *l = l0;
    return l0 >= k0;
}
int Rbwt_exact_match_forward(const rbwt_t *rbwt, const uint8_t *query, const int qlen, uint32_t *k, uint32_t *l)
{
    //unsigned int *in_sp, *in_ep, *out_sp, *out_ep; 
    uint32_t k0, l0;
#ifdef DEBUG
    assert(rbwt && query && k && l);
#endif
    if(qlen <=0){
        return 0;
    }

    k0 = *k;
    l0 = *l;
    
    int step = 0;
    while(k0 <= l0 && step < qlen)
    {
#ifdef DEBUG        
        fprintf(stderr, "[iter %d][sp, ep] = [%u, %u]\n", step, k0, l0);
#endif
        uint8_t c = query[step];
        k0 = rbwt->cumulativeFreq[c] + Rbwt_BWTOccValue(rbwt, k0, c) +1;
        l0 =  rbwt->cumulativeFreq[c] + Rbwt_BWTOccValue(rbwt, l0 +1, c);
        ++step;
    }
    *k = k0;
    *l = l0;
    return l0 >= k0;
}

void Rbwt_destroy(rbwt_t* bwt)
{
    free(bwt->cumulativeFreq);
    free(bwt->bwtCode);
    free(bwt->occValue);
    free(bwt->occValueMajor);
    free(bwt);
}


