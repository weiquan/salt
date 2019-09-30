/*

   BWTConstruct.c		BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/
/*
 * This program is modified in order to construct 5-alphabet bwt.
 * modified by Quan Wei <wquanhit@gmail.com>
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "4bit_bwt_gen.h"
#include "QSufSort.h"

static unsigned int TextLengthFromBytePacked(unsigned int bytePackedLength, unsigned int bitPerChar,
											 unsigned int lastByteLength)
{
	if (bytePackedLength > ALL_ONE_MASK / (BITS_IN_BYTE / bitPerChar)) {
		fprintf(stderr, "TextLengthFromBytePacked(): text length > 2^32!\n");
		exit(1);
	}
	return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;
}

static void initializeVAL(unsigned int *startAddr, const unsigned int length, const unsigned int initValue)
{
	unsigned int i;
	for (i=0; i<length; ++i) startAddr[i] = initValue;
}

#define COUNT_BITS_PER_CHAR 8
#define MASK_CHAR_IN_WORD 0x000F
static void GenerateDNAOccCountTable(DECODE_TABLE_T *dnaDecodeTable)
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
// for BWTIncCreate()
static unsigned int BWTOccValueMajorSizeInWord(const unsigned int numChar)
{
	unsigned int numOfOccValue;
	unsigned int numOfOccIntervalPerMajor;
	numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1; // Value at both end for bi-directional encoding
	numOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
	return (numOfOccValue + numOfOccIntervalPerMajor - 1) / numOfOccIntervalPerMajor * ALPHABET_SIZE;
}
// for BWTIncCreate()
static unsigned int BWTOccValueMinorSizeInWord(const unsigned int numChar)
{
	unsigned int numOfOccValue;
	numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;		// Value at both end for bi-directional encoding
	return (numOfOccValue + OCC_VALUE_PER_WORD - 1) / OCC_VALUE_PER_WORD * ALPHABET_SIZE;
}
// for BWTIncCreate()
static unsigned int BWTResidentSizeInWord(const unsigned int numChar) {

	unsigned int numCharRoundUpToOccInterval;

	// The $ in BWT at the position of inverseSa0 is not encoded
	numCharRoundUpToOccInterval = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;

	return (numCharRoundUpToOccInterval + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

static void BWTIncSetBuildSizeAndTextAddr(BWTInc *bwtInc)
{
	unsigned int maxBuildSize;

	if (bwtInc->bwt->textLength == 0) {
		// initial build
		// Minus 2 because n+1 entries of seq and rank needed for n char
		maxBuildSize = (bwtInc->availableWord - 2 - OCC_INTERVAL / CHAR_PER_WORD)/ (2 * CHAR_PER_WORD + 1) * CHAR_PER_WORD;
		if (bwtInc->initialMaxBuildSize > 0) {
			bwtInc->buildSize = min(bwtInc->initialMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
	} else {
		// Minus 3 because n+1 entries of sorted rank, seq and rank needed for n char
		// Minus numberOfIterationDone because bwt slightly shift to left in each iteration
		maxBuildSize = (bwtInc->availableWord - bwtInc->bwt->bwtSizeInWord - bwtInc->bwt->occSizeInWord - 3
							 - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR)
							 / 3;
		if (maxBuildSize < CHAR_PER_WORD) {
			fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
			exit(1);
		}
		if (bwtInc->incMaxBuildSize > 0) {
            bwtInc->buildSize = min(bwtInc->incMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
		if (bwtInc->buildSize < CHAR_PER_WORD) {
			bwtInc->buildSize = CHAR_PER_WORD;
		}
	}

	if (bwtInc->buildSize < CHAR_PER_WORD) {
		fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
		exit(1);
	}

	bwtInc->buildSize = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;

	bwtInc->packedText = bwtInc->workingMemory + 2 * (bwtInc->buildSize + 1);
	bwtInc->textBuffer = (unsigned char*)(bwtInc->workingMemory + bwtInc->buildSize + 1);

}

// for ceilLog2()
static unsigned int leadingZero(const unsigned int input)
{
	unsigned int l;
	const static unsigned int leadingZero8bit[256] = {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
											 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
											 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
											 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	if (input & 0xFFFF0000) {
		if (input & 0xFF000000) {
			l = leadingZero8bit[input >> 24];
		} else {
			l = 8 + leadingZero8bit[input >> 16];
		}
	} else {
		if (input & 0x0000FF00) {
			l = 16 + leadingZero8bit[input >> 8];
		} else {
			l = 24 + leadingZero8bit[input];
		}
	}
	return l;

}
// for BitPerBytePackedChar()
static unsigned int ceilLog2(const unsigned int input)
{
	if (input <= 1) return 0;
	return BITS_IN_WORD - leadingZero(input - 1);

}
// for ConvertBytePackedToWordPacked()
static unsigned int BitPerBytePackedChar(const unsigned int alphabetSize)
{
	unsigned int bitPerChar;
	bitPerChar = ceilLog2(alphabetSize);
	// Return the largest number of bit that does not affect packing efficiency
	if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar)
		bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
	return bitPerChar;
}
// for ConvertBytePackedToWordPacked()
static unsigned int BitPerWordPackedChar(const unsigned int alphabetSize)
{
	unsigned int bitPerWordPackedChar;
    bitPerWordPackedChar = ceilLog2(alphabetSize);
    if (BITS_IN_WORD / (BITS_IN_WORD / BIT_PER_CHAR) > bitPerWordPackedChar)
		bitPerWordPackedChar = BITS_IN_WORD / (BITS_IN_WORD / BIT_PER_CHAR);
    return bitPerWordPackedChar;
}

static void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize,
										  const unsigned int textLength)
{
	unsigned int i, j, k;
	unsigned int c;
	unsigned int bitPerBytePackedChar;
	unsigned int bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned int bytePerIteration;
	unsigned int byteProcessed = 0;
	unsigned int wordProcessed = 0;
	unsigned int mask, shift;

	unsigned int buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
	charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = (unsigned int)input[byteProcessed] << shift;
			for (j=0; j<charPerByte; j++) {
				buffer[k] = c & mask;
				c <<= bitPerBytePackedChar;
				k++;
			}
			byteProcessed++;
		}

		c = 0;
		for (i=0; i<charPerWord; i++) {
			c |= buffer[i] >> bitPerWordPackedChar * i;
		}
		output[wordProcessed] = c;
		wordProcessed++;

	}

	k = 0;
	for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
		c = (unsigned int)input[byteProcessed] << shift;
		for (j=0; j<charPerByte; j++) {
			buffer[k] = c & mask;
			c <<= bitPerBytePackedChar;
			k++;
		}
		byteProcessed++;
	}

	c = 0;
	for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
		c |= buffer[i] >> bitPerWordPackedChar * i;
	}
	output[wordProcessed] = c;
}

static BWT *BWTCreate(const unsigned int textLength, DECODE_TABLE_T *decodeTable)
{
	BWT *bwt;

	bwt = (BWT*)calloc(1, sizeof(BWT));

	bwt->textLength = 0;
	bwt->inverseSa = 0;

	bwt->cumulativeFreq = (unsigned*)calloc((ALPHABET_SIZE + 1), sizeof(unsigned int*));
	initializeVAL(bwt->cumulativeFreq, ALPHABET_SIZE + 1, 0);

	bwt->bwtSizeInWord = 0;
	bwt->saValueOnBoundary = NULL;

	// Generate decode tables
	if (decodeTable == NULL) {
		bwt->decodeTable = (DECODE_TABLE_T *)calloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD, sizeof(DECODE_TABLE_T));
		GenerateDNAOccCountTable(bwt->decodeTable);
	} else {
		bwt->decodeTable = decodeTable;
	}

	bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord(textLength);
	bwt->occValueMajor = (unsigned*)calloc(bwt->occMajorSizeInWord, sizeof(unsigned int));

	bwt->occSizeInWord = 0;
	bwt->occValue = NULL;

	bwt->saInterval = ALL_ONE_MASK;
	bwt->saValueSize = 0;
	bwt->saValue = NULL;

	bwt->inverseSaInterval = ALL_ONE_MASK;
	bwt->inverseSaSize = 0;
	bwt->inverseSa = NULL;

	return bwt;
}
#define MIN_AVAILABLE_WORD 0x10000
static BWTInc *BWTIncCreate(const unsigned int textLength, const float targetNBit,
					 const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize)
{
	BWTInc *bwtInc;
	unsigned int i;

	if (targetNBit == 0) {
		fprintf(stderr, "BWTIncCreate() : targetNBit = 0!\n");
		exit(1);
	}

	bwtInc = (BWTInc*)calloc(1, sizeof(BWTInc));
	bwtInc->numberOfIterationDone = 0;
	bwtInc->bwt = BWTCreate(textLength, NULL);
	bwtInc->initialMaxBuildSize = initialMaxBuildSize;
	bwtInc->incMaxBuildSize = incMaxBuildSize;
	bwtInc->targetNBit = targetNBit;
	bwtInc->cumulativeCountInCurrentBuild = (unsigned*)calloc((ALPHABET_SIZE + 1), sizeof(unsigned int));
	initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	// Build frequently accessed data
	bwtInc->packedShift = (unsigned*)calloc(CHAR_PER_WORD, sizeof(unsigned int));
	for (i=0; i<CHAR_PER_WORD; i++) {
		bwtInc->packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
	}

	bwtInc->targetTextLength = textLength;
	bwtInc->availableWord = (unsigned int)((textLength + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL / BITS_IN_WORD * bwtInc->targetNBit);
    if (bwtInc->availableWord < MIN_AVAILABLE_WORD) bwtInc->availableWord = MIN_AVAILABLE_WORD;	
    if (bwtInc->availableWord < BWTResidentSizeInWord(textLength) + BWTOccValueMinorSizeInWord(textLength)) {
		fprintf(stderr, "BWTIncCreate() : targetNBit is too low!\n");
		exit(1);
	}
	bwtInc->workingMemory = (unsigned*)calloc(bwtInc->availableWord, BYTES_IN_WORD);

	return bwtInc;

}
// for BWTIncConstruct()
static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned int* __restrict rank,
									  unsigned int* __restrict cumulativeCount, const unsigned int numChar)
{
	unsigned int i, j;
	unsigned int c, t;
	unsigned int packedMask;
	unsigned int rankIndex;
	unsigned int lastWord, numCharInLastWord;

	lastWord = (numChar - 1) / CHAR_PER_WORD;
	numCharInLastWord = numChar - lastWord * CHAR_PER_WORD;

	packedMask = ALL_ONE_MASK >> (BITS_IN_WORD - BIT_PER_CHAR);
	rankIndex = numChar - 1;

	t = packedText[lastWord] >> (BITS_IN_WORD - numCharInLastWord * BIT_PER_CHAR);
	for (i=0; i<numCharInLastWord; i++) {
		c = t & packedMask;
		cumulativeCount[c+1]++;
		rank[rankIndex] = c;
		rankIndex--;
		t >>= BIT_PER_CHAR;
	}

	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
		t = packedText[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t & packedMask;
			cumulativeCount[c+1]++;
			rank[rankIndex] = c;
			rankIndex--;
			t >>= BIT_PER_CHAR;
		}
	}

	// Convert occurrence to cumulativeCount
	for (i=2; i<=ALPHABET_SIZE; i++) { 
        cumulativeCount[i] += cumulativeCount[i-1];
    }

}

#define CHAR_PER_ITERATION 256
static void ForwardDNAAllOccCountNoLimit(const unsigned int*  dna, const unsigned int index,
										 unsigned int* __restrict occCount, const DECODE_TABLE_T*  dnaDecodeTable)
{   
    static const unsigned int truncateRightMask[8] = { 0x00000000, 0xF0000000,
											   0xFF000000, 0xFFF00000, 
											   0xFFFF0000, 0xFFFFF000,
											   0xFFFFFF00, 0xFFFFFFF0};

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	DECODE_TABLE_T sum, check;
    
    for( j =0; j < ALPHABET_SIZE; j++) {occCount[j] = 0;}
    
    iteration = index / CHAR_PER_ITERATION;
	wordToCount = (index - iteration * CHAR_PER_ITERATION) / CHAR_PER_WORD;
	charToCount = index - iteration * CHAR_PER_ITERATION - wordToCount * CHAR_PER_WORD;

	for (i=0; i<iteration; i++) {
		sum = 0;		
		for(j=0; j< CHAR_PER_ITERATION/CHAR_PER_WORD; ++j){
			sum += dnaDecodeTable[*dna >> 16];
			sum += dnaDecodeTable[*dna & 0x0000FFFF];
			++dna;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
	        for (j=0;j<ALPHABET_SIZE;++j) {
				occCount[j] += sum & 0x000000FF;
				sum >>= 8;
			}
		} else {
			//check which character is 256.
	        check=0x0000000000000100;
			for (j=0;j<ALPHABET_SIZE;j++) {
				if (sum == check) {
					occCount[j] += 256;
				}
				check<<=8;
			}
        }
	}

	sum = 0;
	for (j=0; j<wordToCount; j++) {
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
		dna++;
	}

	if (charToCount > 0) {
		c = *dna & truncateRightMask[charToCount];	// increase count of 'a' by char_per_word - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum -= (DECODE_TABLE_T)charToCount - CHAR_PER_WORD;	// decrease count of 'a' by char_per_word - positionToProcess
	}
    
    for (j=0;j<ALPHABET_SIZE;++j) {
		occCount[j] += sum & 0x000000FF; sum >>= 8;
	}
}

static void BWTIncBuildPackedBwt(const unsigned int *relativeRank, unsigned int* __restrict bwt, const unsigned int numChar,
								 const unsigned int *cumulativeCount, const unsigned int *packedShift) {

	unsigned int i, c, r, k;
	unsigned int previousRank, currentRank;
	unsigned int wordIndex, charIndex;
	unsigned int inverseSa0;

	inverseSa0 = previousRank = relativeRank[0];

	for (i=1; i<=numChar; i++) {
		currentRank = relativeRank[i];
		// previousRank > cumulativeCount[c] because $ is one of the char
		
        c = 0;
        for(k = 1 ; k < ALPHABET_SIZE; ++k)
        {
            c += (previousRank>cumulativeCount[k]); 
        }
        // set bwt for currentRank
		if (c > 0) {
			// c <> 'a'
			r = currentRank;
			if (r > inverseSa0) {
				// - 1 because $ at inverseSa0 is not encoded
				r--;
			}
			wordIndex = r / CHAR_PER_WORD;
			charIndex = r - wordIndex * CHAR_PER_WORD;
			bwt[wordIndex] |= c << packedShift[charIndex];
		}
		previousRank = currentRank;
	}
}

static inline unsigned int BWTOccValueExplicit(const BWT *bwt, const unsigned int occIndexExplicit,
											   const unsigned int character)
{
	unsigned int occIndexMajor;

	occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

	if (occIndexExplicit % OCC_VALUE_PER_WORD == 0) {
		return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
			   (bwt->occValue[occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] >> 16);

	} else {
		return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
			   (bwt->occValue[occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] & 0x0000FFFF);
	}
}


static unsigned int ForwardDNAOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character,
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
		c = dna[i] & truncateRightMask[charToCount];	// increase count of 'a' by CHAR_PER_WORD - c;
	    sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0x0000FFFF];
		sum -= (DECODE_TABLE_T)(CHAR_PER_WORD - charToCount);	// decrease count of 'a' by CHAR_PER_WORD - positionToProcess
    }
	sum >>= character * COUNT_BITS_PER_CHAR;
    sum &= 0x000000FF; 
    return (unsigned int)sum;
}

static unsigned int BackwardDNAOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character,
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
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by CHAR_PER_WORD - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0x0000FFFF];
		sum -=  (DECODE_TABLE_T)(CHAR_PER_WORD -charToCount) ;	// decrease count of 'a' by CHAR_PER_WORD - positionToProcess
    }

	for (i=0; i<wordToCount; ++i) {
		++dna;
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
	}
    sum >>= character * COUNT_BITS_PER_CHAR;
    sum &= 0x000000FF; 
	return (unsigned int)sum;
}

static unsigned int BWTOccValue(const BWT *bwt, unsigned int index, const unsigned int character) 
{
	unsigned int occValue;
	unsigned int occExplicitIndex, occIndex;

	// $ is supposed to be positioned at inverseSa0 but it is not encoded
	// therefore index is subtracted by 1 for adjustment
	if (index > bwt->inverseSa0) {
		index--;
	}

	occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex = occExplicitIndex * OCC_INTERVAL;
	occValue = BWTOccValueExplicit(bwt, occExplicitIndex, character);

	if (occIndex == index) {
		return occValue;
	}
    if (occIndex < index) {
   	    return occValue + ForwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, index - occIndex, character, bwt->decodeTable);
	} else {
		return occValue - BackwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, occIndex - index, character, bwt->decodeTable);
	}

}

static unsigned int BWTIncGetAbsoluteRank(BWT *bwt, unsigned int* __restrict absoluteRank, unsigned int* __restrict seq,
										  const unsigned int *packedText, const unsigned int numChar,
										  const unsigned int* cumulativeCount, const unsigned int firstCharInLastIteration)
{
	unsigned int saIndex;
	unsigned int lastWord;
	unsigned int packedMask;
	unsigned int i, j;
	unsigned int c, t;
	unsigned int rankIndex;
	unsigned int shift;
	unsigned int seqIndexFromStart[ALPHABET_SIZE];
	unsigned int seqIndexFromEnd[ALPHABET_SIZE];

	for (i=0; i<ALPHABET_SIZE; i++) {
		seqIndexFromStart[i] = cumulativeCount[i];
		seqIndexFromEnd[i] = cumulativeCount[i+1] - 1;
	}

	shift = BITS_IN_WORD - BIT_PER_CHAR;
	packedMask = ALL_ONE_MASK >> shift;
	saIndex = bwt->inverseSa0;
	rankIndex = numChar - 1;

	lastWord = numChar / CHAR_PER_WORD;
	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
		t = packedText[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t & packedMask;
/*	
            if (saIndex == 192){
                fprintf(stderr, "cumu%u\n", bwt->cumulativeFreq[c]);
                fprintf(stderr, "BWTOCC%u\n", BWTOccValue(bwt, saIndex, c) +1);
             
            }
*/            
            saIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndex, c) + 1;
			// A counting sort using the first character of suffix is done here
			// If rank > inverseSa0 -> fill seq from end, otherwise fill seq from start -> to leave the right entry for inverseSa0
/*          
            if(saIndex > bwt->textLength){
                fprintf(stderr, "saIndex > textLength!!!!\n");
                fprintf(stderr, "%u\n", saIndex);
                fprintf(stderr, "%u\n", bwt->cumulativeFreq[c]);

                exit(1); 
            }
*/
            if (saIndex > bwt->inverseSa0) {
				seq[seqIndexFromEnd[c]] = rankIndex;
				absoluteRank[seqIndexFromEnd[c]] = saIndex;
				seqIndexFromEnd[c]--;
			} else {
				seq[seqIndexFromStart[c]] = rankIndex;
				absoluteRank[seqIndexFromStart[c]] = saIndex;
				seqIndexFromStart[c]++;
			}
			rankIndex--;
			t >>= BIT_PER_CHAR;
		}
	}

	absoluteRank[seqIndexFromStart[firstCharInLastIteration]] = bwt->inverseSa0;	// representing the substring of all preceding characters
	seq[seqIndexFromStart[firstCharInLastIteration]] = numChar;

	return seqIndexFromStart[firstCharInLastIteration];
}

static void BWTIncSortKey(unsigned int* __restrict key, unsigned int* __restrict seq, const unsigned int numItem)
{
	#define EQUAL_KEY_THRESHOLD	4	// Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

	int lowIndex, highIndex, midIndex;
	int lowPartitionIndex, highPartitionIndex;
	int lowStack[32], highStack[32];
	int stackDepth;
	int i, j;
	unsigned int tempSeq, tempKey;
	int numberOfEqualKey;

	if (numItem < 2) return;

	stackDepth = 0;

    lowIndex = 0;
    highIndex = numItem - 1;

	for (;;) {

		for (;;) {

			// Sort small array of data
			if (highIndex - lowIndex < BWTINC_INSERT_SORT_NUM_ITEM) {	 // Insertion sort on smallest arrays
				for (i=lowIndex+1; i<=highIndex; i++) {
					tempSeq = seq[i];
					tempKey = key[i];
					for (j = i; j > lowIndex && key[j-1] > tempKey; j--) {
						seq[j] = seq[j-1];
						key[j] = key[j-1];
					}
					if (j != i) {
						seq[j] = tempSeq;
						key[j] = tempKey;
					}
				}
				break;
			}
            // Choose pivot as median of the lowest, middle, and highest data; sort the three data
            midIndex = average(lowIndex, highIndex);
			if (key[lowIndex] > key[midIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[midIndex];
				key[lowIndex] = key[midIndex];
				seq[midIndex] = tempSeq;
				key[midIndex] = tempKey;
			}
			if (key[lowIndex] > key[highIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[highIndex];
				key[lowIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}
			if (key[midIndex] > key[highIndex]) {
				tempSeq = seq[midIndex];
				tempKey = key[midIndex];
				seq[midIndex] = seq[highIndex];
				key[midIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}
            // Partition data

			numberOfEqualKey = 0;

			lowPartitionIndex = lowIndex + 1;
			highPartitionIndex = highIndex - 1;

			for (;;) {
				while (lowPartitionIndex <= highPartitionIndex && key[lowPartitionIndex] <= key[midIndex]) {
					numberOfEqualKey += (key[lowPartitionIndex] == key[midIndex]);
					lowPartitionIndex++;
				}
				while (lowPartitionIndex < highPartitionIndex) {
					if (key[midIndex] >= key[highPartitionIndex]) {
						numberOfEqualKey += (key[midIndex] == key[highPartitionIndex]);
						break;
					}
					highPartitionIndex--;
				}
				if (lowPartitionIndex >= highPartitionIndex) {
					break;
				}
				tempSeq = seq[lowPartitionIndex];
				tempKey = key[lowPartitionIndex];
				seq[lowPartitionIndex] = seq[highPartitionIndex];
				key[lowPartitionIndex] = key[highPartitionIndex];
				seq[highPartitionIndex] = tempSeq;
				key[highPartitionIndex] = tempKey;
				if (highPartitionIndex == midIndex) {
					// partition key has been moved
					midIndex = lowPartitionIndex;
				}
				lowPartitionIndex++;
				highPartitionIndex--;
			}

			// Adjust the partition index
			highPartitionIndex = lowPartitionIndex;
			lowPartitionIndex--;

			// move the partition key to end of low partition
			tempSeq = seq[midIndex];
			tempKey = key[midIndex];
			seq[midIndex] = seq[lowPartitionIndex];
			key[midIndex] = key[lowPartitionIndex];
			seq[lowPartitionIndex] = tempSeq;
			key[lowPartitionIndex] = tempKey;

			if (highIndex - lowIndex + BWTINC_INSERT_SORT_NUM_ITEM <= EQUAL_KEY_THRESHOLD * numberOfEqualKey) {

				// Many keys = partition key; separate the equal key data from the lower partition
                midIndex = lowIndex;

				for (;;) {
					while (midIndex < lowPartitionIndex && key[midIndex] < key[lowPartitionIndex]) {
						midIndex++;
					}
					while (midIndex < lowPartitionIndex && key[lowPartitionIndex] == key[lowPartitionIndex - 1]) {
						lowPartitionIndex--;
					}
					if (midIndex >= lowPartitionIndex) {
						break;
					}
					tempSeq = seq[midIndex];
					tempKey = key[midIndex];
					seq[midIndex] = seq[lowPartitionIndex - 1];
					key[midIndex] = key[lowPartitionIndex - 1];
					seq[lowPartitionIndex - 1] = tempSeq;
					key[lowPartitionIndex - 1] = tempKey;
					midIndex++;
					lowPartitionIndex--;
				}

			}

			if (lowPartitionIndex - lowIndex > highIndex - highPartitionIndex) {
				// put the larger partition to stack
				lowStack[stackDepth] = lowIndex;
				highStack[stackDepth] = lowPartitionIndex - 1;
				stackDepth++;
				// sort the smaller partition first
				lowIndex = highPartitionIndex;
			} else {
				// put the larger partition to stack
				lowStack[stackDepth] = highPartitionIndex;
				highStack[stackDepth] = highIndex;
				stackDepth++;
				// sort the smaller partition first
				if (lowPartitionIndex > lowIndex) {
					highIndex = lowPartitionIndex - 1;
				} else {
					// all keys in the partition equals to the partition key
					break;
				}
			}
			continue;
		}

		// Pop a range from stack
		if (stackDepth > 0) {
			stackDepth--;
			lowIndex = lowStack[stackDepth];
			highIndex = highStack[stackDepth];
			continue;
		} else return;
	}
}

//from end to begin to check whether sortedRank[i] = sortedRank[i+1],if equall relativerank[seq[i]]=index[i+1],else sortedRank[i] = index[i].
static void BWTIncBuildRelativeRank(unsigned int* __restrict sortedRank, unsigned int* __restrict seq,
									unsigned int* __restrict relativeRank, const unsigned int numItem,
									unsigned int oldInverseSa0, const unsigned int *cumulativeCount)
{
	unsigned int i, c;
	unsigned int s, r;
	unsigned int lastRank, lastIndex;
	unsigned int oldInverseSa0RelativeRank = 0;
	unsigned int freq;

	lastIndex = numItem;
	lastRank = sortedRank[numItem];
	if (lastRank > oldInverseSa0) {
		sortedRank[numItem]--;	// to prepare for merging; $ is not encoded in bwt
	}
	s = seq[numItem];
	relativeRank[s] = numItem;
	if (lastRank == oldInverseSa0) {
		oldInverseSa0RelativeRank = numItem;
		oldInverseSa0++;	// so that this segment of code is not run again
		lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
	}

	c = ALPHABET_SIZE - 1;
	freq = cumulativeCount[c];

	for (i=numItem; i--;) {	// from numItem - 1 to 0
		r = sortedRank[i];
		if (r > oldInverseSa0) {
			sortedRank[i]--;	// to prepare for merging; $ is not encoded in bwt
		}
		s = seq[i];
		if (i < freq) {
			if (lastIndex >= freq) {
				lastRank++;	// to trigger the group across alphabet boundary to be split
			}
			c--;
			freq = cumulativeCount[c];
		}
		if (r == lastRank) {
			relativeRank[s] = lastIndex;
		} else {
			if (i == lastIndex - 1) {
				if (lastIndex < numItem && (int)seq[lastIndex + 1] < 0) {
					seq[lastIndex] = seq[lastIndex + 1] - 1;
				} else {
					seq[lastIndex] = (unsigned int)-1;
				}
			}
			lastIndex = i;
			lastRank = r;
			relativeRank[s] = i;
			if (r == oldInverseSa0) {
				oldInverseSa0RelativeRank = i;
				oldInverseSa0++;	// so that this segment of code is not run again
				lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
			}
		}
	}

}

static void BWTIncBuildBwt(unsigned int*  seq, const unsigned int *relativeRank, const unsigned int numChar,
						   const unsigned int *cumulativeCount)
{
	unsigned int i,j, c;
	unsigned int previousRank, currentRank;

	previousRank = relativeRank[0];

	for (i=1; i<=numChar; ++i) {
		currentRank = relativeRank[i];
		c = 0;
		for(j = 1; j < ALPHABET_SIZE; ++j){
			c+= (previousRank>= cumulativeCount[j]);
		}
		seq[currentRank] = c;
		previousRank = currentRank;
	}
}

static void BWTIncMergeBwt(const unsigned int *sortedRank, const unsigned int* oldBwt, const unsigned int *insertBwt,
						   unsigned int* __restrict mergedBwt, const unsigned int numOldBwt, const unsigned int numInsertBwt)
{
	unsigned int bitsInWordMinusBitPerChar;
	unsigned int leftShift, rightShift;
	unsigned int o;
	unsigned int oIndex, iIndex, mIndex;
	unsigned int mWord, mChar, oWord, oChar;
	unsigned int numInsert;

	bitsInWordMinusBitPerChar = BITS_IN_WORD - BIT_PER_CHAR;

	oIndex = 0;
	iIndex = 0;
	mIndex = 0;

	mWord = 0;
	mChar = 0;

	mergedBwt[0] = 0;	// this can be cleared as merged Bwt slightly shift to the left in each iteration

	while (oIndex < numOldBwt) {

		// copy from insertBwt
		while (iIndex <= numInsertBwt && sortedRank[iIndex] <= oIndex) {
			if (sortedRank[iIndex] != 0) {	// special value to indicate that this is for new inverseSa0
				mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
				mIndex++;
				mChar++;
				if (mChar == CHAR_PER_WORD) {
					mChar = 0;
					mWord++;
					mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
				}
			}
			iIndex++;
		}

		// Copy from oldBwt to mergedBwt
		if (iIndex <= numInsertBwt) {
			o = sortedRank[iIndex];
		} else {
			o = numOldBwt;
		}
		numInsert = o - oIndex;

		oWord = oIndex / CHAR_PER_WORD;
		oChar = oIndex - oWord * CHAR_PER_WORD;
		if (oChar > mChar) {
			leftShift = (oChar - mChar) * BIT_PER_CHAR;
			rightShift = (CHAR_PER_WORD + mChar - oChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord]
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR))
								| (oldBwt[oWord+1] >> rightShift);
			oIndex += min(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord] << leftShift) | (oldBwt[oWord+1] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else if (oChar < mChar) {
			rightShift = (mChar - oChar) * BIT_PER_CHAR;
			leftShift = (CHAR_PER_WORD + oChar - mChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord]
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR));
			oIndex += min(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord-1] << leftShift) | (oldBwt[oWord] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else { // oChar == mChar
			mergedBwt[mWord] = mergedBwt[mWord] | truncateLeft(oldBwt[oWord], mChar * BIT_PER_CHAR);
			oIndex += min(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = oldBwt[oWord];
				oIndex += CHAR_PER_WORD;
			}
		}
		oIndex = o;
		mIndex += numInsert;

		// Clear the trailing garbage in mergedBwt
		mWord = mIndex / CHAR_PER_WORD;
		mChar = mIndex - mWord * CHAR_PER_WORD;
		if (mChar == 0) {
			mergedBwt[mWord] = 0;
		} else {
			mergedBwt[mWord] = truncateRight(mergedBwt[mWord], (BITS_IN_WORD - mChar * BIT_PER_CHAR));
		}

	}

	// copy from insertBwt
	while (iIndex <= numInsertBwt) {
		if (sortedRank[iIndex] != 0) {
			mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
			mIndex++;
			mChar++;
			if (mChar == CHAR_PER_WORD) {
				mChar = 0;
				mWord++;
				mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
			}
		}
		iIndex++;
	}
}

static void BWTClearTrailingBwtCode(BWT *bwt)
{
	unsigned int bwtResidentSizeInWord;
	unsigned int wordIndex, offset;
	unsigned int i;

	bwtResidentSizeInWord = BWTResidentSizeInWord(bwt->textLength);

	wordIndex = bwt->textLength / CHAR_PER_WORD;
	offset = (bwt->textLength - wordIndex * CHAR_PER_WORD) * BIT_PER_CHAR;
	if (offset > 0) {
		bwt->bwtCode[wordIndex] = truncateRight(bwt->bwtCode[wordIndex], BITS_IN_WORD - offset);
	} else {
		if (wordIndex < bwtResidentSizeInWord) {
			bwt->bwtCode[wordIndex] = 0;
		}
	}

	for (i=wordIndex+1; i<bwtResidentSizeInWord; i++) {
		bwt->bwtCode[i] = 0;
	}
}


static void BWTGenerateOccValueFromBwt(const unsigned int*  bwt, unsigned int* __restrict occValue,
								unsigned int* __restrict occValueMajor,
								const unsigned int textLength, const DECODE_TABLE_T*  decodeTable)
{
	unsigned int numberOfOccValueMajor, numberOfOccValue;
	unsigned int wordBetweenOccValue;
	unsigned int numberOfOccIntervalPerMajor;
	unsigned int c;
	unsigned int i, j, k;
	unsigned int occMajorIndex;
	unsigned int occIndex, bwtIndex;
	DECODE_TABLE_T sum;
	unsigned int tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];

	wordBetweenOccValue = OCC_INTERVAL / CHAR_PER_WORD;

	// Calculate occValue
	// [lh3] by default: OCC_INTERVAL_MAJOR=65536, OCC_INTERVAL=256
	numberOfOccValue = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;				// Value at both end for bi-directional encoding
	numberOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
	numberOfOccValueMajor = (numberOfOccValue + numberOfOccIntervalPerMajor - 1) / numberOfOccIntervalPerMajor;

	for (i=0;i<ALPHABET_SIZE;i++) {
        tempOccValue0[i] = 0;
        occValueMajor[i] = 0;
    }

	occIndex = 0;
	bwtIndex = 0;
	for (occMajorIndex=1; occMajorIndex<numberOfOccValueMajor; occMajorIndex++) {

		for (i=0; i<numberOfOccIntervalPerMajor/2; i++) {

            for (j=0;j<ALPHABET_SIZE;j++) {
                tempOccValue1[j] = tempOccValue0[j];
            }
            
            for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum = decodeTable[c >> 16];
 				sum += decodeTable[c & 0x0000FFFF];
                for (k=0;k<ALPHABET_SIZE;k++) {
				    tempOccValue1[k] += (unsigned int)(sum & 0x00000000000000FF);	sum >>= 8;
                }
				bwtIndex++;
			}
            for (k=0;k<ALPHABET_SIZE;k++) {
			    occValue[occIndex * ALPHABET_SIZE + k] = (tempOccValue0[k] << 16) | tempOccValue1[k];
			    tempOccValue0[k] = tempOccValue1[k];
            }

            occIndex++;
            for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum = decodeTable[c >> 16];
     	        sum += decodeTable[c & 0x0000FFFF];
                for (k=0;k<ALPHABET_SIZE;k++) {
				    tempOccValue0[k] += (unsigned int)(sum & 0x00000000000000FF);	sum >>= 8;
                }
				bwtIndex++;
			}

		}
        for (k=0;k<ALPHABET_SIZE;k++) {
		    occValueMajor[occMajorIndex * ALPHABET_SIZE + k] = occValueMajor[(occMajorIndex - 1) * 
                                                               ALPHABET_SIZE + k] + tempOccValue0[k];
			//printf("DIU : occValueMajor[%u] = %u\n",occMajorIndex * ALPHABET_SIZE + k,occValueMajor[occMajorIndex * ALPHABET_SIZE + k] );
		    tempOccValue0[k] = 0;
        }
    }

	while (occIndex < (numberOfOccValue-1)/2) {

        for (k=0;k<ALPHABET_SIZE;k++) {
		    tempOccValue1[k] = tempOccValue0[k];
        }
        for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum = decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
			    tempOccValue1[k] += (unsigned int)(sum & 0x00000000000000FF);	sum >>= 8;
            }
			bwtIndex++;
		}			sum = 0;
        for (k=0;k<ALPHABET_SIZE;k++) {
		    occValue[occIndex * ALPHABET_SIZE + k] = (tempOccValue0[k] << 16) | tempOccValue1[k];
		    tempOccValue0[k] = tempOccValue1[k];
        }

		occIndex++;
        for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum = decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
			    tempOccValue0[k] += (unsigned int)(sum & 0x00000000000000FF);	sum >>= 8;
            }
			bwtIndex++;
		}
    }

    for (k=0;k<ALPHABET_SIZE;k++) {
	    tempOccValue1[k] = tempOccValue0[k];
    }
    if (occIndex * 2 < numberOfOccValue - 1) {
		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum = decodeTable[c >> 16];
            sum += decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
			    tempOccValue1[k] += (unsigned int)(sum & 0x00000000000000FF);	sum >>= 8;
            }
			bwtIndex++;
		}
    }
    for (k=0;k<ALPHABET_SIZE;k++) {
	    occValue[occIndex * ALPHABET_SIZE + k] = (tempOccValue0[k] << 16) | tempOccValue1[k];
    }


}

static void BWTIncConstruct(BWTInc *bwtInc, const unsigned int numChar)
{
	unsigned int i;
	unsigned int mergedBwtSizeInWord, mergedOccSizeInWord;
	unsigned int firstCharInThisIteration;

	unsigned int *relativeRank, *seq, *sortedRank, *insertBwt, *mergedBwt;
	unsigned int newInverseSa0RelativeRank, oldInverseSa0RelativeRank, newInverseSa0;

	#ifdef DEBUG
	if (numChar > bwtInc->buildSize) {
		fprintf(stderr, "BWTIncConstruct(): numChar > buildSize!\n");
		exit(1);
	}
	#endif

	mergedBwtSizeInWord = BWTResidentSizeInWord(bwtInc->bwt->textLength + numChar);
	mergedOccSizeInWord = BWTOccValueMinorSizeInWord(bwtInc->bwt->textLength + numChar);

	initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	if (bwtInc->bwt->textLength == 0) {		// Initial build
        // Set address
		seq = bwtInc->workingMemory;
		relativeRank = seq + bwtInc->buildSize + 1;
		mergedBwt = insertBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord;	// build in place

		BWTIncPutPackedTextToRank(bwtInc->packedText, relativeRank, bwtInc->cumulativeCountInCurrentBuild, numChar);

		firstCharInThisIteration = relativeRank[0];
		relativeRank[numChar] = 0;

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)ALPHABET_SIZE - 1, 0, FALSE);
		newInverseSa0 = relativeRank[0];

		// Clear BWT area
		initializeVAL(insertBwt, mergedBwtSizeInWord, 0);

		// Build BWT
		BWTIncBuildPackedBwt(relativeRank, insertBwt, numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->packedShift);

		// so that the cumulativeCount is not deducted
		bwtInc->firstCharInLastIteration = ALPHABET_SIZE;

	} else {		// Incremental build
		// Set address
		sortedRank = bwtInc->workingMemory;
		seq = sortedRank + bwtInc->buildSize + 1;
		insertBwt = seq;
		relativeRank = seq + bwtInc->buildSize + 1;

		// Store the first character of this iteration
		firstCharInThisIteration = bwtInc->packedText[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

		// Count occurrence of input text
		ForwardDNAAllOccCountNoLimit(bwtInc->packedText, numChar, bwtInc->cumulativeCountInCurrentBuild + 1, bwtInc->bwt->decodeTable);
		// Add the first character of the previous iteration to represent the inverseSa0 of the previous iteration
		bwtInc->cumulativeCountInCurrentBuild[bwtInc->firstCharInLastIteration + 1]++;
        
        for (i = 2; i <=ALPHABET_SIZE; ++i){
            bwtInc->cumulativeCountInCurrentBuild[i] += bwtInc->cumulativeCountInCurrentBuild[i-1];
        }

		// Get rank of new suffix among processed suffix
		// The seq array is built into ALPHABET_SIZE + 2 groups; ALPHABET_SIZE groups + 1 group divided into 2 by inverseSa0 + inverseSa0 as 1 group
		oldInverseSa0RelativeRank = BWTIncGetAbsoluteRank(bwtInc->bwt, sortedRank, seq, bwtInc->packedText,
														  numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->firstCharInLastIteration);

		// Sort rank by ALPHABET_SIZE + 2 groups (or ALPHABET_SIZE + 1 groups when inverseSa0 sit on the border of a group)
		for (i=0; i<ALPHABET_SIZE; i++) {
			if (bwtInc->cumulativeCountInCurrentBuild[i] > oldInverseSa0RelativeRank ||
				bwtInc->cumulativeCountInCurrentBuild[i+1] <= oldInverseSa0RelativeRank) {
				BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], bwtInc->cumulativeCountInCurrentBuild[i+1] - bwtInc->cumulativeCountInCurrentBuild[i]);
			} else {
				if (bwtInc->cumulativeCountInCurrentBuild[i] < oldInverseSa0RelativeRank) {
					BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], oldInverseSa0RelativeRank - bwtInc->cumulativeCountInCurrentBuild[i]);
				}
				if (bwtInc->cumulativeCountInCurrentBuild[i+1] > oldInverseSa0RelativeRank + 1) {
					BWTIncSortKey(sortedRank + oldInverseSa0RelativeRank + 1, seq + oldInverseSa0RelativeRank + 1, bwtInc->cumulativeCountInCurrentBuild[i+1] - oldInverseSa0RelativeRank - 1);
				}
			}
		}

		// build relative rank; sortedRank is updated for merging to cater for the fact that $ is not encoded in bwt
		// the cumulative freq information is used to make sure that inverseSa0 and suffix beginning with different characters are kept in different unsorted groups)
		BWTIncBuildRelativeRank(sortedRank, seq, relativeRank, numChar, bwtInc->bwt->inverseSa0, bwtInc->cumulativeCountInCurrentBuild);
#ifdef DEBUG
		if (relativeRank[numChar] != oldInverseSa0RelativeRank) {
			fprintf(stderr, "BWTIncConstruct(): relativeRank[numChar] != oldInverseSa0RelativeRank!\n");
			exit(1);
		}
#endif

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)numChar, 1, TRUE);

		newInverseSa0RelativeRank = relativeRank[0];
		newInverseSa0 = sortedRank[newInverseSa0RelativeRank] + newInverseSa0RelativeRank;

		sortedRank[newInverseSa0RelativeRank] = 0;	// a special value so that this is skipped in the merged bwt

		// Build BWT
		BWTIncBuildBwt(seq, relativeRank, numChar, bwtInc->cumulativeCountInCurrentBuild);

		// Merge BWT
		mergedBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord
				    - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR;
					// minus numberOfIteration * occInterval to create a buffer for merging
		BWTIncMergeBwt(sortedRank, bwtInc->bwt->bwtCode, insertBwt, mergedBwt, bwtInc->bwt->textLength, numChar);

	}

	// Build auxiliary structure and update info and pointers in BWT
	bwtInc->bwt->textLength += numChar;
	bwtInc->bwt->bwtCode = mergedBwt;
	bwtInc->bwt->bwtSizeInWord = mergedBwtSizeInWord;
	bwtInc->bwt->occSizeInWord = mergedOccSizeInWord;
	if (mergedBwt < bwtInc->workingMemory + mergedOccSizeInWord) {
		fprintf(stderr, "BWTIncConstruct() : Not enough memory allocated!\n");
		exit(1);
	}

	bwtInc->bwt->occValue = mergedBwt - mergedOccSizeInWord;

	BWTClearTrailingBwtCode(bwtInc->bwt);
	BWTGenerateOccValueFromBwt(bwtInc->bwt->bwtCode, bwtInc->bwt->occValue, bwtInc->bwt->occValueMajor,
							   bwtInc->bwt->textLength, bwtInc->bwt->decodeTable);

	bwtInc->bwt->inverseSa0 = newInverseSa0;
   
	for (i=1; i<=ALPHABET_SIZE; i++) {
        bwtInc->bwt->cumulativeFreq[i] += bwtInc->cumulativeCountInCurrentBuild[i] - (bwtInc->firstCharInLastIteration <= i-1);
    }

	bwtInc->firstCharInLastIteration = firstCharInThisIteration;

	// Set build size and text address for the next build
	BWTIncSetBuildSizeAndTextAddr(bwtInc);
	bwtInc->numberOfIterationDone++;

}

static BWTInc *BWTIncConstructFromPacked(const char *inputFileName, const float targetNBit,
								  const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize)
{

	FILE *packedFile;
	unsigned int packedFileLen;
	unsigned int totalTextLength;
	unsigned int textToLoad, textSizeInByte;
	unsigned int processedTextLength;
	unsigned char lastByteLength;

	BWTInc *bwtInc;

	packedFile = (FILE*)fopen(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "BWTIncConstructFromPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftell(packedFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "BWTIncConstructFromPacked: Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);
	totalTextLength = TextLengthFromBytePacked(packedFileLen, BIT_PER_CHAR, lastByteLength);

	bwtInc = BWTIncCreate(totalTextLength, targetNBit, initialMaxBuildSize, incMaxBuildSize);

	BWTIncSetBuildSizeAndTextAddr(bwtInc);

	if (bwtInc->buildSize > totalTextLength) {
		textToLoad = totalTextLength;
	} else {
		textToLoad = totalTextLength - ((totalTextLength - bwtInc->buildSize + CHAR_PER_WORD - 1) / CHAR_PER_WORD * CHAR_PER_WORD);
	}
	textSizeInByte = textToLoad / CHAR_PER_BYTE;	// excluded the odd byte

	fseek(packedFile, -2, SEEK_CUR);
	fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
	fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte + 1, packedFile);
	fseek(packedFile, -((int)textSizeInByte + 1), SEEK_CUR);

	ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
	BWTIncConstruct(bwtInc, textToLoad);

	processedTextLength = textToLoad;

	while (processedTextLength < totalTextLength) {
		textToLoad = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;
		if (textToLoad > totalTextLength - processedTextLength) {
			textToLoad = totalTextLength - processedTextLength;
		}
		textSizeInByte = textToLoad / CHAR_PER_BYTE;
		fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
		fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte, packedFile);
		fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
		ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
		BWTIncConstruct(bwtInc, textToLoad);
		processedTextLength += textToLoad;
		if (bwtInc->numberOfIterationDone % 10 == 0) {
			printf("[BWTIncConstructFromPacked] %u iterations done. %u characters processed.\n",
				   bwtInc->numberOfIterationDone, processedTextLength);
		}
	}
	return bwtInc;
}

static void BWTFree(BWT *bwt)
{
	if (bwt == 0) return;
	free(bwt->cumulativeFreq);
	free(bwt->bwtCode);
	free(bwt->occValue);
	free(bwt->occValueMajor);
	free(bwt->saValue);
	free(bwt->inverseSa);
	free(bwt->decodeTable);
	free(bwt->saIndexRange);
	free(bwt->saValueOnBoundary);
	free(bwt);
}

static void BWTIncFree(BWTInc *bwtInc)
{
	if (bwtInc == 0) return;
	free(bwtInc->bwt);
	free(bwtInc->workingMemory);
	free(bwtInc);
}

static unsigned int BWTFileSizeInWord(const unsigned int numChar)
{
	// The $ in BWT at the position of inverseSa0 is not encoded
	return (numChar + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
}

static void BWTSaveBwtCodeAndOcc(const BWT *bwt, const char *bwtFileName, const char *occValueFileName)
{
	FILE *bwtFile;
	FILE *occValueFile; 
	unsigned int bwtLength;

	bwtFile = (FILE*)fopen(bwtFileName, "wb");
	if (bwtFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open BWT code file!\n");
		exit(1);
	}
    fwrite(&bwt->textLength, sizeof(uint32_t), 1, bwtFile);
	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, bwtFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, bwtFile);
	//bwtLength = BWTFileSizeInWord(bwt->textLength);
	//fwrite(&bwtLength, sizeof(uint32_t), 1, bwtFile);
    //fwrite(bwt->bwtCode, sizeof(unsigned int), bwtLength, bwtFile);
	fwrite(&bwt->bwtSizeInWord, sizeof(uint32_t), 1, bwtFile);
    fwrite(bwt->bwtCode, sizeof(unsigned int), bwt->bwtSizeInWord, bwtFile);
    fclose(bwtFile);

	occValueFile = (FILE*)fopen(occValueFileName, "wb");
	if (occValueFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open occ value file!\n");
		exit(1);
	}

//	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, occValueFile);
// 	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, occValueFile);
	fwrite(&bwt->occSizeInWord, sizeof(uint32_t), 1, occValueFile);
    fwrite(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
	fwrite(&bwt->occMajorSizeInWord, sizeof(uint32_t), 1, occValueFile);
	fwrite(bwt->occValueMajor, sizeof(unsigned int), bwt->occMajorSizeInWord, occValueFile);
	fclose(occValueFile);

}

void Rbwt_bwt_bwtgen(const char *fn_pac, const char *fn)
{
	double targetNBit = 6.0;//for 4bit  targetNBit should be seted 2.5
    int initMaxBuildSize = 10000000;
    int incMaxBuildSize = 10000000;
    char fn_bwt[128], fn_Occ[128];
    strcpy(fn_bwt, fn);strcat(fn_bwt, ".bwt");
    strcpy(fn_Occ, fn);strcat(fn_Occ, ".occ");
    BWTInc *bwtInc;
    bwtInc = BWTIncConstructFromPacked(fn_pac, targetNBit, initMaxBuildSize, incMaxBuildSize);
	printf("[bwt_gen] Finished constructing BWT in %u iterations.\n", bwtInc->numberOfIterationDone);
	BWTSaveBwtCodeAndOcc(bwtInc->bwt, fn_bwt, fn_Occ);
	BWTIncFree(bwtInc);
}

int Rbwt_bwt_bwtgen_main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: bwtgen <in.pac> <out.prefix>\n");
		return 1;
	}
	Rbwt_bwt_bwtgen(argv[1], argv[2]);
	return 0;
}

#ifdef MAIN_BWT_GEN
int main(int argc, char *argv[])
{
	return Rbwt_bwt_bwtgen_main(argc, argv);
}
#endif
