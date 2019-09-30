/*
 * =====================================================================================
 *
 *       Filename:  fwrite.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2015年01月19日 16时11分09秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
int main ( int argc, char *argv[] )
{
    FILE *fp = fopen("index.snpaln.R.seedLen", "w");
    int seedLen = 25;
    fwrite(&seedLen, sizeof(int), 1, fp);
    fclose(fp); 
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

