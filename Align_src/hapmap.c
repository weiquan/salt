/*
 * =====================================================================================
 *
 *       Filename:  hapmap.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/14/2012 08:04:55 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <string.h>
#include "hapmap.h"
#define CHRID_SIZE 32  
#define TMP_SIZE 128
static unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

hapmap_t *hapmap_init(FILE *fp)
{
    char tmp[TMP_SIZE];
    hapmap_t *hm = calloc(1, sizeof(hapmap_t));
    hm->chrID = malloc(CHRID_SIZE); 
    hm->fp_hapmap = fp;  
    if(hm->fp_hapmap == NULL){
        fprintf(stderr, "[hapmap_init] : fp_hapmap == NULL!\n");
        exit(1);
    }
    fgets(tmp, TMP_SIZE, hm->fp_hapmap);//skip 1st line  
    return hm;
}
void hapmap_destroy(hapmap_t *hm)
{
    free(hm->chrID);
    free(hm->snp_pos);
    free(hm->snp_type);
}
int hapmap_get_snpnum(hapmap_t *hm)
{
    char tmp[TMP_SIZE];
    char *chrID_cur;
    char chrID_last[CHRID_SIZE] = "";
    int snp_num = 0;
    long pos_origin = ftell(hm->fp_hapmap);
    
    if (fgets(tmp, TMP_SIZE, hm->fp_hapmap) != NULL ) { 
        chrID_cur = (char *)strtok(tmp, "\t");
        ++snp_num;
        strncpy(chrID_last, chrID_cur, CHRID_SIZE); 
    } else {
        return snp_num;       
    }
    while ( fgets(tmp, TMP_SIZE, hm->fp_hapmap) != NULL )
    {
        chrID_cur = strtok(tmp, "\t");        
        if (strcmp(chrID_cur, chrID_last) != 0 ){
            break;
        }
        ++snp_num;
        strncpy(chrID_last, chrID_cur, CHRID_SIZE);
    }
    fseek(hm->fp_hapmap, pos_origin, SEEK_SET);
    return snp_num;
}
int hapmap_readhm(hapmap_t *hm)
{
    int i, j;
    char tmp[TMP_SIZE];
    const char *str_buf;
    int char_nt;
    uint8_t snp_type ;

    int snp_num = hapmap_get_snpnum(hm);
    if (snp_num > hm->snp_num){
        hm->snp_pos = (uint32_t *)realloc(hm->snp_pos, snp_num *sizeof(uint32_t));
        if(hm->snp_pos == NULL){
            fprintf(stderr, "[hapmap_readhm]:realloc snp_pos fail!\n");
            exit(1);
        }
        hm->snp_type = (uint8_t *)realloc(hm->snp_type, snp_num*sizeof(uint8_t)); 
        if(hm->snp_type == NULL){
            fprintf(stderr, "[hapmap_readhm]:realloc snp_pos fail!\n");
            exit(1);
        }
    }
    hm->snp_num = snp_num;    
    if (snp_num == 0){
        fprintf(stderr, "[hapmap_readhm] :snp_num == 0\n");
        return -1;
    }
    //get chrID and 1st snp_type   
    snp_type = 0;
    if (fgets(tmp, TMP_SIZE, hm->fp_hapmap) == NULL){
        fprintf(stderr, "[hammap_readhm]: fgets == NULL\n");
        exit(0);
    }
    strncpy(hm->chrID, strtok(tmp, "\t"), CHRID_SIZE);
    hm->snp_pos[0] = atoi(strtok(NULL, "\t")) -1;
    str_buf = strtok(NULL, "\t");//snp
    for(j = 0; j < strlen(str_buf); j += 2)
    {
        char_nt = nst_nt4_table[(uint32_t)str_buf[j]];
        if(char_nt > 3)
        {
            fprintf(stderr, "#$^^\n");
            exit(1);
        }
        snp_type = snp_type|(1<<char_nt);
    } 
    str_buf = strtok(NULL, "\t");//ref
    char_nt = nst_nt4_table[(uint32_t)str_buf[0]];
    snp_type = snp_type|(char_nt<<4);
    hm->snp_type[0] = snp_type;
    //rest snp_type
    for(i = 1; i < snp_num; ++i)
    {
        snp_type = 0;
        if (fgets(tmp, TMP_SIZE, hm->fp_hapmap) == NULL){
            fprintf(stderr, "[hammap_readhm]: fgets == NULL\n");
            exit(0);
        }
        strtok(tmp, "\t"); //skip chrID
        hm->snp_pos[i] = atoi(strtok(NULL, "\t")) -1;
        str_buf = strtok(NULL, "\t");//snp
        for(j = 0; j < strlen(str_buf); j += 2)
        {
            char_nt = nst_nt4_table[(uint32_t)str_buf[j]];
            snp_type = snp_type|(1<<char_nt);
        } 
        str_buf = strtok(NULL, "\t");//ref
        char_nt = nst_nt4_table[(uint32_t)str_buf[0]];
        snp_type = snp_type|(char_nt<<4);
        hm->snp_type[i] = snp_type;
    }   
    return 0;
}
/*
int main()
{
    hapmap_t *hm;
    FILE *fp = fopen("mart_.txt", "r");
    hm = hapmap_init(fp);
    //read 1st Chromosome
    hapmap_readhm(hm);
    printf("chr_name:%s\n",hm->chrID);
    int i;
    for(i = 0; i < hm->snp_num; ++i)
    {
        printf("[%d] :%d\t", i, hm->snp_pos[i]);
        printf("%d\n", hm->snp_type[i]);
    }
    //read 2nd Chromosome
    hapmap_readhm(hm);
    printf("chr_name:%s\n",hm->chrID);
    for(i = 0; i < hm->snp_num; ++i)
    {
        printf("[%d] :%d\n", i, hm->snp_pos[i]);
        printf("%d\n", hm->snp_type[i]);

    }
    hapmap_destroy(hm);
    // 
    return 0;
}*/
