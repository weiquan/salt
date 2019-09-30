/*
 * =====================================================================================
 *
 *       Filename:  variant.c
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
#include "variant.h"
#define RNAME_SIZE 32  
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

variant_t *variant_init(FILE *fp)
{
    char tmp[TMP_SIZE];
    variant_t *v = calloc(1, sizeof(variant_t));
    v->rname = malloc(RNAME_SIZE); 
    v->fp_variant = fp;  
    if(v->fp_variant == NULL){
        fprintf(stderr, "[variant_init] : fp_variant == NULL!\n");
        exit(1);
    }
    fgets(tmp, TMP_SIZE, v->fp_variant);//skip 1st line  
    return v;
}
void variant_destroy(variant_t *v)
{
    free(v->rname);
    free(v->snp_pos);
    free(v->snp_type);
}
int variant_get_snpnum(variant_t *v)
{
    char tmp[TMP_SIZE];
    char *rname_cur;
    char rname_last[RNAME_SIZE] = "";
    int snp_num = 0;
    long pos_origin = ftell(v->fp_variant);
    
    if (fgets(tmp, TMP_SIZE, v->fp_variant) != NULL ) { 
        rname_cur = (char *)strtok(tmp, "\t");
        ++snp_num;
        strncpy(rname_last, rname_cur, RNAME_SIZE); 
    } else {
        return snp_num;       
    }
    while ( fgets(tmp, TMP_SIZE, v->fp_variant) != NULL )
    {
        rname_cur = strtok(tmp, "\t");        
        if (strcmp(rname_cur, rname_last) != 0 ){
            break;
        }
        ++snp_num;
        strncpy(rname_last, rname_cur, RNAME_SIZE);
    }
    fseek(v->fp_variant, pos_origin, SEEK_SET);
    return snp_num;
}
int variant_readv(variant_t *v)
{
    int i, j;
    char tmp[TMP_SIZE];
    const char *str_buf;
    int char_nt;
    uint8_t snp_type ;

    int snp_num = variant_get_snpnum(v);
    if (snp_num > v->snp_num){
        v->snp_pos = (uint32_t *)realloc(v->snp_pos, snp_num *sizeof(uint32_t));
        if(v->snp_pos == NULL){
            fprintf(stderr, "[variant_readv]:realloc snp_pos fail!\n");
            exit(1);
        }
        v->snp_type = (uint8_t *)realloc(v->snp_type, snp_num*sizeof(uint8_t)); 
        if(v->snp_type == NULL){
            fprintf(stderr, "[variant_readv]:realloc snp_pos fail!\n");
            exit(1);
        }
    }
    v->snp_num = snp_num;    
    if (snp_num == 0){
        fprintf(stderr, "[variant_readv] :snp_num == 0\n");
        return -1;
    }
    //get rname and 1st snp_type   
    snp_type = 0;
    if (fgets(tmp, TMP_SIZE, v->fp_variant) == NULL){
        fprintf(stderr, "[variant_readv]: fgets == NULL\n");
        exit(0);
    }
    strncpy(v->rname, strtok(tmp, "\t"), RNAME_SIZE);
    v->snp_pos[0] = atoi(strtok(NULL, "\t")) -1;
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
    v->snp_type[0] = snp_type;
    //rest snp_type
    for(i = 1; i < snp_num; ++i)
    {
        snp_type = 0;
        if (fgets(tmp, TMP_SIZE, v->fp_variant) == NULL){
            fprintf(stderr, "[variant_readv]: fgets == NULL\n");
            exit(0);
        }
        strtok(tmp, "\t"); //skip rname
        v->snp_pos[i] = atoi(strtok(NULL, "\t")) -1;
        str_buf = strtok(NULL, "\t");//snp
        for(j = 0; j < strlen(str_buf); j += 2)
        {
            char_nt = nst_nt4_table[(uint32_t)str_buf[j]];
            snp_type = snp_type|(1<<char_nt);
        } 
        str_buf = strtok(NULL, "\t");//ref
        char_nt = nst_nt4_table[(uint32_t)str_buf[0]];
        snp_type = snp_type|(char_nt<<4);
        v->snp_type[i] = snp_type;
    }   
    return 0;
}
/*
int main()
{
    variant_t *v;
    FILE *fp = fopen("mart_.txt", "r");
    v = variant_init(fp);
    //read 1st Chromosome
    variant_readv(v);
    printf("chr_name:%s\n",v->rname);
    int i;
    for(i = 0; i < v->snp_num; ++i)
    {
        printf("[%d] :%d\t", i, v->snp_pos[i]);
        printf("%d\n", v->snp_type[i]);
    }
    //read 2nd Chromosome
    variant_readv(v);
    printf("chr_name:%s\n",v->rname);
    for(i = 0; i < v->snp_num; ++i)
    {
        printf("[%d] :%d\n", i, v->snp_pos[i]);
        printf("%d\n", v->snp_type[i]);

    }
    variant_destroy(v);
    // 
    return 0;
}*/
