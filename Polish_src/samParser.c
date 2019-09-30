/*
 * =====================================================================================
 *
 *       Filename:  samParser.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/27/2014 04:18:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <string.h>
#include "samParser.h"
#include "kstring.h"
static unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
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


void sam_skipHeader(FILE *fp)
{
    char buffer[1024] = {0};
    while(fgets(buffer, 1024, fp) != NULL){
        if(buffer[0] != '@'){
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        } 
    
    }


}
char *readline(FILE *fp)
{
    int l_line = 0;
    int m_line = 1024;
    char *line = calloc(m_line, sizeof(char));
    int c = fgetc(fp);
    while( c != EOF && c != '\n'){
        if(l_line == m_line) {
            m_line <<= 1; 
        }
        line = realloc(line, m_line);
        line[l_line++] = c;
        c = fgetc(fp);
    }
    if(l_line == m_line) { 
        m_line <<= 1;
        line = realloc(line, m_line);
    }
    line[l_line] = '\0';

    return line;
}
const char *DELIM_TAB = "\t";
const char *DELIM_COMMA = ",";
const char *DELIM_SEMICOLON = ";";
const char *DELIM_COLON = ":";


SAM_t *sam_readline(FILE *fp)
{
    int i; 

    char *line = readline(fp);
    if(strlen(line) == 0){
        free(line);
        return NULL;
    }
    SAM_t *sam = calloc(1, sizeof(SAM_t)); 
    sam->buf = line;
    multi_hits_t *h0 = sam->multi_hits[0] = calloc(1, sizeof(multi_hits_t));
    multi_hits_t *h1 = sam->multi_hits[1] = calloc(1, sizeof(multi_hits_t));
    kv_resize(hit_t, *h0, 32); 
    kv_resize(hit_t, *h1, 32); 

        
    sam->seq_name = strtok(line, DELIM_TAB);// QNAME 
    sam->flag = (unsigned int)atoi(strtok(NULL, DELIM_TAB)); // FLAG
    char *chrom = strtok(NULL, DELIM_TAB); //RNAME
    unsigned int pos = (unsigned int)strtoul(strtok(NULL, DELIM_TAB), NULL, 10); //POS
    if((sam->flag &0x4) == 0 && strcmp(chrom, "*") != 0){ 
        if((sam->flag&0x0010) == 0 ){//forward
            if(h0->n == h0->m) {
                h0->m <<=1;
                h0->a = realloc(h0->a, h0->m * sizeof(hit_t));
            }
            h0->a[h0->n].pos = pos;
            h0->a[h0->n].chrom = chrom;
            ++h0->n;
    
        } else{//backward
            if(h1->n == h1->m) {
                h1->m <<=1;
                h1->a = realloc(h1->a, h1->m * sizeof(hit_t));
            }
            h1->a[h1->n].pos = pos;
            h1->a[h1->n].chrom = chrom;
            ++h1->n;
        }
    }
    sam->mapq = atoi(strtok(NULL, DELIM_TAB));//MAPQ
    //sam->cigar = strtok(NULL, DELIM_TAB);//CIGAR
    sam->b0 = -100000;
    sam->b1 = -100000;
    char *cigar = strtok(NULL, DELIM_TAB);//CIGAR
    sam->cigar.l = 0;
    sam->cigar.m = 1024;
    sam->cigar.s = calloc(1024, 1); 
    sam->mate_seqname = strtok(NULL, DELIM_TAB); //MRNM
    sam->mate_pos = (unsigned int)strtoul(strtok(NULL, DELIM_TAB), NULL, 10); //MPOS
    sam->isize = (unsigned int)strtoul(strtok(NULL, DELIM_TAB), NULL, 10); //ISIZE
    sam->seq = strtok(NULL, DELIM_TAB); sam->l_seq = strlen(sam->seq);
    sam->nst_seq = calloc((sam->l_seq+15)/8*8, 1);
    sam->nst_rseq = calloc((sam->l_seq+15)/8*8, 1);
    for(i=0; i< sam->l_seq; ++i) sam->nst_seq[i] = nst_nt4_table[(unsigned char)sam->seq[i]];
    for(i=0; i < sam->l_seq; ++i) { 
        sam->nst_rseq[i] = 3 - sam->nst_seq[sam->l_seq-i-1];
    }
    if(sam->flag &0x0010) {
        unsigned char *tmp = sam->nst_seq;
        sam->nst_seq = sam->nst_rseq;
        sam->nst_rseq = tmp;

    } 
    sam->qual = strtok(NULL, DELIM_TAB);
    char *opt = strtok(NULL, DELIM_TAB);
    while(opt != NULL){
        if(strstr(opt, "XA") != NULL){
            char *multi_aln = strtok(opt, DELIM_COLON);
            multi_aln = strtok(NULL, DELIM_COLON);
            multi_aln = strtok(NULL, DELIM_COLON);
            
            char *aln = strtok(multi_aln, DELIM_SEMICOLON);
            while(aln != NULL) {
                int l_aln = strlen(aln);
		if(l_aln == 0){ break;}
                char *aln_chrom = strtok(aln, DELIM_COMMA);
                char *aln_pos = strtok(NULL, DELIM_COMMA);
		//fprintf(stderr, "[Aln]:%s %s\n", aln_chrom, aln_pos);
		
		if(aln_pos[0] != '-'){
                    if(h0->n == h0->m) {
                        


			h0->m <<=1;
                        h0->a = realloc(h0->a, h0->m * sizeof(hit_t));
                    }
                    h0->a[h0->n].pos = (unsigned int)strtoul(aln_pos, NULL, 10);
                    h0->a[h0->n].chrom = aln_chrom;
                    ++h0->n;
            
                } else{
                    if(h1->n == h1->m) {
                        h1->m <<=1;
                        h1->a = realloc(h1->a, h1->m * sizeof(hit_t));
                    }
                    h1->a[h1->n].pos = (unsigned int)strtoul(aln_pos+1, NULL, 10);
                    h1->a[h1->n].chrom = aln_chrom;
                    ++h1->n;

                } 
                multi_aln += l_aln +1; 
                aln = strtok(multi_aln, DELIM_SEMICOLON);
            }
        }
        
        opt = strtok(NULL, DELIM_TAB);
    } 
     
    return sam;    

}
void sam_destroy(SAM_t *sam)
{
    free(sam->cigar.s);
    free(sam->multi_hits[0]->a);
    free(sam->multi_hits[1]->a);
    free(sam->multi_hits[0]);
    free(sam->multi_hits[1]);
    free(sam->nst_seq);
    free(sam->nst_rseq);
    free(sam->buf);
    free(sam);
}
