/*
 * =====================================================================================
 *
 *       Filename:  snpSeg.c
 *
 *    Description:  用来生成所有中间位置为snp的segment。
 *                  segment长49bp，中间元素为snp。  
 *
 *        Version:  1.0
 *        Created:  09/22/2012 02:04:11 PM
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
#include <zlib.h>
#include <time.h>
#include "kseq.h"
#include "hapmap.h"
KSEQ_INIT(gzFile, gzread)
#define WIN_MAX_SNP_NUM 5
int ss_core(const char *fn_ref, const char *fn_hm, int l_seed, const char *fn_out)
{
    clock_t ss_startTime, ss_endTime;  
    ss_startTime = clock(); 
    int i, j;
    gzFile fp_ref = Z_NULL;FILE *fp_hm = NULL, *fp_out= NULL;
    kseq_t *seq; hapmap_t *hm;
    const int WIN_SNP_DISTANCE = l_seed-1;//窗口中核心snp距离最左端的距离,也是最右端距离
    //init 
    fp_ref = gzopen(fn_ref, "r");
    if (fp_ref == Z_NULL){
        fprintf(stderr, "gzopen %s fail!\n", fn_ref);
        exit(1);
    }
    fp_hm = fopen(fn_hm, "r");
    if (fp_hm == NULL){
        fprintf(stderr, "fopen %s fail!\n", fn_hm);
        exit(1);
    }
    fp_out = fopen(fn_out, "w");
    if (fp_out == NULL){
        fprintf(stderr, "fopen %s fail!\n", fn_out);
        exit(1);
    }
    
    seq = kseq_init(fp_ref);
    hm = hapmap_init(fp_hm); 
     
    //core loop
    int l; uint32_t tot_l = 0;
    uint32_t snp_tot_num = 0;
    //uint32_t segment_refid;
    while( kseq_read(seq) > 0 )
    {
        l = seq->seq.l;
        uint32_t *snp_pos, snp_num; uint8_t *snp_type;

        int win_start, win_end;
        int win_snp_start,win_snp_end, win_snp_num, win_snp_mid;       
        int win_segment_num, segment_num_factor1, segment_num_factor2; 
        char snp_type_cur;
        
        fprintf(stderr, "[GLP]: generate local pattern in %s.\n", seq->name.s);

        hapmap_readhm(hm);
        snp_num = hm->snp_num;
        if(snp_num == 0){
            fprintf(stderr, "[GLP]: no snp in %s.\n", seq->name.s); 
            continue;
        }
        snp_pos = hm->snp_pos;
        snp_type = hm->snp_type;  
        if ( strcmp(hm->chrID, seq->name.s) != 0 ){
            fprintf(stderr, "[GLP]: ref chrID not match hapmap chrID!\n"); 
            continue;
        } 

        win_snp_mid = 0;
        while(win_snp_mid < snp_num){
            //compute snp num in current window
            // [win_snp_start, win_snp_end)
            win_snp_start = win_snp_mid-1;
            win_snp_end = win_snp_mid+1;
            while(win_snp_start>=0&& snp_pos[win_snp_mid] - snp_pos[win_snp_start] <= WIN_SNP_DISTANCE){
                --win_snp_start;
            }
            win_snp_start++;
            while( win_snp_end <=snp_num && snp_pos[win_snp_end]- snp_pos[win_snp_mid] <= WIN_SNP_DISTANCE ){
                ++win_snp_end;
            }

            win_snp_num = win_snp_end - win_snp_start;
            if(win_snp_num > WIN_MAX_SNP_NUM){
                fprintf(stderr, "[skip variant]: %u\n", snp_pos[win_snp_mid]);
                win_snp_mid++;
                continue; 
            }
            //compute window boundery
            //win_start = snp_pos[win_snp_start] - WIN_SNP_DISTANCE;
            //win_end = snp_pos[win_snp_end-1] +WIN_SNP_DISTANCE; 
            win_start = (snp_pos[win_snp_mid] > WIN_SNP_DISTANCE)? (snp_pos[win_snp_mid]- WIN_SNP_DISTANCE) : 0;
            win_end = snp_pos[win_snp_mid] + WIN_SNP_DISTANCE < l ?
                      snp_pos[win_snp_mid] + WIN_SNP_DISTANCE : l-1;                  ; 
            //compute segment num in current window

            win_segment_num = hapmap_get_snptypenum(snp_type[win_snp_start]); 
            for(i = 1; i < win_snp_num; ++i){
                win_segment_num *= hapmap_get_snptypenum(snp_type[win_snp_start+i]);
            }
//            fprintf(fp_out, ">%d\t%u\t%u\n", snp_tot_num++, snp_pos[win_snp_mid]+tot_l -24, snp_pos[win_snp_mid]);         
            fprintf(fp_out, ">%d\t%u\n", snp_tot_num++, snp_pos[win_snp_mid]+tot_l +WIN_SNP_DISTANCE);         
            
            if(snp_tot_num == 1){
                fputc('#', fp_out); 
            }
            //遍历该window中每个segment 
           
             
            for(i = 0; i < win_segment_num; ++i){
                int snp_type_i;
                int k = i; 
                
            //if (i == segment_refid) continue; 
                //generate cur segment sequence 
                segment_num_factor1 = 1; 
                for(j = 0; j < win_snp_num; ++j){

                    snp_type_cur = snp_type[win_snp_start+j];
                    segment_num_factor1 *= hapmap_get_snptypenum(snp_type_cur);
                    segment_num_factor2 = win_segment_num / segment_num_factor1; 
                    snp_type_i = k/segment_num_factor2;
                    //skip ref nt
                    if( win_snp_mid - win_snp_start == j &&hapmap_nt2snptypei(snp_type_cur, snp_type_cur>>4) == snp_type_i){
                            break;
                    }
                    k -= snp_type_i * segment_num_factor2;
                    seq->seq.s[snp_pos[win_snp_start+j]] = "ACGTN"[hapmap_get_snptype(snp_type_cur, snp_type_i)]; 
                }
                //output current segment
                //fprintf(fp_out, ">Seg%d\t%d\t%s\n", snp_tot_num++, hm->snp_pos[win_snp_start], hm->chrID);
                if( win_snp_mid - win_snp_start == j &&hapmap_nt2snptypei(snp_type_cur, snp_type_cur>>4) == snp_type_i){
                    continue;
                }
                for(j = win_start; j <= win_end; ++j){
                    fputc(seq->seq.s[j], fp_out);
                }
                fputc('#', fp_out);fputc('\n', fp_out);
            } 
            ++win_snp_mid;
        }//end while snp oteration    
        tot_l += l;
    }//end while chrome iteration
    //destroy 
    kseq_destroy(seq);
    hapmap_destroy(hm);
    {
        gzclose(fp_ref);
        fclose(fp_hm);
        fclose(fp_out);
    } 
    ss_endTime = clock();
    fprintf(stderr, "[GLP]:time escaped %lu sec.\n", (ss_endTime - ss_startTime)/CLOCKS_PER_SEC);
 
}
int ss_core_alt(const char *fn_ref, const char *fn_hm, int l_seed, const char *fn_out)
{
    clock_t ss_startTime, ss_endTime;  
    ss_startTime = clock(); 
    int i, j;
    gzFile fp_ref = Z_NULL;FILE *fp_hm = NULL, *fp_out= NULL;
    kseq_t *seq; hapmap_t *hm;
    const int WIN_SNP_DISTANCE = l_seed-1;//窗口中核心snp距离最左端的距离,也是最右端距离
    //init 
    fp_ref = gzopen(fn_ref, "r");
    if (fp_ref == Z_NULL){
        fprintf(stderr, "gzopen %s fail!\n", fn_ref);
        exit(1);
    }
    fp_hm = fopen(fn_hm, "r");
    if (fp_hm == NULL){
        fprintf(stderr, "fopen %s fail!\n", fn_hm);
        exit(1);
    }
    fp_out = fopen(fn_out, "w");
    if (fp_out == NULL){
        fprintf(stderr, "fopen %s fail!\n", fn_out);
        exit(1);
    }
    
    seq = kseq_init(fp_ref);
    hm = hapmap_init(fp_hm); 
     
    //core loop
    int l; uint32_t tot_l = 0;
    uint32_t snp_tot_num = 0;
    //uint32_t segment_refid;
    while( kseq_read(seq) > 0 )
    {
        l = seq->seq.l;
        uint32_t *snp_pos, snp_num; uint8_t *snp_type;

        int win_start, win_end;
        int win_snp_start,win_snp_end, win_snp_num, win_snp_mid;       
        int win_segment_num, segment_num_factor1, segment_num_factor2; 
        char snp_type_cur;
        
        fprintf(stderr, "[GLP]: generate local pattern in %s.\n", seq->name.s);

        hapmap_readhm(hm);
        snp_num = hm->snp_num;
        if(snp_num == 0){
            fprintf(stderr, "[GLP]: no snp in %s.\n", seq->name.s); 
            continue;
        }
        snp_pos = hm->snp_pos;
        snp_type = hm->snp_type;  
        if ( strcmp(hm->chrID, seq->name.s) != 0 ){
            fprintf(stderr, "[GLP]: ref chrID not match hapmap chrID!\n"); 
            continue;
        } 

        win_snp_mid = 0;
        while(win_snp_mid < snp_num){
            //compute snp num in current window
            // [win_snp_start, win_snp_end)
            win_snp_start = win_snp_mid;
            win_snp_end = win_snp_mid+1;
            /*
	    while(win_snp_start>=0&& snp_pos[win_snp_mid] - snp_pos[win_snp_start] <= WIN_SNP_DISTANCE){
                --win_snp_start;
            }

            win_snp_start++;
            */
	    while( win_snp_end <=snp_num && snp_pos[win_snp_end]- snp_pos[win_snp_mid] <= WIN_SNP_DISTANCE ){
                ++win_snp_end;
            }

            win_snp_num = win_snp_end - win_snp_start;
            if(win_snp_num > WIN_MAX_SNP_NUM){
                fprintf(stderr, "[skip variant]: %u\n", snp_pos[win_snp_mid]);
                win_snp_mid++;
                continue; 
            }
            //compute window boundery
            //win_start = snp_pos[win_snp_start] - WIN_SNP_DISTANCE;
            //win_end = snp_pos[win_snp_end-1] +WIN_SNP_DISTANCE; 
            win_start = (snp_pos[win_snp_start] > WIN_SNP_DISTANCE)? (snp_pos[win_snp_start]- WIN_SNP_DISTANCE) : 0;
            if(win_snp_start > 0 && snp_pos[win_snp_start] - snp_pos[win_snp_start-1] <= WIN_SNP_DISTANCE){
	    	win_start = snp_pos[win_snp_start-1] + 1;
	    }
	    win_end = snp_pos[win_snp_mid] + WIN_SNP_DISTANCE < l ?
                      snp_pos[win_snp_mid] + WIN_SNP_DISTANCE : l-1;                  ; 
            //compute segment num in current window

            win_segment_num = hapmap_get_snptypenum(snp_type[win_snp_start]); 
            for(i = 1; i < win_snp_num; ++i){
                win_segment_num *= hapmap_get_snptypenum(snp_type[win_snp_start+i]);
            }
//            fprintf(fp_out, ">%d\t%u\t%u\n", snp_tot_num++, snp_pos[win_snp_mid]+tot_l -24, snp_pos[win_snp_mid]);         
            fprintf(fp_out, ">%d_%u\t%u\n", snp_tot_num++, win_segment_num,  snp_pos[win_snp_mid]+tot_l +WIN_SNP_DISTANCE);         
            
            if(snp_tot_num == 1){
                fputc('#', fp_out); 
            }
            //遍历该window中每个segment 
           
            uint32_t segment_refid; 
            for(i = 0; i < win_segment_num; ++i){
                int snp_type_i;
                int k = i; 
                
            //if (i == segment_refid) continue; 
                //generate cur segment sequence 
                segment_num_factor1 = 1; 
                for(j = 0; j < win_snp_num; ++j){

                    snp_type_cur = snp_type[win_snp_start+j];
                    segment_num_factor1 *= hapmap_get_snptypenum(snp_type_cur);
                    segment_num_factor2 = win_segment_num / segment_num_factor1; 
                    snp_type_i = k/segment_num_factor2;
                    //skip ref nt
                    /*
		    if( win_snp_mid - win_snp_start == j &&hapmap_nt2snptypei(snp_type_cur, snp_type_cur>>4) == snp_type_i){
                            break;
                    }
		    */
                    k -= snp_type_i * segment_num_factor2;
                    seq->seq.s[snp_pos[win_snp_start+j]] = "ACGTN"[hapmap_get_snptype(snp_type_cur, snp_type_i)]; 
                }
                //output current segment
                //fprintf(fp_out, ">Seg%d\t%d\t%s\n", snp_tot_num++, hm->snp_pos[win_snp_start], hm->chrID);
                /*
		if( win_snp_mid - win_snp_start == j &&hapmap_nt2snptypei(snp_type_cur, snp_type_cur>>4) == snp_type_i){
                    continue;
                }
		*/
                for(j = win_start; j <= win_end; ++j){
                    fputc(seq->seq.s[j], fp_out);
                }
                fputc('#', fp_out);fputc('\n', fp_out);
            } 
            ++win_snp_mid;
        }//end while snp oteration    
        tot_l += l;
    }//end while chrome iteration
    //destroy 
    kseq_destroy(seq);
    hapmap_destroy(hm);
    {
        gzclose(fp_ref);
        fclose(fp_hm);
        fclose(fp_out);
    } 
    ss_endTime = clock();
    fprintf(stderr, "[GLP]:time escaped %lu sec.\n", (ss_endTime - ss_startTime)/CLOCKS_PER_SEC);
 
}

int ss_main(int argc, char *argv[])
{
    int l_seed = 25; 
    if (argc < 4){
        fprintf(stderr, "Usage: ss3 <in.fasta> <in.hapmap> <out>\n");
        return 1;
    }
    ss_core_alt(argv[1], argv[2], l_seed, argv[3]);
    return 0;
}
#ifdef MAIN_SS
int main(int argc, char *argv[])
{
    return ss_main(argc, argv);
}

#endif 
