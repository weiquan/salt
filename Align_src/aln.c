/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/05/2012 02:53:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * ===============================================================m======================
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include "aln.h"
#define PACKAGE_VERSION "beta0.1"
#define CONTACT "Quan Wei <wquanhit@gmail.com>"



opt_t *opt_init(){

    fprintf(stderr, "[opt_init]:opt init!\n"); 
    opt_t *opt = (opt_t *)calloc(1, sizeof(opt_t));  
    opt->n_diff = -1;
   
    opt->seed = time(0);
    opt->n_threads = 1;
    opt->fn_index = NULL;
    opt->fn_read1 = NULL;
    opt->fn_read2 = NULL;
    opt->rg_id = NULL;
    opt->use_sw_extend = 0;//only works in PE mode singleton pair
    opt->l_seed = 25;
    opt->l_overlap = -1;
    opt->max_tlen = 550;
    opt->min_tlen = 250;
    opt->cmd = (kstring_t *)calloc(1, sizeof(kstring_t));
    opt->max_seed = 50;//max locate number per bwt range
    opt->max_locate = 1000;//max locate number per bwt range
    opt->max_hits = 50; 
    opt->max_walk = 1000;    
    opt->se = 1;
    opt->l_read = 100;
    opt->print_xa_cigar = 0;
    opt->print_nm_md = 0;
	return opt;
}
void opt_destroy(opt_t *opt)
{

    if(opt->rg_id!=NULL) free(opt->rg_id);
    free(opt->cmd->s); free(opt->cmd);
    free(opt);
}
int usage()
{

    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   snpaln\n"); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
    fprintf(stderr, "Usage:     snpaln [Options] <Index.prefix> <Read_mate1> [Read_mate2]\n\n"); 
    fprintf(stderr, "Options:   -h, --help                   help\n"); 
    fprintf(stderr, "           -t, --threads       <int>    thread\n"); 
    fprintf(stderr, "           -n, --num           <int>    diff number [-1]\n"); 
    fprintf(stderr, "           -g, --group         <str>    read group id\n"); 
    fprintf(stderr, "           -l, --read_length   <int>    read length [100]\n");  
    fprintf(stderr, "           -c, --xa_cigar               print cigar in XA fields [False]\n");    
    fprintf(stderr, "           -d, --md                     print tag NM and MD [False]\n");    
    
    fprintf(stderr, "           -r, --overlap       <int>    overlap length [non-overlap seeding]\n");  
    fprintf(stderr, "           -v, --ref           <int>    only seeding on primary reference  [0]\n");  
    fprintf(stderr, "           -s, --max_seed      <int>    max seed occ per [50]\n");    
    fprintf(stderr, "           -m, --max_locate    <int>    max locate number per bwt range [200]\n");    
    fprintf(stderr, "           -R, --max_best_hits <int>    max best hits  [50]\n");    
    fprintf(stderr, "           -X, --extend        <int>    extend algorithm [lv|sw]\n"); 
    
    
    fprintf(stderr, "           -p, --pe                     paired end mode [SE]\n"); 
    fprintf(stderr, "           -a, --min_tlen      <int>    min insert size [250]\n"); 
    fprintf(stderr, "           -b, --max_tlen      <int>    max insert size [550]\n"); 


    fprintf(stderr, "           -e, --sw                     use smith-watermon when singleton in PE mode [False]\n"); 
    

   
   

    fprintf(stderr, "\n"); 
    return EXIT_SUCCESS;
}
const char * short_options = "t:n:hpa:b:g:em:s:l:cdr:vM:O:E:X:";
struct option long_options[] = {
    { "threads",     1,   NULL,    't'   },
    { "num",     1,   NULL,    'n' },
    { "help",     0,   NULL,    'h'   },
    { "pe",     0,   NULL,    'p'   },
    { "min_tlen",     1,   NULL,    'a'   },
    { "max_tlen",     1,   NULL,    'b'   },
    { "group",     1,   NULL,    'g'   },
    { "sw",     0,   NULL,    'e'   },
    { "max_locate",     1,   NULL,    'm'   },
    { "max_seed",     1,   NULL,    's'   },
    { "read_length",     1,   NULL,    'l'   },
    { "overlap",     1,   NULL,    'r'   },
    { "xa_cigar",     0,   NULL,    'c'   },
    { "md",     0,   NULL,    'd'   },
    { "ref",     0,   NULL,    'v'   },
    { "mismatch",     1,   NULL,    'M'   },
    { "gapop",     1,   NULL,    'O'   },
    { "gapex",     1,   NULL,    'E'   },
    { "extend",     1,   NULL,    'X'   },
    { 0,     0,   0,    0   }
};
int opt_parse(int argc, char *argv[], opt_t* opt){

    int i;
    fprintf(stderr, "[opt_parse]:parse opt!\n"); 

    kstring_t *s = opt->cmd;

    for(i = 0; i < argc-1; ++i){
	    ksprintf(s, "%s ", argv[i]);
    }
    ksprintf(s, "%s", argv[i]);

    int c; int option_index=0;
    while((c = getopt_long(argc, argv, short_options, long_options, &option_index))>=0)
    {
        switch(c){
            case 't':
                opt->n_threads = atoi(optarg);
                break;
            case 'n':
                opt->n_diff = atoi(optarg);
                break;
            case 'h':
                return usage();
                break;
            case 'p':
                opt->se = 0;
                break;
            case 'a'://paired end only
                //"min_tlen:max_tlen"
                opt->min_tlen = atoi(optarg);
                break;
            case 'b':
                opt->max_tlen = atoi(optarg);
                break;
            case 'g'://rg
                opt->rg_id = strdup(optarg);
                break;
            case 'e':
                opt->use_sw_extend = 1; 
                break;
            case 's':
                opt->max_seed = atoi(optarg); 
                break;
            case 'm':
                opt->max_locate = atoi(optarg); 
                break;
 
            case 'l':
                opt->l_read = atoi(optarg);
                break;
            case 'c':
                opt->print_xa_cigar = 1; 
                break;
            case 'd':
                opt->print_nm_md = 1;
                break;
            case 'v':
                opt->ref = 1;
                break;

            case 'r':
                opt->l_overlap = atoi(optarg); 
                break;
            case 'M':
                opt->mismatch_penalty = atoi(optarg);
                break;
            case 'O':
                opt->gapop_penalty = atoi(optarg);
                break;
            case 'E':
                opt->gapext_penalty = atoi(optarg);
                break;
            case 'X':
                opt->extend_algo = atoi(optarg);
                break;
            default:
                fprintf(stderr, "[ERROR]: no arg %c\n", c); 
                return EXIT_FAILURE;
              
        }
    
    }
    if(optind >= argc){
        fprintf(stderr, "[opt_parse]: index prefix and read file can't be omited!\n"); 
        exit(1); 
    }
    opt->fn_index = argv[optind++];
    opt->fn_read1 = argv[optind++];
    if(opt->se != 1) opt->fn_read2 = argv[optind++];   
    char fn[1024];
    strncpy(fn, opt->fn_index, 1024); strcat(fn, ".R.seedLen");
    FILE *fp = fopen(fn, "rb");
    if(NULL == fp) {
        fprintf(stderr, "[Error]: Seed length can't be parsed!\n");
        exit(1);
    }    
    fread(&opt->l_seed, sizeof(int), 1, fp);
    if(opt->l_overlap <= 0) opt->l_overlap = opt->l_seed;
    fclose(fp);

    return EXIT_SUCCESS;  

}
void opt_log(const opt_t *opt){
    
    fprintf(stderr, "\n**********************************************************************************\n");
    fprintf(stderr, "Options information!\n");
    fprintf(stderr, "Comand line:       %s\n", opt->cmd->s);
    fprintf(stderr, "Thread num:        %d\n", opt->n_threads);
    fprintf(stderr, "Max diff :         %d\n", opt->n_mismatch);
    if(opt->rg_id != NULL) 
    fprintf(stderr, "Rg_id:             %s\n", opt->rg_id);
    fprintf(stderr, "Index prefix:      %s\n", opt->fn_index);
    fprintf(stderr, "Kmer length:       %d\n", opt->l_seed);
    if(opt->se) {
    fprintf(stderr, "SE/PE:             Single end.\n");
    fprintf(stderr, "Read file name:    %s\n", opt->fn_read1);
    
    } else{
    fprintf(stderr, "SE/PE:             Paired end.\n");
    fprintf(stderr, "Read file name1:   %s\n", opt->fn_read1);
    fprintf(stderr, "Read file name2:   %s\n", opt->fn_read2);
    fprintf(stderr, "Insert size:       [%d, %d]\n", opt->min_tlen, opt->max_tlen);

    }

    fprintf(stderr, "Max locate number: [%d]\n", opt->max_locate);
    fprintf(stderr, "\n**********************************************************************************\n");
}



int aln_main(int argc, char *argv[])
{

    //int i;
    int ret;    

    if(argc < 2){
        usage();
        exit(1);
    }
    fprintf(stderr, "\n************Parse opt*********************************\n");
    opt_t *opt = opt_init();
     
    if(opt_parse( argc, argv, opt) == EXIT_FAILURE){
        fprintf(stderr, "[%s]: opt parse fail, please make sure arguments is right!\n", __func__);
        return EXIT_FAILURE;
    }
    opt_log(opt);
    
    fprintf(stderr, "\n************Parse opt end*****************************\n");
    fprintf(stderr, "\n************Begin to align!***************************\n");
   



    if(opt->se) { ret = alnse_core(opt);} 
    else { ret = alnpe_core(opt);}

    fprintf(stderr, "\n************Finish to align!**************************\n");
    
    opt_destroy(opt);
    return ret;
}

int main ( int argc, char *argv[] )
{
     return aln_main(argc, argv);
}
/* ----------  end of function main  ---------- */


