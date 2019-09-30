#ifndef __SAM_H
#define __SAM_H
#include <stdio.h>
#include "query.h"
#include "aln.h"

void aln_samhead(const opt_t *opt, bntseq_t *bntseq);
void aln_samse(index_t *index, query_t *query, const aln_opt_t *opt);
void alnpe_sam(index_t *index, query_t *query, const aln_opt_t *opt);

#endif 
