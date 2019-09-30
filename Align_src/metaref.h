#include <stdint.h>
typedef struct{
    uint32_t *seq;
    uint32_t l;
} mixRef_t;

mixRef_t *mixRef_restore(const char *fn);
void mixRef_destroy(mixRef_t *mixRef);
int main_gen_mixRef(const char *fn_fa, const char *fn_snp, const char *fn_mixRef);
