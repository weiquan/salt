#ifndef MIXREF_H
#define MIXREF_H
#include <stdint.h>
typedef struct{
    uint32_t *seq;
    uint32_t l;
} mixRef_t;

mixRef_t *mixRef_restore(const char *fn);
void mixRef_destroy(mixRef_t *mixRef);
int build_mixRef(const char *fn_fa, const char *fn_hapmap, const char *fn_mixRef);
#endif
