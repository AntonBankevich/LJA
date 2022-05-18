#ifndef DR_YAK_LIB_H
#define DR_YAK_LIB_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ketopt.h"
#include "bseq.h"
#include "yak-priv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include "kthread.h"
#include "ketopt.h"
#include "yak-priv.h"


#ifdef __cplusplus
extern "C" {
#endif

int lib_triobin(int t, const char* paternal_yak, const char* maternal_yak, const char* contigs);

int lib_count(int32_t bf_shift, int32_t n_thread, const char *out_string, const char *in_string);

#ifdef __cplusplus
}
#endif

#endif //DR_YAK_LIB_H
