#include "yak_lib.h"


static inline int64_t mm_parse_num(const char *str)
{
    double x;
    char *p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g') x *= 1e9;
    else if (*p == 'M' || *p == 'm') x *= 1e6;
    else if (*p == 'K' || *p == 'k') x *= 1e3;
    return (int64_t)(x + .499);
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define CHUNK_SIZE 200000000

typedef struct {
    int c[16];
    int sc[2];
    int nk;
} tb_cnt_t;

typedef struct {
    int max;
    uint32_t *s;
} tb_buf_t;

typedef struct {
    int k, n_threads, print_diff;
    double ratio_thres;
    bseq_file_t *fp;
    const yak_ch_t *ch;
    tb_buf_t *buf;
} tb_shared_t;

typedef struct {
    int n_seq;
    tb_shared_t *aux;
    bseq1_t *seq;
    tb_cnt_t *cnt;
} tb_step_t;

static void tb_worker(void *_data, long k, int tid)
{
    tb_step_t *t = (tb_step_t*)_data;
    tb_shared_t *aux = t->aux;
    bseq1_t *s = &t->seq[k];
    tb_buf_t *b = &aux->buf[tid];
    uint64_t x[4], mask;
    int i, l, shift;
    if (aux->ch->k < 32) {
        mask = (1ULL<<2*aux->ch->k) - 1;
        shift = 2 * (aux->ch->k - 1);
    } else {
        mask = (1ULL<<aux->ch->k) - 1;
        shift = aux->ch->k - 1;
    }
    if (s->l_seq > b->max) {
        b->max = s->l_seq;
        kroundup32(b->max);
        b->s = (uint32_t*)realloc(b->s, b->max * sizeof(uint32_t));
    }
    memset(b->s, 0, s->l_seq * sizeof(uint32_t));
    for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->l_seq; ++i) {
        int flag, c = seq_nt4_table[(uint8_t)s->seq[i]];
        if (c < 4) {
            if (aux->ch->k < 32) {
                x[0] = (x[0] << 2 | c) & mask;
                x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
            } else {
                x[0] = (x[0] << 1 | (c&1))  & mask;
                x[1] = (x[1] << 1 | (c>>1)) & mask;
                x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
                x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
            }
            if (++l >= aux->k) {
                int type = 0, c1, c2;
                uint64_t y;
                ++t->cnt[k].nk;
                if (aux->ch->k < 32)
                    y = yak_hash64(x[0] < x[1]? x[0] : x[1], mask);
                else
                    y = yak_hash_long(x);
                flag = yak_ch_get(aux->ch, y);
                if (flag < 0) flag = 0;
                c1 = flag&3, c2 = flag>>2&3;
                if (c1 == 2 && c2 == 0) type = 1;
                else if (c2 == 2 && c1 == 0) type = 2;
                b->s[i] = type;
                ++t->cnt[k].c[flag];
                if (aux->print_diff && (flag>>2&3) != (flag&3))
                    printf("D\t%s\t%d\t%d\t%d\n", s->name, i, flag&3, flag>>2&3);
            }
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
    }
    for (l = 0, i = 1; i <= s->l_seq; ++i) {
        if (i == s->l_seq || b->s[i] != b->s[l]) {
            if (b->s[l] > 0 && i - l >= aux->k - 4)
                t->cnt[k].sc[b->s[l] - 1] += i - l;
            l = i;
        }
    }
}

static char tb_classify(const int sc[2], const int *c, int k, double ratio_thres)
{
    char type;
    if (sc[0] == 0 && sc[1] == 0) {
        if (c[0<<2|2] == c[2<<2|0]) type = '0';
        else if (c[0<<2|2] >= k - 4 + c[2<<2|0] && (c[2<<2|0] <= 1 || c[0<<2|2] * 0.05 > c[2<<2|0])) type = 'p';
        else if (c[2<<2|0] >= k - 4 + c[0<<2|2] && (c[0<<2|2] <= 1 || c[2<<2|0] * 0.05 > c[0<<2|2])) type = 'm';
        else type = '0';
    } else if (sc[0] > k && sc[1] > k) {
        type = 'a';
    } else if (sc[0] >= k - 4 + sc[1] && sc[0] * 0.05 >= sc[1] && c[0<<2|2] * ratio_thres > c[2<<2|0]) {
        type = 'p';
    } else if (sc[1] >= k - 4 + sc[0] && sc[1] * 0.05 >= sc[0] && c[2<<2|0] * ratio_thres > c[0<<2|2]) {
        type = 'm';
    } else {
        type = 'a';
    }
    return type;
}

static void *tb_pipeline(void *shared, int step, void *_data)
{
    tb_shared_t *aux = (tb_shared_t*)shared;
    if (step == 0) {
        tb_step_t *s;
        s = (tb_step_t*)calloc(1, sizeof(tb_step_t));
        s->seq = bseq_read(aux->fp, CHUNK_SIZE, 0, &s->n_seq);
        s->aux = aux;
        if (s->n_seq) {
            s->cnt = (tb_cnt_t*)calloc(s->n_seq, sizeof(tb_cnt_t));
            fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
            return s;
        } else free(s);
    } else if (step == 1) {
        int i;
        tb_step_t *s = (tb_step_t*)_data;
        kt_for(aux->n_threads, tb_worker, s, s->n_seq);
        for (i = 0; i < s->n_seq; ++i) {
            int *c = s->cnt[i].c;
            char type;
            type = tb_classify(s->cnt[i].sc, c, aux->k, aux->ratio_thres);
            printf("%s\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", s->seq[i].name, type, s->cnt[i].sc[0], s->cnt[i].sc[1],
                   c[0<<2|2], c[2<<2|0], c[0<<2|1], c[1<<2|0], s->cnt[i].nk, c[0]);
            free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
        }
        free(s->seq); free(s->cnt); free(s);
    }
    return 0;
}

int lib_triobin(int t, const char* paternal_yak, const char* maternal_yak, const char* contigs) {
    ketopt_t o = KETOPT_INIT;
    int i, c, min_cnt = 2, mid_cnt = 5;
    yak_ch_t *ch;
    tb_shared_t aux;

    memset(&aux, 0, sizeof(tb_shared_t));
    aux.n_threads = 8, aux.print_diff = 0;
    aux.ratio_thres = 0.33;
    aux.n_threads = t;
    ch = yak_ch_restore_core(0,  paternal_yak, YAK_LOAD_TRIOBIN1, min_cnt, mid_cnt);
    ch = yak_ch_restore_core(ch, maternal_yak, YAK_LOAD_TRIOBIN2, min_cnt, mid_cnt);

    aux.k = ch->k;
    aux.fp = bseq_open(contigs);
    if (aux.fp == 0) {
        fprintf(stderr, "ERROR: fail to open file '%s'\n", contigs);
        exit(1);
    }
    aux.ch = ch;
    aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
    kt_pipeline(2, tb_pipeline, &aux, 2);
    bseq_close(aux.fp);
    yak_ch_destroy(ch);
    for (i = 0; i < aux.n_threads; ++i) free(aux.buf[i].s);
    free(aux.buf);
    return 0;
}

int lib_count(int32_t bf_shift, int32_t n_thread, const char* out_string, const char* in_string ) {
    yak_ch_t *h;
    int c;
//    char *fn_out = 0;
    yak_copt_t opt;
    ketopt_t o = KETOPT_INIT;
    yak_copt_init(&opt);
    opt.bf_shift = bf_shift;
    opt.n_thread = n_thread;
    const char* fn_out = out_string;

    h = yak_count(in_string, &opt, 0);
    if (opt.bf_shift > 0) {
        yak_ch_destroy_bf(h);
        yak_ch_clear(h, opt.n_thread);
        h = yak_count(in_string, &opt, h);
        yak_ch_shrink(h, 2, YAK_MAX_COUNT, opt.n_thread);
        fprintf(stderr, "[M::%s] %ld distinct k-mers after shrinking\n", __func__, (long)h->tot);
    }
    if (fn_out) yak_ch_dump(h, fn_out);
    yak_ch_destroy(h);
    return 0;
}
