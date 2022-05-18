#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "kthread.h"
#include "yak-priv.h"
#include "bseq.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef struct {
	int64_t c[YAK_N_COUNTS];
	int32_t max;
	uint64_t *s;
} qv_cntbuf_t;

typedef struct {
	int k;
	const yak_qopt_t *opt;
	bseq_file_t *ks;
	const yak_ch_t *ch;
	qv_cntbuf_t *buf;
} qv_shared_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	qv_shared_t *qs;
} qv_step_t;

static void worker_qv(void *_data, long k, int tid)
{
	qv_step_t *data = (qv_step_t*)_data;
	qv_shared_t *qs = data->qs;
	bseq1_t *s = &data->seqs[k];
	qv_cntbuf_t *b = &qs->buf[tid];
	int i, l, tot, non0, shift = 2 * (qs->ch->k - 1);
	uint64_t x[2], mask = (1ULL<<2*qs->ch->k) - 1;

	assert(qs->ch->k < 32);
	if (s->l_seq < qs->opt->min_len) return;
	if (b->max < s->l_seq) {
		b->max = s->l_seq;
		kroundup32(b->max);
		b->s = (uint64_t*)realloc(b->s, b->max * sizeof(uint64_t));
	}
	for (i = l = 0, tot = non0 = 0, x[0] = x[1] = 0; i < s->l_seq; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[i]];
		if (c < 4) {
			x[0] = (x[0] << 2 | c) & mask;
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			if (++l >= qs->k) {
				int t;
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				y = yak_hash64(y, mask);
				t = yak_ch_get(qs->ch, y);
				if (t < 0) t = 0;
				if (t > 0) ++non0;
				b->s[tot++] = (uint64_t)i<<32 | t;
			}
		} else l = 0, x[0] = x[1] = 0;
	}

	if (qs->opt->print_each) {
		double qv = -1.0;
		if (tot > 0) {
			if (non0 > 0) {
				if (tot > non0) {
					qv = log((double)tot / non0) / qs->k;
					qv = -4.3429448190325175 * log(qv);
				} else qv = 99.0;
			} else qv = 0.0;
		}
		printf("SQ\t%s\t%d\t%d\t%d\t%.2f\n", s->name, s->l_seq, tot, non0, qv);
	}

	if (non0 < tot * qs->opt->min_frac) return;
	for (i = 0; i < tot; ++i)
		++b->c[(int32_t)b->s[i]];
}

static void *yak_qv_cb(void *shared, int step, void *_data)
{
	qv_shared_t *qs = (qv_shared_t*)shared;
	if (step == 0) {
		qv_step_t *ret;
		ret = calloc(1, sizeof(qv_step_t));
		ret->seqs = bseq_read(qs->ks, qs->opt->chunk_size, 0, &ret->n_seqs);
		ret->qs = qs;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		int i;
		double rt, eff;
		qv_step_t *data = (qv_step_t*)_data;
		kt_for(qs->opt->n_threads, worker_qv, data, data->n_seqs);
		rt = yak_realtime();
		eff = yak_cputime() / (rt + 1e-6);
		fprintf(stderr, "[M::%s@%.2f*%.2f] processed %d sequences\n", __func__, rt, eff, data->n_seqs);
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			free(s->seq); free(s->qual); free(s->comment); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

void yak_qv(const yak_qopt_t *opt, const char *fn, const yak_ch_t *ch, int64_t *cnt)
{
	qv_shared_t qs;
	int i, j, n_cnt = 1<<YAK_COUNTER_BITS;
	memset(&qs, 0, sizeof(qv_shared_t));
	qs.k = ch->k;
	qs.opt = opt;
	qs.ch = ch;
	qs.buf = calloc(opt->n_threads, sizeof(qv_cntbuf_t));
	qs.ks = bseq_open(fn);
	kt_pipeline(2, yak_qv_cb, &qs, 2);
	bseq_close(qs.ks);
	memset(cnt, 0, n_cnt * sizeof(int64_t));
	for (i = 0; i < opt->n_threads; ++i) {
		for (j = 0; j < n_cnt; ++j)
			cnt[j] += qs.buf[i].c[j];
		free(qs.buf[i].s);
	}
	free(qs.buf);
}

void yak_qopt_init(yak_qopt_t *opt)
{
	memset(opt, 0, sizeof(yak_qopt_t));
	opt->chunk_size = 1000000000;
	opt->n_threads = 4;
	opt->min_frac = 0.5;
	opt->fpr = 0.00004;
}

int yak_qv_solve(const int64_t *hist, const int64_t *cnt, int kmer, double fpr, yak_qstat_t *qs)
{
	extern int gjdn(double *a, double *b, int n, int m);
	const int max_pow = 2;
	int32_t i, j, k, c, n_cnt = YAK_N_COUNTS, max_c, max_cnt, min_c, min_cnt, n_ext;
	double adj_sum, x[YAK_N_COUNTS], y[YAK_N_COUNTS], A[5 * 5], B[5]; // max power: 4
	double *xp;

	memset(qs, 0, sizeof(yak_qstat_t));
	qs->qv = -1.0, qs->err = cnt[0];
	for (c = 0, qs->tot = 0; c < n_cnt; ++c)
		qs->tot += cnt[c], qs->adj_cnt[c] = cnt[c];
	qs->qv_raw = qs->tot > 0 && qs->tot > cnt[0]? -4.3429448190325175 * log(log((double)qs->tot / (qs->tot - cnt[0])) / kmer) : -1.0;

	// find the max and the min
	for (c = 2, max_cnt = 0, max_c = -1; c < n_cnt - 1; ++c)
		if (max_cnt < cnt[c]) max_cnt = cnt[c], max_c = c;
	for (c = 2, min_cnt = max_cnt, min_c = -1; c < max_c; ++c)
		if (min_cnt > cnt[c]) min_cnt = cnt[c], min_c = c;
	qs->cov = (double)cnt[max_c] / hist[max_c];

	// find the upper fpr
	qs->fpr_upper = 1.0;
	for (c = 2; c < max_c; ++c) {
		double e = cnt[c] / (qs->cov * hist[c]);
		if (qs->fpr_upper > e) qs->fpr_upper = e;
	}
	if (fpr > qs->fpr_upper) fpr = qs->fpr_upper * 0.5;

	// find the lower bound (if possible)
	qs->fpr_lower = 0.0;
	if (min_c > 2 && hist[2] > hist[min_c]) {
		double e = (cnt[2] - cnt[min_c]) / (qs->cov * (hist[2] - hist[min_c]));
		if (qs->fpr_lower < e) qs->fpr_lower = e;
	}
	if (fpr < qs->fpr_lower) fpr = qs->fpr_lower;
	if (qs->fpr_lower >= qs->fpr_upper)
		fprintf(stderr, "Warning: the FPR upper bound is smaller than the lower bound. Trust the lower bound.\n");

	// don't compute adjusted qv if the k-mer histogram is not derived from high-coverage data.
	if (max_c <= 4) return -1;
	n_ext = max_c - min_c + 1 < 8? max_c - min_c + 1 : 8;
	if (n_ext < 3) return -1;

	// compute adj_cnt[] in the range of [min_c, max_c)
	for (c = max_c - 1; c >= min_c; --c) {
		double err = (hist[c] - cnt[c] / qs->cov) / (1.0 - fpr);
		qs->adj_cnt[c] = cnt[c] - err * qs->cov * fpr;
		if (qs->adj_cnt[c] < 0.0) qs->adj_cnt[c] = 0.0;
		//if (qs->adj_cnt[c] > qs->adj_cnt[c+1]) qs->adj_cnt[c] = qs->adj_cnt[c+1] * 0.99;
	}

	// fit the tail
	for (k = 0; k < n_ext; ++k) {
		x[k] = min_c + k;
		y[k] = qs->adj_cnt[min_c + k + 1] / qs->adj_cnt[min_c + k];
	}
	xp = (double*)calloc(n_ext * (max_pow * 2 + 1), sizeof(double));
	for (k = 0; k < n_ext; ++k) {
		double t = 1.0;
		for (i = 0; i <= max_pow * 2; ++i)
			xp[i * n_ext + k] = t, t *= x[k];
	}
	for (i = 0; i <= max_pow; ++i) {
		double sum;
		for (j = 0; j <= i; ++j) {
			for (k = 0, sum = 0.0; k < n_ext; ++k)
				sum += xp[(i + j) * n_ext + k];
			A[i * (max_pow + 1) + j] = A[j * (max_pow + 1) + i] = sum;
		}
		for (k = 0, sum = 0.0; k < n_ext; ++k)
			sum += xp[i * n_ext + k] * y[k];
		B[i] = sum;
	}
	gjdn(A, B, max_pow + 1, 1);
	free(xp);

	// extrapolate to the rest
	for (c = min_c - 1; c >= 0; --c) {
		double r = 0.0, t = 1.0;
		for (i = 0; i <= max_pow; ++i)
			r += B[i] * t, t *= c;
		if (r < 1.01) r = 1.01;
		qs->adj_cnt[c] = qs->adj_cnt[c + 1] / r;
	}

	// compute adjusted qv
	for (c = 0, adj_sum = 0.0; c < n_cnt; ++c)
		adj_sum += qs->adj_cnt[c];
	if (adj_sum <= (double)qs->tot) {
		qs->err = qs->tot - adj_sum;
		qs->qv = -4.3429448190325175 * log(log(qs->tot / adj_sum) / kmer);
	} else {
		fprintf(stderr, "WARNING: failed to estimate the calibrated QV\n");
		qs->err = 0;
		qs->qv = qs->qv_raw;
	}
	return 0;
}
