// A42 SA74 T2: the open lead at large N. Transitive percolation (or originary percolation) grown to millions of
// elements, read on the largest interval from a root below all elements to the element with the largest past:
//   longest chain exactly (dynamic programming over the direct causes);
//   ordering fraction from sampled ancestor counts (each sample counted exactly by search);
//   midpoint from sampled elements, then climbed along direct causes and effects to the best one.
// Also reads sprinkled intervals of flat d-dimensional spacetime with the same three tests, for the bands at large N.
// usage: tp_big tp N p seed samples        -> one JSON line
//        tp_big op N p seed samples        (originary percolation: every element has at least one direct cause)
//        tp_big mm N p seed samples        (ordering fraction only, for tuning p)
//        tp_big sprinkle d N seed          (sprinkled interval of d-dimensional Minkowski space)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

static uint64_t S[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t nxt(void) { uint64_t r = rotl(S[1] * 5, 7) * 9, t = S[1] << 17; S[2] ^= S[0]; S[3] ^= S[1]; S[1] ^= S[2]; S[0] ^= S[3]; S[2] ^= t; S[3] = rotl(S[3], 45); return r; }
static double unif(void) { return (nxt() >> 11) * 0x1.0p-53; }
static void seedrng(uint64_t x) { for (int i = 0; i < 4; i++) { x += 0x9e3779b97f4a7c15ULL; uint64_t z = x; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; S[i] = z ^ (z >> 31); } }

static double CM[7] = {0, 0, 1.9544, 2.0567, 2.1443, 2.2492, 2.3};
static double mval(double d) { if (d <= 2) return CM[2]; if (d >= 6) return CM[6]; int lo = (int)floor(d); if (lo == 6) return CM[6]; return CM[lo] + (CM[lo + 1] - CM[lo]) * (d - lo); }
static double chain_dim(double n, double L) { double lo = 1, hi = 20; if (mval(lo) * pow(n, 1 / lo) - L <= 0) return 1; if (mval(hi) * pow(n, 1 / hi) - L >= 0) return 20;
  for (int i = 0; i < 100; i++) { double m = (lo + hi) / 2; if (mval(m) * pow(n, 1 / m) - L > 0) lo = m; else hi = m; } return (lo + hi) / 2; }
static double fofd(double d) { return exp(lgamma(d + 1) + lgamma(d / 2) - log(2.0) - lgamma(3 * d / 2)); }
static double mm_dim(double f) { if (f >= 1) return 1; if (f <= fofd(20)) return 20; double lo = 1, hi = 20; for (int i = 0; i < 100; i++) { double m = (lo + hi) / 2; if (fofd(m) > f) lo = m; else hi = m; } return (lo + hi) / 2; }

// graph: elements 1..N; element 0 is the root below all (implicit); parents of j are among 1..j-1
static int N; static int64_t *ps, *cs; static int32_t *pl, *cl; static int64_t E;
static int32_t *stamp; static int32_t curst = 0; static int32_t *stk; static uint8_t *inI;

static void grow(double p, int originary) {
  int64_t cap = (int64_t)(p * (double)N * (double)N / 2 * 1.1) + 1024;
  pl = malloc(sizeof(int32_t) * cap); ps = malloc(sizeof(int64_t) * (N + 2)); E = 0; ps[1] = 0; ps[0] = 0;
  double lq = log(1 - p);
  for (int j = 1; j <= N; j++) {
    ps[j] = E;
    for (int tries = 0;; tries++) {
      int64_t start = E;
      if (j > 1) {
        int64_t i = 1 + (int64_t)floor(log(1 - unif()) / lq);
        while (i <= j - 1) { if (E >= cap) { cap = cap * 3 / 2; pl = realloc(pl, sizeof(int32_t) * cap); } pl[E++] = (int32_t)i; i += 1 + (int64_t)floor(log(1 - unif()) / lq); }
      }
      if (!originary || E > start || j == 1) break;
    }
  }
  ps[N + 1] = E;
  // children
  cs = calloc(N + 2, sizeof(int64_t));
  for (int64_t e = 0; e < E; e++) cs[pl[e] + 1]++;
  for (int j = 1; j <= N + 1; j++) cs[j] += cs[j - 1];
  cl = malloc(sizeof(int32_t) * (E + 1)); int64_t *fill = malloc(sizeof(int64_t) * (N + 2)); memcpy(fill, cs, sizeof(int64_t) * (N + 2));
  for (int j = 1; j <= N; j++) for (int64_t e = ps[j]; e < ps[j + 1]; e++) cl[fill[pl[e]]++] = j;
  free(fill);
  stamp = calloc(N + 2, sizeof(int32_t)); stk = malloc(sizeof(int32_t) * (N + 2)); inI = calloc(N + 2, 1);
}
// number of ancestors of z among 1..N (the root not counted); marks them with the current stamp
static int64_t anc_count(int z) {
  curst++; int64_t top = 0, cnt = 0; stk[top++] = z; stamp[z] = curst;
  while (top) { int u = stk[--top]; for (int64_t e = ps[u]; e < ps[u + 1]; e++) { int v = pl[e]; if (stamp[v] != curst) { stamp[v] = curst; cnt++; stk[top++] = v; } } }
  return cnt;
}
// number of descendants of z inside the interval (marked inI), z not counted
static int64_t desc_in(int z) {
  curst++; int64_t top = 0, cnt = 0; stk[top++] = z; stamp[z] = curst;
  while (top) { int u = stk[--top]; for (int64_t e = cs[u]; e < cs[u + 1]; e++) { int v = cl[e]; if (inI[v] && stamp[v] != curst) { stamp[v] = curst; cnt++; stk[top++] = v; } } }
  return cnt;
}

static void read_tp(const char *mode, int originary, double p, uint64_t seed, int samples) {
  clock_t t0 = clock(); seedrng(seed); grow(p, originary);
  double tg = (double)(clock() - t0) / CLOCKS_PER_SEC;
  if (getenv("DUMP")) { FILE *fp = fopen(getenv("DUMP"), "w"); for (int j = 1; j <= N; j++) { for (int64_t e = ps[j]; e < ps[j + 1]; e++) fprintf(fp, "%d ", pl[e]); fprintf(fp, "\n"); } fclose(fp); }
  // the element with the largest past among the last 40
  int y = N; int64_t best = -1;
  for (int j = N; j > N - 40 && j >= 1; j--) { int64_t a = anc_count(j); if (a > best) { best = a; y = j; } }
  anc_count(y);
  int32_t *mem = malloc(sizeof(int32_t) * (best + 1)); int64_t m = 0;
  for (int j = 1; j < y; j++) if (stamp[j] == curst) { inI[j] = 1; mem[m++] = j; }
  inI[y] = 1;
  int64_t n = m + 2;           // root, the past of y, and y
  // ordering fraction from sampled ancestor counts
  double sum = 0, sum2 = 0;
  for (int s = 0; s < samples; s++) { int z = mem[(int64_t)(unif() * m)]; double a = (double)anc_count(z) + 1; sum += a; sum2 += a * a; }
  double mean = sum / samples, sd = sqrt(fmax(0, sum2 / samples - mean * mean));
  double R = (double)(n - 1) + (double)m * mean;
  double f = R / ((double)n * (double)(n - 1) / 2), fse = (double)m * sd / sqrt((double)samples) / ((double)n * (double)(n - 1) / 2);
  double dmm = mm_dim(f), dlo = mm_dim(f + 2 * fse), dhi = mm_dim(f - 2 * fse);
  if (!strcmp(mode, "mm")) {
    printf("{\"N\": %d, \"p\": %.10g, \"seed\": %llu, \"n\": %lld, \"mm\": %.4f}\n", N, p, (unsigned long long)seed, (long long)n, dmm);
    return;
  }
  // longest chain: dynamic programming over direct causes inside the interval
  int32_t *L = calloc(N + 2, sizeof(int32_t)); int Ly = 0;
  for (int64_t k = 0; k < m; k++) { int j = mem[k]; int b = 0; for (int64_t e = ps[j]; e < ps[j + 1]; e++) if (inI[pl[e]] && L[pl[e]] > b) b = L[pl[e]]; L[j] = b + 1; }
  { int b = 0; for (int64_t e = ps[y]; e < ps[y + 1]; e++) if (inI[pl[e]] && L[pl[e]] > b) b = L[pl[e]]; Ly = b + 1; }
  int64_t chain = (int64_t)Ly + 1;   // plus the root
  double dch = chain_dim((double)n, (double)chain);
  // midpoint: candidates on the balanced slice (longest chain from the root equal to the longest chain to y, within 1),
  // each counted exactly, then climbed along direct causes and effects
  int32_t *H = calloc(N + 2, sizeof(int32_t)); H[y] = 1;
  for (int64_t k = m - 1; k >= 0; k--) { int j = mem[k]; int b = 0; for (int64_t e = cs[j]; e < cs[j + 1]; e++) { int c = cl[e]; if (inI[c] && H[c] > b) b = H[c]; } H[j] = b + 1; }
  int32_t *cand = malloc(sizeof(int32_t) * (m + 1)); int64_t nc = 0;
  int K = 32; int use_sketch = getenv("SKETCH") != NULL;
  if (use_sketch) {
    // bottom-k sketches of ancestor sets and of descendant sets inside the interval, to rank every element by an
    // estimate of min(a, b); the best-ranked are then counted exactly
    int32_t *pos = malloc(sizeof(int32_t) * (N + 2)); for (int64_t k = 0; k < m; k++) pos[mem[k]] = (int32_t)k;
    float *rk = malloc(sizeof(float) * (N + 2)); for (int j = 0; j <= N; j++) rk[j] = (float)unif();
    float *SA = malloc(sizeof(float) * K * m), *SD = malloc(sizeof(float) * K * (m + 1)); uint8_t *la = calloc(m, 1), *ld = calloc(m + 1, 1);
    float buf[64], tmp[64];
    #define MERGE(dst, dl, src, sl, extra) { int i_ = 0, j_ = 0, o_ = 0; float ex_ = (extra); int ue_ = 0; \
      while (o_ < K) { float c_ = 2.0f; int w_ = -1; \
        if (i_ < (dl) && (dst)[i_] < c_) { c_ = (dst)[i_]; w_ = 0; } if (j_ < (sl) && (src)[j_] < c_) { c_ = (src)[j_]; w_ = 1; } if (!ue_ && ex_ < c_) { c_ = ex_; w_ = 2; } \
        if (w_ < 0) break; if (w_ == 0) i_++; else if (w_ == 1) j_++; else ue_ = 1; \
        if (o_ == 0 || tmp[o_ - 1] != c_) tmp[o_++] = c_; } \
      memcpy((dst), tmp, sizeof(float) * o_); (dl) = o_; }
    for (int64_t k = 0; k < m; k++) { int j = mem[k]; int l = 0;
      for (int64_t e = ps[j]; e < ps[j + 1]; e++) { int p_ = pl[e]; if (!inI[p_]) continue; int64_t q = pos[p_]; int pl_ = la[q]; MERGE(buf, l, &SA[q * K], pl_, rk[p_]); }
      memcpy(&SA[k * K], buf, sizeof(float) * l); la[k] = (uint8_t)l; }
    for (int64_t k = m - 1; k >= 0; k--) { int j = mem[k]; int l = 0;
      for (int64_t e = cs[j]; e < cs[j + 1]; e++) { int c_ = cl[e]; if (!inI[c_]) continue; if (c_ == y) { float ex2 = rk[y]; int z0 = 0; MERGE(buf, l, tmp, z0, ex2); continue; } int64_t q = pos[c_]; int cl_ = ld[q]; MERGE(buf, l, &SD[q * K], cl_, rk[c_]); }
      memcpy(&SD[k * K], buf, sizeof(float) * l); ld[k] = (uint8_t)l; }
    // score every element; keep the best T by estimated min(a, b)
    int T = getenv("TOPT") ? atoi(getenv("TOPT")) : 600; double *sc = malloc(sizeof(double) * m);
    for (int64_t k = 0; k < m; k++) { double a = la[k] < K ? la[k] : (K - 1) / SA[k * K + K - 1]; double b = ld[k] < K ? ld[k] : (K - 1) / SD[k * K + K - 1]; sc[k] = a < b ? a : b; }
    for (int t = 0; t < T && t < m; t++) { int64_t bi = -1; double bv = -1; for (int64_t k = 0; k < m; k++) if (sc[k] > bv) { bv = sc[k]; bi = k; } cand[nc++] = mem[bi]; sc[bi] = -2; }
    free(pos); free(rk); free(SA); free(SD); free(la); free(ld); free(sc);
  } else
  { for (int64_t k = 0; k < m; k++) { int j = mem[k]; if (abs(L[j] - H[j]) <= (getenv("SLICE") ? atoi(getenv("SLICE")) : 1)) cand[nc++] = j; } }
  int64_t bestmin = 0; int bz = mem[0];
  int msamp = getenv("MIDALL") ? (int)nc : samples;
  int ns = nc < msamp ? (int)nc : msamp;
  for (int s = 0; s < ns; s++) {
    int z = (nc <= msamp) ? cand[s] : cand[(int64_t)(unif() * nc)];
    int64_t a = anc_count(z) + 2, b = desc_in(z) + 1; int64_t v = a < b ? a : b;
    if (v > bestmin) { bestmin = v; bz = z; }
  }
  if (getenv("DBG")) fprintf(stderr, "nc %lld ns %d best-before-climb %lld\n", (long long)nc, ns, (long long)bestmin);
  for (int step = 0; step < 400; step++) {
    int moved = 0;
    for (int64_t e = ps[bz]; e < ps[bz + 1] && !moved; e++) { int z = pl[e]; if (!inI[z]) continue; int64_t a = anc_count(z) + 2, b = desc_in(z) + 1, v = a < b ? a : b; if (v > bestmin) { bestmin = v; bz = z; moved = 1; } }
    for (int64_t e = cs[bz]; e < cs[bz + 1] && !moved; e++) { int z = cl[e]; if (!inI[z] || z == y) continue; int64_t a = anc_count(z) + 2, b = desc_in(z) + 1, v = a < b ? a : b; if (v > bestmin) { bestmin = v; bz = z; moved = 1; } }
    if (!moved) break;
  }
  double dmid = log2((double)n / (double)bestmin);
  double tt = (double)(clock() - t0) / CLOCKS_PER_SEC;
  printf("{\"model\": \"%s\", \"N\": %d, \"p\": %.10g, \"pN\": %.4f, \"seed\": %llu, \"edges\": %lld, \"y\": %d, \"n\": %lld, \"chain\": %lld, \"chain_d\": %.4f, "
         "\"f\": %.6f, \"f_se\": %.6f, \"mm\": %.4f, \"mm_lo\": %.4f, \"mm_hi\": %.4f, \"mid_best\": %lld, \"mid\": %.4f, \"slice\": %lld, \"samples\": %d, \"grow_sec\": %.2f, \"sec\": %.2f}\n",
         originary ? "originary percolation" : "transitive percolation", N, p, p * N, (unsigned long long)seed, (long long)E, y, (long long)n, (long long)chain, dch,
         f, fse, dmm, dlo, dhi, (long long)bestmin, dmid, (long long)nc, samples, tg, tt);
}

// sprinkled interval of d-dimensional Minkowski space: all three tests (chain exact, relations exact, midpoint by search)
typedef struct { double t, x[5]; } Pt;
static int cmpt(const void *a, const void *b) { double u = ((const Pt *)a)->t, v = ((const Pt *)b)->t; return (u > v) - (u < v); }
static int rel(const Pt *a, const Pt *b, int d) { double dt = b->t - a->t; if (dt <= 0) return 0; double s = 0; for (int k = 0; k < d - 1; k++) { double q = b->x[k] - a->x[k]; s += q * q; } return dt * dt >= s; }
static void read_sprinkle(int d, int n, uint64_t seed) {
  clock_t t0 = clock(); seedrng(seed);
  Pt *P = malloc(sizeof(Pt) * (n + 2)); int have = 0;
  while (have < n) { Pt q; q.t = unif() - 0.5; double s = 0; for (int k = 0; k < d - 1; k++) { q.x[k] = unif() - 0.5; s += q.x[k] * q.x[k]; } if (sqrt(s) <= 0.5 - fabs(q.t)) P[have++] = q; }
  qsort(P, n, sizeof(Pt), cmpt);
  int M = n + 2; Pt *Q = malloc(sizeof(Pt) * M); memset(&Q[0], 0, sizeof(Pt)); Q[0].t = -0.5; memcpy(&Q[1], P, sizeof(Pt) * n); memset(&Q[M - 1], 0, sizeof(Pt)); Q[M - 1].t = 0.5;
  int32_t *L = malloc(sizeof(int32_t) * M); int32_t *na = calloc(M, sizeof(int32_t)); int32_t *nd = calloc(M, sizeof(int32_t)); double R = 0;
  for (int j = 0; j < M; j++) { int b = 0; for (int i = 0; i < j; i++) if (rel(&Q[i], &Q[j], d)) { R += 1; na[j]++; nd[i]++; if (L[i] > b) b = L[i]; } L[j] = b + 1; }
  double f = R / ((double)M * (M - 1) / 2);
  // midpoint: largest smaller half over all inner points (counts include the point and the tip)
  int64_t best = 0; for (int z = 1; z < M - 1; z++) { int64_t a = na[z] + 1, b = nd[z] + 1, v = a < b ? a : b; if (v > best) best = v; }
  printf("{\"model\": \"sprinkled\", \"d\": %d, \"N\": %d, \"seed\": %llu, \"n\": %d, \"chain\": %d, \"chain_d\": %.4f, \"f\": %.6f, \"mm\": %.4f, \"mid\": %.4f, \"sec\": %.2f}\n",
         d, n, (unsigned long long)seed, M, L[M - 1], chain_dim((double)M, (double)L[M - 1]), f, mm_dim(f), log2((double)M / (double)best), (double)(clock() - t0) / CLOCKS_PER_SEC);
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: see header\n"); return 1; }
  if (!strcmp(argv[1], "sprinkle")) { read_sprinkle(atoi(argv[2]), atoi(argv[3]), strtoull(argv[4], 0, 10)); return 0; }
  N = atoi(argv[2]); double p = atof(argv[3]); uint64_t seed = strtoull(argv[4], 0, 10); int samples = argc > 5 ? atoi(argv[5]) : 600;
  read_tp(argv[1], !strcmp(argv[1], "op"), p, seed, samples);
  return 0;
}
