/* regge_tab.c -- 4D Regge lattice (the model and numbering of regge_r2_par.c) with every neighbor found by formula (rk2.h):
 * the state is the squared edge lengths plus three numbers per triangle, and the only tables are two small ones built below.
 *   kd  : CPU, double precision        kf : CPU, single precision        gpu : Metal, single precision (-DUSE_METAL)
 *   check L k lambda a eps seed cfg|-      (needs -DWITH_ORIG: every formula against the original engine's tables, and chains
 *                                           against the table-driven update of rk.h)
 *   run KERNEL L k lambda a eps nsweeps seed out cfg|-
 * Environment: DCHK=n  every n sweeps (and at the end) recompute in double precision on the CPU and report the largest difference
 *              SAVE=file  write the final configuration (int L, double eps, double q[NE])
 * Build (Linux): gcc -O2 [-DWITH_ORIG] -o regge_tab regge_tab.c -lm -lpthread
 * Build (Mac):   clang -x objective-c -fobjc-arc -O2 -DUSE_METAL -o regge_tab regge_tab.c -framework Metal -framework Foundation */
#ifdef WITH_ORIG
#define main regge_orig_main
#include "regge_r2_par.c"
#undef main
#else
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
static int L, N0, NE, NT, NS;
static int PIDX[5][5], PA[10], PB[10], TRI_OF_PAIR[10][3], tricode[16][16], perms[24][4];
static double eps = 0.5, lam = 1.0, kk = 0.0, aa = 0.0;
static uint64_t rs[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next_u64(void) {
    const uint64_t result = rotl(rs[1] * 5, 7) * 9; const uint64_t t = rs[1] << 17;
    rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45);
    return result;
}
static void seed_rng(uint64_t seed) {
    uint64_t z = seed;
    for (int i = 0; i < 4; i++) { z += 0x9e3779b97f4a7c15ULL; uint64_t x = z; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; rs[i] = x ^ (x >> 31); }
}
static void init_tables(void) {
    int k = 0;
    for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) { PIDX[i][j] = PIDX[j][i] = k; PA[k] = i; PB[k] = j; k++; }
    for (int k2 = 0; k2 < 10; k2++) { int n = 0; for (int v = 0; v < 5; v++) if (v != PA[k2] && v != PB[k2]) TRI_OF_PAIR[k2][n++] = v; }
    int c = 0;
    for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) tricode[a][b] = -1;
    for (int a = 1; a < 16; a++) for (int b = 1; b < 16; b++) if ((a & b) == 0) tricode[a][b] = c++;
    int n = 0;
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) for (int cc = 0; cc < 4; cc++) for (int d = 0; d < 4; d++)
        if (a != b && a != cc && a != d && b != cc && b != d && cc != d) { perms[n][0] = a; perms[n][1] = b; perms[n][2] = cc; perms[n][3] = d; n++; }
}
static double wall_sec(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
#endif

#define rk_uint unsigned int
#define rk_uchar unsigned char
#define DEV
#define THR
#define KTAB static const
#define real double
#define RC(x) x
#define RKN(n) n##_td
#define RK_SQRT sqrt
#define RK_ACOS acos
#define RK_EXP exp
#define RK_TOL 1e-12
#include "rk2.h"
#ifdef WITH_ORIG
#undef RKN
#define RKN(n) n##_kd
#include "rk.h"
#endif
#undef real
#undef RC
#undef RKN
#undef RK_SQRT
#undef RK_ACOS
#undef RK_EXP
#undef RK_TOL
#define real float
#define RC(x) x##f
#define RKN(n) n##_tf
#define RK_SQRT sqrtf
#define RK_ACOS acosf
#define RK_EXP expf
#define RK_TOL 1e-6f
#include "rk2.h"
#undef real
#undef RC
#undef RKN

#ifdef WITH_ORIG
/* the table-driven update of rk.h on the original engine's arrays, for the comparison in "check" */
static int nmax = 0; static double *sV, *sth, *sA; static int *sstat;
static void sweep_kd(unsigned int k0, unsigned int k1) {
    for (int c = 0; c < ncol; c++) { const int b0 = col_off[c], n = col_off[c + 1] - b0;
        for (int i = 0; i < n; i++) update_edge_kd(col_list[b0 + i], i, eps, lam, kk, aa, k0, k1, max_eh, q, V, th, thsum, area, Vh, dacc, dvh,
                                                   es_off, es_list, sedge, shinge, eh_off, eh_list, eh_cont, tedge, sV, sth, sA, sstat); }
}
static double take_acc(int *st) {
    long a = 0, t = 0;
    for (int i = 0; i < nmax; i++) { a += st[4*i]; t += st[4*i + 3]; st[4*i] = st[4*i + 1] = st[4*i + 2] = st[4*i + 3] = 0; }
    return t ? (double)a / (double)t : 0.0;
}
#endif

/* ---------------- the two small tables ---------------- */
static unsigned char TAB[RK2_TAB_BYTES], TAB2[RK2_R_BYTES];
static int mch[24][5];
static int code_of(int base_mask_neg, int plus_mask) {        /* neighbor code of the offset -bits(base) + bits(plus) */
    int code = 0, w = 1;
    for (int i = 0; i < 4; i++) { int o = ((plus_mask >> i) & 1) - ((base_mask_neg >> i) & 1); code += (o + 1) * w; w *= 3; }
    return code;
}
static void build_small_tables(void) {
    for (int p = 0; p < 24; p++) { mch[p][0] = 0; for (int a = 1; a <= 4; a++) mch[p][a] = mch[p][a-1] | (1 << perms[p][a-1]); }
    memset(TAB, 0, sizeof TAB); memset(TAB2, 0, sizeof TAB2);
    for (int t = 0; t < 15; t++) {
        const int M = t + 1; int ns = 0, nh = 0, hl[RK2_MAXH][2];
        for (int p = 0; p < 24; p++) for (int k = 0; k < 10; k++) {
            const int a = PA[k], b = PB[k];
            if ((mch[p][b] ^ mch[p][a]) != M) continue;
            const int base = mch[p][a], i = ns++;                       /* the simplex's base site is v - bits(base) */
            if (i >= 24) { fprintf(stderr, "table: more than 24 simplices around an edge\n"); exit(1); }
            unsigned char *es = TAB + RK2_OFF_ES + (t * 24 + i) * 3; es[0] = (unsigned char)code_of(base, 0); es[1] = (unsigned char)p; es[2] = (unsigned char)k;
            for (int k2 = 0; k2 < 10; k2++) {
                const int a2 = PA[k2], b2 = PB[k2];
                TAB[RK2_OFF_QE + (t * 24 + i) * 20 + 2*k2] = (unsigned char)code_of(base, mch[p][a2]);
                TAB[RK2_OFF_QE + (t * 24 + i) * 20 + 2*k2 + 1] = (unsigned char)((mch[p][b2] ^ mch[p][a2]) - 1);
            }
            for (int k2 = 0; k2 < 10; k2++) {
                const int a3 = TRI_OF_PAIR[k2][0], b3 = TRI_OF_PAIR[k2][1], c3 = TRI_OF_PAIR[k2][2];
                const int S1 = mch[p][b3] ^ mch[p][a3], S2 = mch[p][c3] ^ mch[p][b3];
                const int oc = code_of(base, mch[p][a3]), tc = tricode[S1][S2];
                int j = -1;
                for (int x = 0; x < nh; x++) if (hl[x][0] == oc && hl[x][1] == tc) { j = x; break; }
                if (j < 0) {
                    if (nh >= RK2_MAXH) { fprintf(stderr, "table: more than %d hinges around an edge\n", RK2_MAXH); exit(1); }
                    j = nh++; hl[j][0] = oc; hl[j][1] = tc;
                    unsigned char *te = TAB + RK2_OFF_TE + (t * RK2_MAXH + j) * 6;
                    te[0] = (unsigned char)oc; te[1] = (unsigned char)(S1 - 1);
                    te[2] = (unsigned char)oc; te[3] = (unsigned char)((S1 | S2) - 1);
                    te[4] = (unsigned char)code_of(base, mch[p][a3] | S1); te[5] = (unsigned char)(S2 - 1);     /* a3's site plus the steps S1 (disjoint bits) */
                    int cont = 0; for (int x = 0; x < 3; x++) if (te[2*x] == 40 && te[2*x + 1] == t) cont = 1;
                    unsigned char *eh = TAB + RK2_OFF_EH + (t * RK2_MAXH + j) * 3; eh[0] = (unsigned char)oc; eh[1] = (unsigned char)tc; eh[2] = (unsigned char)cont;
                }
                TAB[RK2_OFF_HJ + (t * 24 + i) * 10 + k2] = (unsigned char)j;
            }
        }
        TAB[RK2_OFF_NS + t] = (unsigned char)ns; TAB[RK2_OFF_NH + t] = (unsigned char)nh;
    }
    for (int p = 0; p < 24; p++) for (int k = 0; k < 10; k++) {
        const int a3 = TRI_OF_PAIR[k][0], b3 = TRI_OF_PAIR[k][1], c3 = TRI_OF_PAIR[k][2];
        const int S1 = mch[p][b3] ^ mch[p][a3], S2 = mch[p][c3] ^ mch[p][b3], tc = tricode[S1][S2], base = mch[p][a3];
        const int i = TAB2[RK2_R_NHT + tc]++;
        if (i >= RK2_MAXHT) { fprintf(stderr, "table: more than %d simplices around a triangle\n", RK2_MAXHT); exit(1); }
        unsigned char *en = TAB2 + RK2_R_ENT + (tc * RK2_MAXHT + i) * 3; en[0] = (unsigned char)code_of(base, 0); en[1] = (unsigned char)p; en[2] = (unsigned char)k;
        for (int k2 = 0; k2 < 10; k2++) {
            TAB2[RK2_R_QE + (tc * RK2_MAXHT + i) * 20 + 2*k2] = (unsigned char)code_of(base, mch[p][PA[k2]]);
            TAB2[RK2_R_QE + (tc * RK2_MAXHT + i) * 20 + 2*k2 + 1] = (unsigned char)((mch[p][PB[k2]] ^ mch[p][PA[k2]]) - 1);
        }
        unsigned char *te = TAB2 + RK2_R_TE + tc * 6;
        te[0] = 40; te[1] = (unsigned char)(S1 - 1); te[2] = 40; te[3] = (unsigned char)((S1 | S2) - 1); te[4] = (unsigned char)code_of(0, S1); te[5] = (unsigned char)(S2 - 1);
    }
}

/* ---------------- state ---------------- */
static double *Dq, *Dts, *Dar, *Dvh;                 /* double state */
static float *Fq, *Fts, *Far, *Fvh;                  /* single state */
static int *stat4, colp2, nper, ncol2, nmax2;
static void colors2(void) {
    colp2 = 0; for (int p = 3; p <= L; p++) if (L % p == 0) { colp2 = p; break; }
    nper = L / colp2; ncol2 = 15 * colp2 * colp2 * colp2 * colp2; nmax2 = nper * nper * nper * nper;
    stat4 = calloc((size_t)nmax2 * 4, sizeof(int));
}
static inline int site_of(int c, int slot, int *t) {     /* color c, slot within the color -> site and edge type */
    *t = c % 15; int r = c / 15, g = slot, v = 0, w = 1;
    for (int i = 0; i < 4; i++) { v += (r % colp2 + colp2 * (g % nper)) * w; r /= colp2; g /= nper; w *= L; }
    return v;
}
static void sweep_td(unsigned k0, unsigned k1) {
    double dI;
    for (int c = 0; c < ncol2; c++) for (int s = 0; s < nmax2; s++) { int t; const int v = site_of(c, s, &t);
        update_edge2_td(v, t, L, eps, lam, kk, aa, k0, k1, Dq, Dts, Dar, Dvh, TAB, stat4 + 4*s, &dI); }
}
static void sweep_tf(unsigned k0, unsigned k1) {
    float dI;
    for (int c = 0; c < ncol2; c++) for (int s = 0; s < nmax2; s++) { int t; const int v = site_of(c, s, &t);
        update_edge2_tf(v, t, L, (float)eps, (float)lam, (float)kk, (float)aa, k0, k1, Fq, Fts, Far, Fvh, TAB, stat4 + 4*s, &dI); }
}
static long recompute_td(void) { long bad = 0; for (int h = 0; h < NT; h++) bad += recompute_hinge_td(h / 50, h % 50, L, Dq, Dts, Dar, Dvh, TAB2); return bad; }
static long recompute_tf(void) { long bad = 0; for (int h = 0; h < NT; h++) bad += recompute_hinge_tf(h / 50, h % 50, L, Fq, Fts, Far, Fvh, TAB2); return bad; }
static double take_acc2(void) {
    long a = 0, t = 0;
    for (int i = 0; i < nmax2; i++) { a += stat4[4*i]; t += stat4[4*i + 3]; stat4[4*i] = stat4[4*i + 1] = stat4[4*i + 2] = stat4[4*i + 3] = 0; }
    return t ? (double)a / (double)t : 0.0;
}
typedef struct { double sv, sda, sr2, sq; } obs_t;
static obs_t observe_d(void) { obs_t o = {0, 0, 0, 0};
    for (int h = 0; h < NT; h++) { const double d = 6.283185307179586476925286766559 - Dts[h]; o.sv += Dvh[h]; o.sda += d * Dar[h]; if (Dvh[h] > 0) o.sr2 += d * d * Dar[h] * Dar[h] / Dvh[h]; }
    for (int e = 0; e < NE; e++) o.sq += Dq[e]; return o; }
static obs_t observe_f(void) { obs_t o = {0, 0, 0, 0};
    for (int h = 0; h < NT; h++) { const double d = 6.283185307179586476925286766559 - (double)Fts[h], A = Far[h], vh = Fvh[h]; o.sv += vh; o.sda += d * A; if (vh > 0) o.sr2 += d * d * A * A / vh; }
    for (int e = 0; e < NE; e++) o.sq += Fq[e]; return o; }

#ifdef USE_METAL
#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
static id<MTLDevice> mdev; static id<MTLCommandQueue> mqueue; static id<MTLComputePipelineState> psweep, precomp;
static id<MTLBuffer> Bq, Bts, Bar, Bvh, Btab, Btab2, Bstat, Bbad;
typedef struct { float eps, lam, kk, aa; unsigned int k0, k1; int L, t, rbase, colp, nper, count; } mparams_t;
static id<MTLBuffer> mbuf(size_t bytes, const void *src) {
    id<MTLBuffer> b = [mdev newBufferWithLength:(bytes ? bytes : 16) options:MTLResourceStorageModeShared];
    if (!b) { fprintf(stderr, "Metal: buffer of %zu bytes not allocated\n", bytes); exit(1); }
    if (src) memcpy([b contents], src, bytes); else memset([b contents], 0, bytes);
    return b;
}
static void metal_setup(const char *dir) {
    mdev = MTLCreateSystemDefaultDevice();
    if (!mdev) { fprintf(stderr, "Metal: no device\n"); exit(1); }
    printf("# Metal device: %s, memory limit %.1f GB\n", [[mdev name] UTF8String], (double)[mdev recommendedMaxWorkingSetSize] / 1e9);
    mqueue = [mdev newCommandQueue];
    char path[1024]; snprintf(path, sizeof path, "%s/rk2.h", dir);
    NSError *err = nil;
    NSString *rk = [NSString stringWithContentsOfFile:[NSString stringWithUTF8String:path] encoding:NSUTF8StringEncoding error:&err];
    if (!rk) { fprintf(stderr, "Metal: cannot read %s\n", path); exit(1); }
    NSString *pre = @"#include <metal_stdlib>\nusing namespace metal;\n#define real float\n#define RC(x) x##f\n#define RKN(n) n\n#define RK_SQRT sqrt\n#define RK_ACOS acos\n"
                     "#define RK_EXP exp\n#define RK_TOL 1e-6f\n#define DEV device\n#define THR thread\n#define KTAB constant\n#define rk_uint uint\n#define rk_uchar uchar\n";
    NSString *post = @"\nstruct Params { float eps, lam, kk, aa; uint k0, k1; int L, t, rbase, colp, nper, count; };\n"
        "kernel void sweep_color2(device float *q [[buffer(0)]], device float *thsum [[buffer(1)]], device float *area [[buffer(2)]], device float *Vh [[buffer(3)]],\n"
        "  device const uchar *tab [[buffer(4)]], device int *stat [[buffer(5)]], constant Params &P [[buffer(6)]], uint gid [[thread_position_in_grid]]) {\n"
        "  if ((int)gid >= P.count) return;\n"
        "  int g = (int)gid, v = P.rbase, w = P.colp;\n"
        "  for (int i = 0; i < 4; i++) { v += (g % P.nper) * w; g /= P.nper; w *= P.L; }\n"
        "  float dI;\n"
        "  update_edge2(v, P.t, P.L, P.eps, P.lam, P.kk, P.aa, P.k0, P.k1, q, thsum, area, Vh, tab, stat + 4 * (int)gid, &dI);\n}\n"
        "kernel void recompute2(device const float *q [[buffer(0)]], device float *thsum [[buffer(1)]], device float *area [[buffer(2)]], device float *Vh [[buffer(3)]],\n"
        "  device const uchar *tab2 [[buffer(4)]], device int *bad [[buffer(5)]], constant Params &P [[buffer(6)]], uint gid [[thread_position_in_grid]]) {\n"
        "  if ((int)gid >= P.count) return;\n"
        "  int b = recompute_hinge((int)gid / 50, (int)gid % 50, P.L, q, thsum, area, Vh, tab2);\n"
        "  if (b > 0) bad[0] = 1;\n}\n";
    NSString *src = [[pre stringByAppendingString:rk] stringByAppendingString:post];
    MTLCompileOptions *opt = [MTLCompileOptions new];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    opt.fastMathEnabled = NO;
#pragma clang diagnostic pop
    id<MTLLibrary> lib = [mdev newLibraryWithSource:src options:opt error:&err];
    if (!lib) { fprintf(stderr, "Metal: shader did not compile:\n%s\n", [[err localizedDescription] UTF8String]); exit(1); }
    psweep = [mdev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"sweep_color2"] error:&err];
    precomp = [mdev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"recompute2"] error:&err];
    if (!psweep || !precomp) { fprintf(stderr, "Metal: pipeline failed: %s\n", [[err localizedDescription] UTF8String]); exit(1); }
    Bq = mbuf(sizeof(float) * NE, NULL); Bts = mbuf(sizeof(float) * NT, NULL); Bar = mbuf(sizeof(float) * NT, NULL); Bvh = mbuf(sizeof(float) * NT, NULL);
    Btab = mbuf(RK2_TAB_BYTES, TAB); Btab2 = mbuf(RK2_R_BYTES, TAB2); Bstat = mbuf(sizeof(int) * nmax2 * 4, NULL); Bbad = mbuf(16, NULL);
    Fq = [Bq contents]; Fts = [Bts contents]; Far = [Bar contents]; Fvh = [Bvh contents]; free(stat4); stat4 = [Bstat contents];
    printf("# Metal: memory on the chip %.2f GB\n", (4.0 * NE + 12.0 * NT) / 1e9);
}
static void sweep_gpu(unsigned k0, unsigned k1) {
    const int per = 256;
    for (int c0 = 0; c0 < ncol2; c0 += per) {
        @autoreleasepool {
            id<MTLCommandBuffer> cb = [mqueue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setComputePipelineState:psweep];
            [enc setBuffer:Bq offset:0 atIndex:0]; [enc setBuffer:Bts offset:0 atIndex:1]; [enc setBuffer:Bar offset:0 atIndex:2]; [enc setBuffer:Bvh offset:0 atIndex:3];
            [enc setBuffer:Btab offset:0 atIndex:4]; [enc setBuffer:Bstat offset:0 atIndex:5];
            for (int c = c0; c < ncol2 && c < c0 + per; c++) {
                int r = c / 15, rbase = 0, w = 1; for (int i = 0; i < 4; i++) { rbase += (r % colp2) * w; r /= colp2; w *= L; }
                mparams_t P = { (float)eps, (float)lam, (float)kk, (float)aa, k0, k1, L, c % 15, rbase, colp2, nper, nmax2 };
                [enc setBytes:&P length:sizeof P atIndex:6];
                [enc dispatchThreads:MTLSizeMake(nmax2, 1, 1) threadsPerThreadgroup:MTLSizeMake(psweep.threadExecutionWidth, 1, 1)];
            }
            [enc endEncoding]; [cb commit]; [cb waitUntilCompleted];
            if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "Metal: sweep failed: %s\n", [[cb.error localizedDescription] UTF8String]); exit(1); }
        }
    }
}
static long recompute_gpu(void) {
    @autoreleasepool {
        id<MTLCommandBuffer> cb = [mqueue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:precomp];
        [enc setBuffer:Bq offset:0 atIndex:0]; [enc setBuffer:Bts offset:0 atIndex:1]; [enc setBuffer:Bar offset:0 atIndex:2]; [enc setBuffer:Bvh offset:0 atIndex:3];
        [enc setBuffer:Btab2 offset:0 atIndex:4]; [enc setBuffer:Bbad offset:0 atIndex:5];
        mparams_t P = { 0, 0, 0, 0, 0, 0, L, 0, 0, colp2, nper, NT };
        [enc setBytes:&P length:sizeof P atIndex:6];
        [enc dispatchThreads:MTLSizeMake(NT, 1, 1) threadsPerThreadgroup:MTLSizeMake(precomp.threadExecutionWidth, 1, 1)];
        [enc endEncoding]; [cb commit]; [cb waitUntilCompleted];
        if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "Metal: recompute failed: %s\n", [[cb.error localizedDescription] UTF8String]); exit(1); }
    }
    int *b = [Bbad contents]; long bad = b[0]; b[0] = 0; return bad;
}
#endif

static void setup2(int Lin, const char *cfg, double eps0, int need_double, int need_float_cpu) {
    L = Lin;
#ifndef WITH_ORIG
    init_tables();
#endif
    N0 = L * L * L * L; NE = 15 * N0; NT = 50 * N0; NS = 24 * N0;
    build_small_tables(); colors2();
    Dq = malloc(sizeof(double) * NE);
    for (int e = 0; e < NE; e++) Dq[e] = (double)__builtin_popcount(e % 15 + 1);       /* flat start */
    if (eps0 > 0) eps = eps0;
    if (cfg && strcmp(cfg, "-")) {
        FILE *fl = fopen(cfg, "rb"); int Lf = 0; double ef;
        if (!fl || fread(&Lf, sizeof(int), 1, fl) != 1 || Lf != L || fread(&ef, sizeof(double), 1, fl) != 1 || fread(Dq, sizeof(double), NE, fl) != (size_t)NE) { fprintf(stderr, "bad config file\n"); exit(1); }
        fclose(fl); if (eps0 <= 0) eps = ef;
    }
    if (need_double) { Dts = malloc(sizeof(double) * NT); Dar = malloc(sizeof(double) * NT); Dvh = malloc(sizeof(double) * NT); }
    if (need_float_cpu) { Fq = malloc(sizeof(float) * NE); Fts = malloc(sizeof(float) * NT); Far = malloc(sizeof(float) * NT); Fvh = malloc(sizeof(float) * NT); }
    printf("# L=%d NE=%d NT=%d colors=%d (p=%d) edges per color=%d; simplices per edge type:", L, NE, NT, ncol2, colp2, nmax2);
    for (int t = 0; t < 15; t++) printf(" %d", TAB[RK2_OFF_NS + t]);
    printf("; hinges per edge type:"); for (int t = 0; t < 15; t++) printf(" %d", TAB[RK2_OFF_NH + t]);
    printf("\n# state in single precision: %.3f GB (4 bytes x (edges + 3 triangles))\n", (4.0 * NE + 12.0 * NT) / 1e9);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: check L k lambda a eps seed cfg|-   |   run kd|kf|gpu L k lambda a eps nsweeps seed out cfg|-\n"); return 1; }
    char dir[1024]; strncpy(dir, argv[0], sizeof dir - 1); dir[sizeof dir - 1] = 0; { char *sl = strrchr(dir, '/'); if (sl) *sl = 0; else strcpy(dir, "."); }
#ifdef WITH_ORIG
    if (!strcmp(argv[1], "check") && argc >= 9) {
        kk = atof(argv[3]); lam = atof(argv[4]); aa = atof(argv[5]);
        L = atoi(argv[2]); init_tables();
        if (!build_lattice(0)) return 1;
        setup2(L, argv[8], atof(argv[6]), 1, 1);
        memcpy(q, Dq, sizeof(double) * NE);
        seed_rng((uint64_t)strtoull(argv[7], NULL, 10));
        /* A. every formula against the original engine's tables */
        long bad = 0; int nb[81];
        for (int e = 0; e < NE; e++) {
            const int v = e / 15, t = e % 15; neigh_td(v, L, nb);
            const int ns = TAB[RK2_OFF_NS + t], nh = TAB[RK2_OFF_NH + t];
            if (ns != es_off[e+1] - es_off[e] || nh != eh_off[e+1] - eh_off[e]) { bad++; continue; }
            for (int i = 0; i < ns; i++) {
                const unsigned char *es = TAB + RK2_OFF_ES + (t * 24 + i) * 3; const int s = nb[es[0]] * 24 + es[1];
                int found = 0; for (int x = es_off[e]; x < es_off[e+1]; x++) if (es_list[x] == s) found = 1;
                if (!found || sedge[s*10 + es[2]] != e) bad++;
                for (int k = 0; k < 10; k++) {
                    if (nb[TAB[RK2_OFF_QE + (t*24 + i)*20 + 2*k]] * 15 + TAB[RK2_OFF_QE + (t*24 + i)*20 + 2*k + 1] != sedge[s*10 + k]) bad++;
                    const int j = TAB[RK2_OFF_HJ + (t*24 + i)*10 + k]; const unsigned char *eh = TAB + RK2_OFF_EH + (t * RK2_MAXH + j) * 3;
                    if (nb[eh[0]] * 50 + eh[1] != shinge[s*10 + k]) bad++;
                }
            }
            for (int j = 0; j < nh; j++) {
                const unsigned char *eh = TAB + RK2_OFF_EH + (t * RK2_MAXH + j) * 3, *te = TAB + RK2_OFF_TE + (t * RK2_MAXH + j) * 6; const int h = nb[eh[0]] * 50 + eh[1];
                int found = -1; for (int x = eh_off[e]; x < eh_off[e+1]; x++) if (eh_list[x] == h) found = x;
                if (found < 0 || eh_cont[found] != eh[2]) { bad++; continue; }
                for (int x = 0; x < 3; x++) if (nb[te[2*x]] * 15 + te[2*x + 1] != tedge[3*h + x]) bad++;
            }
        }
        for (int h = 0; h < NT; h++) {
            const int u = h / 50, c = h % 50; neigh_td(u, L, nb);
            const int n = TAB2[RK2_R_NHT + c]; int cnt = 0;
            for (int i = 0; i < n; i++) { const unsigned char *en = TAB2 + RK2_R_ENT + (c * RK2_MAXHT + i) * 3; const int s = nb[en[0]] * 24 + en[1];
                if (shinge[s*10 + en[2]] != h) bad++;
                for (int k = 0; k < 10; k++) if (nb[TAB2[RK2_R_QE + (c*RK2_MAXHT + i)*20 + 2*k]] * 15 + TAB2[RK2_R_QE + (c*RK2_MAXHT + i)*20 + 2*k + 1] != sedge[s*10 + k]) bad++;
                cnt++; }
            for (int x = 0; x < 3; x++) if (nb[TAB2[RK2_R_TE + c*6 + 2*x]] * 15 + TAB2[RK2_R_TE + c*6 + 2*x + 1] != tedge[3*h + x]) bad++;
            (void)cnt;
        }
        { long tot = 0; for (int c = 0; c < 50; c++) tot += TAB2[RK2_R_NHT + c]; if (tot != 240) bad++; }
        printf("CHECK A  formulas against the engine's tables, every edge, simplex and triangle of the %d^4 lattice: %ld mismatches  (%s)\n", L, bad, bad ? "FAIL" : "PASS");
        /* B. the colors: same edges, same order as the engine's colored sweep */
        setenv("PARFINE", "1", 1); par_on = 1; nthr = 1; par_start();
        long badc = (ncol != ncol2);
        for (int c = 0; c < ncol && !badc; c++) { if (col_off[c+1] - col_off[c] != nmax2) { badc++; break; }
            for (int s = 0; s < nmax2; s++) { int t; const int v = site_of(c, s, &t); if (col_list[col_off[c] + s] != 15 * v + t) badc++; } }
        printf("CHECK B  colors by formula against the engine's color lists: %ld mismatches  (%s)\n", badc, badc ? "FAIL" : "PASS");
        /* C. recomputation by formula against the engine's full recomputation */
        double sv, sda, sa, md; const double aa_keep = aa; if (aa == 0.0) aa = 1e-30;        /* the engine fills hinge volumes only when a is not zero */
        if (recompute_all(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid start\n"); return 1; }
        aa = aa_keep;
        long rb = recompute_td(); double mc = 0;
        for (int h = 0; h < NT; h++) { mc = fmax(mc, fabs(Dts[h] - thsum[h])); mc = fmax(mc, fabs(Dar[h] - area[h])); mc = fmax(mc, fabs(Dvh[h] - Vh[h])); }
        obs_t o = observe_d();
        printf("CHECK C  recomputation by formula against the engine: largest difference %.3e, sums differ by %.1e and %.1e (relative), invalid %ld  (%s)\n",
               mc, fabs(o.sv - sv) / sv, fabs(o.sda - sda) / (fabs(sda) + 1e-300), rb, mc < 1e-10 && !rb ? "PASS" : "FAIL");
        /* D. chains: formula update (double) against the table update of rk.h (double), same random numbers */
        for (int c = 0; c < ncol; c++) if (col_off[c + 1] - col_off[c] > nmax) nmax = col_off[c + 1] - col_off[c];
        const int NSW = 12; sV = malloc(sizeof(double) * nmax * 32); sth = malloc(sizeof(double) * nmax * 320); sA = malloc(sizeof(double) * nmax * max_eh); sstat = calloc(nmax * 4, sizeof(int));
        printf("CHECK D  formula update against table update, double precision, same start and random numbers:\n  sweep  R(table)    R(formula)  acc(t)  acc(f)  edges that differ (of %d)\n", NE);
        long totdiff = 0;
        for (int i = 0; i < NSW; i++) {
            const uint64_t key = next_u64();
            sweep_kd((unsigned)key, (unsigned)(key >> 32)); const double acct = take_acc(sstat);
            if (aa == 0.0) aa = 1e-30; if (recompute_all(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid (table chain)\n"); return 1; } aa = aa_keep;
            sweep_td((unsigned)key, (unsigned)(key >> 32)); const double accf = take_acc2();
            if (recompute_td()) { fprintf(stderr, "invalid (formula chain)\n"); return 1; }
            o = observe_d(); long nd = 0; for (int e = 0; e < NE; e++) if (q[e] != Dq[e]) nd++;
            totdiff += nd;
            printf("  %5d  %+.8f  %+.8f  %.4f  %.4f  %ld\n", i + 1, 2 * sda / sv, 2 * o.sda / o.sv, acct, accf, nd);
        }
        printf("CHECK D  %s\n", totdiff ? "FAIL" : "PASS: the two updates made the same decision on every edge in every sweep");
        /* E. single against double precision with the formula update; single-precision state recomputed in single every sweep */
        double *q0 = malloc(sizeof(double) * NE); memcpy(q0, Dq, sizeof(double) * NE);
        double Rd[64], accd[64]; uint64_t keys[64]; double *snap = malloc(sizeof(double) * NE * NSW);
        for (int i = 0; i < NSW; i++) { keys[i] = next_u64(); sweep_td((unsigned)keys[i], (unsigned)(keys[i] >> 32)); accd[i] = take_acc2(); recompute_td(); o = observe_d(); Rd[i] = 2 * o.sda / o.sv; memcpy(snap + (long)i * NE, Dq, sizeof(double) * NE); }
        for (int e = 0; e < NE; e++) Fq[e] = (float)q0[e];
        recompute_tf();
        printf("CHECK E  single against double precision (formula update):\n  sweep  R(double)   R(single)   acc(d)  acc(s)  edges that took a different decision  largest difference among the rest\n");
        for (int i = 0; i < NSW; i++) {
            sweep_tf((unsigned)keys[i], (unsigned)(keys[i] >> 32)); const double accf = take_acc2(); const long fb = recompute_tf(); obs_t of = observe_f();
            long nd = 0; double mrest = 0; const double *sn = snap + (long)i * NE;
            for (int e = 0; e < NE; e++) { const double d = fabs((double)Fq[e] - sn[e]); if (d > 1e-4 * eps) nd++; else if (d > mrest) mrest = d; }
            printf("  %5d  %+.6f  %+.6f  %.4f  %.4f  %8ld   %.2e%s\n", i + 1, Rd[i], 2 * of.sda / of.sv, accd[i], accf, nd, mrest, fb ? "  INVALID" : "");
        }
        return 0;
    }
#endif
    if (!strcmp(argv[1], "run") && argc >= 12) {
        const char *kern = argv[2]; const int mode = !strcmp(kern, "kd") ? 0 : !strcmp(kern, "kf") ? 1 : 2;
        kk = atof(argv[4]); lam = atof(argv[5]); aa = atof(argv[6]);
        const int dchk = getenv("DCHK") ? atoi(getenv("DCHK")) : 0;
        setup2(atoi(argv[3]), argv[11], atof(argv[7]), mode == 0 || dchk > 0, mode == 1);
        const int nsw = atoi(argv[8]); seed_rng((uint64_t)strtoull(argv[9], NULL, 10));
        if (mode == 2) {
#ifdef USE_METAL
            metal_setup(dir);
#else
            fprintf(stderr, "this build has no Metal\n"); return 1;
#endif
        }
        long bad = 0;
        if (mode == 0) bad = recompute_td();
        else { for (int e = 0; e < NE; e++) Fq[e] = (float)Dq[e];
            if (mode == 1) bad = recompute_tf();
#ifdef USE_METAL
            else bad = recompute_gpu();
#endif
        }
        if (bad) { fprintf(stderr, "invalid start\n"); return 1; }
        FILE *fo = fopen(argv[10], "w");
        fprintf(fo, "# kernel=%s (formula lookups) L=%d k=%.6f lambda=%.6f a=%.6f eps=%.5f seed=%s colors=%d edges per color=%d\n# cols: sweep sumV sumDA R meanq acc sumR2 sweep_sec recompute_sec\n", kern, L, kk, lam, aa, eps, argv[9], ncol2, nmax2);
        double tsw = 0, tre = 0, maxd = 0; obs_t o = {0, 0, 0, 0};
        for (int i = 0; i < nsw; i++) {
            const uint64_t key = next_u64(); const unsigned k0 = (unsigned)key, k1 = (unsigned)(key >> 32); const double w0 = wall_sec(); double w1;
            if (mode == 0) { sweep_td(k0, k1); w1 = wall_sec(); bad = recompute_td(); o = observe_d(); }
            else if (mode == 1) { sweep_tf(k0, k1); w1 = wall_sec(); bad = recompute_tf(); o = observe_f(); }
            else {
#ifdef USE_METAL
                sweep_gpu(k0, k1); w1 = wall_sec(); bad = recompute_gpu(); o = observe_f();
#else
                w1 = w0;
#endif
            }
            const double acc = take_acc2(), w2 = wall_sec(); tsw += w1 - w0; tre += w2 - w1;
            if (bad) { fprintf(stderr, "invalid geometry after sweep %d\n", i); return 2; }
            fprintf(fo, "%d %.10g %.10g %.8f %.8f %.5f %.10g %.4f %.4f\n", i, o.sv, o.sda, 2 * o.sda / o.sv, o.sq / NE, acc, o.sr2, w1 - w0, w2 - w1);
            if ((i + 1) % 50 == 0) fflush(fo);
            if (mode != 0 && dchk > 0 && ((i + 1) % dchk == 0 || i == nsw - 1)) {       /* the same configuration recomputed in double precision */
                for (int e = 0; e < NE; e++) Dq[e] = (double)Fq[e];
                const long db = recompute_td(); const obs_t od = observe_d();
                const double d1 = fabs(od.sv - o.sv) / od.sv, d2 = fabs(2 * od.sda / od.sv - 2 * o.sda / o.sv);
                if (d1 > maxd) maxd = d1;
                fprintf(fo, "# double check after sweep %d: invalid %ld, sumV relative difference %.2e, R difference %.2e\n", i, db, d1, d2);
                printf("double check after sweep %d: invalid %ld, sumV relative difference %.2e, R difference %.2e\n", i, db, d1, d2);
            }
        }
        fclose(fo);
        printf("kernel %s (formula) L=%d: %d sweeps, %.4f s per sweep (update) + %.4f s (recompute and measure); last R=%.6f\n", kern, L, nsw, tsw / nsw, tre / nsw, 2 * o.sda / o.sv);
        if (getenv("SAVE")) { if (mode != 0) for (int e = 0; e < NE; e++) Dq[e] = (double)Fq[e];
            FILE *fs = fopen(getenv("SAVE"), "wb"); fwrite(&L, sizeof(int), 1, fs); fwrite(&eps, sizeof(double), 1, fs); fwrite(Dq, sizeof(double), NE, fs); fclose(fs); }
        return 0;
    }
    fprintf(stderr, "bad arguments\n"); return 1;
}
