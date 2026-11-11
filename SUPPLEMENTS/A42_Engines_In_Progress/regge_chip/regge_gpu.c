/* regge_gpu.c -- the colored sweep of regge_r2_par.c with its single-edge update taken from rk.h, in three builds:
 *   kd  : CPU, double precision (reference)
 *   kf  : CPU, single precision (what the graphics chip computes, run on the CPU)
 *   gpu : Metal, single precision (only when built with -DUSE_METAL on a Mac)
 * After every sweep the squared edge lengths are copied to double precision and volumes, angles and sums are recomputed there
 * (recompute_all of the original engine), so measurements are in double precision in all three builds, and the single-precision
 * copies start each sweep fresh. The lattice, the colors and the measurements are the original engine's, included unchanged.
 *
 *   check L k lambda a eps seed cfg|-            tests of the update itself and single against double precision
 *   run KERNEL L k lambda a eps nsweeps seed out cfg|-    a run with one kernel; writes sweep, sums, R, acceptance, seconds
 * Build (Linux):  gcc -O2 -o regge_gpu regge_gpu.c -lm -lpthread
 * Build (Mac):    clang -x objective-c -fobjc-arc -O2 -DUSE_METAL -o regge_gpu regge_gpu.c -framework Metal -framework Foundation */
#define main regge_orig_main
#include "regge_r2_par.c"
#undef main
#define rk_uint unsigned int
#define rk_uchar unsigned char
#define DEV
#define THR
#define KTAB static const

#define real double
#define RC(x) x
#define RKN(n) n##_kd
#define RK_SQRT sqrt
#define RK_ACOS acos
#define RK_EXP exp
#define RK_TOL 1e-12
#include "rk.h"
#undef real
#undef RC
#undef RKN
#undef RK_SQRT
#undef RK_ACOS
#undef RK_EXP
#undef RK_TOL

#define real float
#define RC(x) x##f
#define RKN(n) n##_kf
#define RK_SQRT sqrtf
#define RK_ACOS acosf
#define RK_EXP expf
#define RK_TOL 1e-6f
#include "rk.h"
#undef real
#undef RC
#undef RKN

static int nmax = 0;                                  /* largest color */
static double *sV, *sth, *sA; static int *sstat;      /* scratch of the double kernel */
typedef struct { float *q, *V, *th, *thsum, *area, *Vh, *dacc, *dvh, *tmpV, *tmpth, *tmpA; int *stat; } fb_t;
static fb_t F;
static long eh_total = 0, es_total = 0;

#ifdef USE_METAL
#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
static id<MTLDevice> mdev; static id<MTLCommandQueue> mqueue; static id<MTLComputePipelineState> mpipe;
static id<MTLBuffer> MB[22];
typedef struct { float eps, lam, kk, aa; unsigned int k0, k1; int max_eh, offset, count, pad; } mparams_t;
static void *mbuf(int idx, size_t bytes, const void *src) {
    MB[idx] = [mdev newBufferWithLength:(bytes ? bytes : 16) options:MTLResourceStorageModeShared];
    if (!MB[idx]) { fprintf(stderr, "Metal: buffer %d of %zu bytes not allocated\n", idx, bytes); exit(1); }
    if (src) memcpy([MB[idx] contents], src, bytes); else memset([MB[idx] contents], 0, bytes);
    return [MB[idx] contents];
}
static void metal_setup(const char *dir) {
    mdev = MTLCreateSystemDefaultDevice();
    if (!mdev) { fprintf(stderr, "Metal: no device\n"); exit(1); }
    printf("# Metal device: %s, memory limit %.1f GB\n", [[mdev name] UTF8String], (double)[mdev recommendedMaxWorkingSetSize] / 1e9);
    mqueue = [mdev newCommandQueue];
    char path[1024]; snprintf(path, sizeof path, "%s/rk.h", dir);
    NSError *err = nil;
    NSString *rk = [NSString stringWithContentsOfFile:[NSString stringWithUTF8String:path] encoding:NSUTF8StringEncoding error:&err];
    if (!rk) { fprintf(stderr, "Metal: cannot read %s\n", path); exit(1); }
    NSString *pre = @"#include <metal_stdlib>\nusing namespace metal;\n#define real float\n#define RC(x) x##f\n#define RKN(n) n\n#define RK_SQRT sqrt\n#define RK_ACOS acos\n"
                     "#define RK_EXP exp\n#define RK_TOL 1e-6f\n#define DEV device\n#define THR thread\n#define KTAB constant\n#define rk_uint uint\n#define rk_uchar uchar\n";
    NSString *post = @"\nstruct Params { float eps, lam, kk, aa; uint k0, k1; int max_eh, offset, count, pad; };\n"
        "kernel void sweep_color(device float *q [[buffer(0)]], device float *V [[buffer(1)]], device float *th [[buffer(2)]], device float *thsum [[buffer(3)]],\n"
        "  device float *area [[buffer(4)]], device float *Vh [[buffer(5)]], device float *dacc [[buffer(6)]], device float *dvh [[buffer(7)]],\n"
        "  device const int *es_off [[buffer(8)]], device const int *es_list [[buffer(9)]], device const int *sedge [[buffer(10)]], device const int *shinge [[buffer(11)]],\n"
        "  device const int *eh_off [[buffer(12)]], device const int *eh_list [[buffer(13)]], device const uchar *eh_cont [[buffer(14)]], device const int *tedge [[buffer(15)]],\n"
        "  device const int *col_list [[buffer(16)]], device float *tmpV [[buffer(17)]], device float *tmpth [[buffer(18)]], device float *tmpA [[buffer(19)]],\n"
        "  device int *stat [[buffer(20)]], constant Params &P [[buffer(21)]], uint gid [[thread_position_in_grid]]) {\n"
        "  if ((int)gid >= P.count) return;\n"
        "  update_edge(col_list[P.offset + (int)gid], (int)gid, P.eps, P.lam, P.kk, P.aa, P.k0, P.k1, P.max_eh, q, V, th, thsum, area, Vh, dacc, dvh,\n"
        "              es_off, es_list, sedge, shinge, eh_off, eh_list, eh_cont, tedge, tmpV, tmpth, tmpA, stat);\n}\n";
    NSString *src = [[pre stringByAppendingString:rk] stringByAppendingString:post];
    MTLCompileOptions *opt = [MTLCompileOptions new];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    opt.fastMathEnabled = NO;                           /* keep IEEE arithmetic: no reordering, no approximate functions */
#pragma clang diagnostic pop
    id<MTLLibrary> lib = [mdev newLibraryWithSource:src options:opt error:&err];
    if (!lib) { fprintf(stderr, "Metal: shader did not compile:\n%s\n", [[err localizedDescription] UTF8String]); exit(1); }
    mpipe = [mdev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"sweep_color"] error:&err];
    if (!mpipe) { fprintf(stderr, "Metal: pipeline failed: %s\n", [[err localizedDescription] UTF8String]); exit(1); }
    F.q = mbuf(0, sizeof(float) * NE, NULL); F.V = mbuf(1, sizeof(float) * NS, NULL); F.th = mbuf(2, sizeof(float) * NS * 10, NULL);
    F.thsum = mbuf(3, sizeof(float) * NT, NULL); F.area = mbuf(4, sizeof(float) * NT, NULL); F.Vh = mbuf(5, sizeof(float) * NT, NULL);
    F.dacc = mbuf(6, sizeof(float) * NT, NULL); F.dvh = mbuf(7, sizeof(float) * NT, NULL);
    /* the tables live once, in memory the chip and the CPU both see: copy, free the CPU copy, point the engine at the shared one */
#define SHARE(idx, ptr, bytes) do { void *n_ = mbuf(idx, bytes, ptr); free(ptr); ptr = n_; } while (0)
    SHARE(8, es_off, sizeof(int) * (NE + 1)); SHARE(9, es_list, sizeof(int) * es_total); SHARE(10, sedge, sizeof(int) * NS * 10); SHARE(11, shinge, sizeof(int) * NS * 10);
    SHARE(12, eh_off, sizeof(int) * (NE + 1)); SHARE(13, eh_list, sizeof(int) * eh_total); SHARE(14, eh_cont, eh_total); SHARE(15, tedge, sizeof(int) * NT * 3);
    SHARE(16, col_list, sizeof(int) * NE);
    F.tmpV = mbuf(17, sizeof(float) * nmax * 32, NULL); F.tmpth = mbuf(18, sizeof(float) * nmax * 320, NULL); F.tmpA = mbuf(19, sizeof(float) * nmax * max_eh, NULL);
    F.stat = mbuf(20, sizeof(int) * nmax * 4, NULL);
    printf("# Metal: thread width %lu, max threads per group %lu\n", (unsigned long)mpipe.threadExecutionWidth, (unsigned long)mpipe.maxTotalThreadsPerThreadgroup);
}
static void sweep_gpu(unsigned int k0, unsigned int k1) {
    const int per = 256;                                  /* colors per command buffer */
    for (int c0 = 0; c0 < ncol; c0 += per) {
        @autoreleasepool {
            id<MTLCommandBuffer> cb = [mqueue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setComputePipelineState:mpipe];
            for (int i = 0; i < 21; i++) [enc setBuffer:MB[i] offset:0 atIndex:i];
            for (int c = c0; c < ncol && c < c0 + per; c++) {
                mparams_t P = { (float)eps, (float)lam, (float)kk, (float)aa, k0, k1, max_eh, col_off[c], col_off[c + 1] - col_off[c], 0 };
                if (P.count <= 0) continue;
                [enc setBytes:&P length:sizeof P atIndex:21];
                [enc dispatchThreads:MTLSizeMake(P.count, 1, 1) threadsPerThreadgroup:MTLSizeMake(mpipe.threadExecutionWidth, 1, 1)];
            }
            [enc endEncoding]; [cb commit]; [cb waitUntilCompleted];
            if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "Metal: command buffer failed: %s\n", [[cb.error localizedDescription] UTF8String]); exit(1); }
        }
    }
}
#endif

static void fb_alloc_cpu(void) {
    F.q = malloc(sizeof(float) * NE); F.V = malloc(sizeof(float) * NS); F.th = malloc(sizeof(float) * NS * 10);
    F.thsum = malloc(sizeof(float) * NT); F.area = malloc(sizeof(float) * NT); F.Vh = malloc(sizeof(float) * NT);
    F.dacc = calloc(NT, sizeof(float)); F.dvh = calloc(NT, sizeof(float));
    F.tmpV = malloc(sizeof(float) * nmax * 32); F.tmpth = malloc(sizeof(float) * nmax * 320); F.tmpA = malloc(sizeof(float) * nmax * max_eh);
    F.stat = calloc(nmax * 4, sizeof(int));
}
static void mirror(void) {       /* single-precision copies of the double-precision state */
    for (int e = 0; e < NE; e++) F.q[e] = (float)q[e];
    for (int s = 0; s < NS; s++) F.V[s] = (float)V[s];
    for (long i = 0; i < (long)NS * 10; i++) F.th[i] = (float)th[i];
    for (int t = 0; t < NT; t++) { F.thsum[t] = (float)thsum[t]; F.area[t] = (float)area[t]; F.Vh[t] = (float)Vh[t]; }
}
static void sweep_kd(unsigned int k0, unsigned int k1) {
    for (int c = 0; c < ncol; c++) { const int b0 = col_off[c], n = col_off[c + 1] - b0;
        for (int i = 0; i < n; i++) update_edge_kd(col_list[b0 + i], i, eps, lam, kk, aa, k0, k1, max_eh, q, V, th, thsum, area, Vh, dacc, dvh,
                                                   es_off, es_list, sedge, shinge, eh_off, eh_list, eh_cont, tedge, sV, sth, sA, sstat); }
}
static void sweep_kf(unsigned int k0, unsigned int k1) {
    for (int c = 0; c < ncol; c++) { const int b0 = col_off[c], n = col_off[c + 1] - b0;
        for (int i = 0; i < n; i++) update_edge_kf(col_list[b0 + i], i, (float)eps, (float)lam, (float)kk, (float)aa, k0, k1, max_eh, F.q, F.V, F.th, F.thsum, F.area, F.Vh, F.dacc, F.dvh,
                                                   es_off, es_list, sedge, shinge, eh_off, eh_list, eh_cont, tedge, F.tmpV, F.tmpth, F.tmpA, F.stat); }
}
static double take_acc(int *st) {
    long a = 0, t = 0;
    for (int i = 0; i < nmax; i++) { a += st[4*i]; t += st[4*i + 3]; st[4*i] = st[4*i + 1] = st[4*i + 2] = st[4*i + 3] = 0; }
    return t ? (double)a / (double)t : 0.0;
}
static int full(double *sv, double *sda, double *sa, double *md) { return par_on ? recompute_all_par(sv, sda, sa, md) : recompute_all(sv, sda, sa, md); }
static double action(void) { double sv, sda, sa, md; if (full(&sv, &sda, &sa, &md)) return NAN; return lam * sv - kk * sda + aa * full_sumR2; }

static void setup(int Lin, const char *cfg, double eps0) {
    L = Lin; init_tables();
    if (!build_lattice(0)) exit(1);
    eps = eps0;
    if (cfg && strcmp(cfg, "-")) {
        FILE *fl = fopen(cfg, "rb"); int Lf = 0; double ef;
        if (!fl || fread(&Lf, sizeof(int), 1, fl) != 1 || Lf != L || fread(&ef, sizeof(double), 1, fl) != 1 || fread(q, sizeof(double), NE, fl) != (size_t)NE) { fprintf(stderr, "bad config file\n"); exit(1); }
        fclose(fl); if (eps0 <= 0) eps = ef;
    }
    setenv("PARFINE", "1", 1);                           /* colors, not blocks: a color is one batch for the chip */
    par_on = 1; { const char *pr = getenv("PAR"); nthr = pr && atoi(pr) > 0 ? atoi(pr) : 1; }
    par_start();
    for (int c = 0; c < ncol; c++) if (col_off[c + 1] - col_off[c] > nmax) nmax = col_off[c + 1] - col_off[c];
    es_total = es_off[NE]; eh_total = eh_off[NE];
    sV = malloc(sizeof(double) * nmax * 32); sth = malloc(sizeof(double) * nmax * 320); sA = malloc(sizeof(double) * nmax * max_eh); sstat = calloc(nmax * 4, sizeof(int));
    printf("# L=%d NE=%d NS=%d NT=%d colors=%d (p=%d) largest color=%d max_eh=%d edge-simplex pairs=%ld edge-hinge pairs=%ld threads=%d\n",
           L, NE, NS, NT, ncol, colp, nmax, max_eh, es_total, eh_total, nthr);
    printf("# single-precision tables for the chip: %.2f GB\n", (4.0 * (NE + NS * 11.0 + NT * 5.0) + 4.0 * (2.0 * NE + es_total + 20.0 * NS + eh_total + 3.0 * NT + NE) + eh_total) / 1e9);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: check L k lambda a eps seed cfg|-   |   run kd|kf|gpu L k lambda a eps nsweeps seed out cfg|-\n"); return 1; }
    char dir[1024]; strncpy(dir, argv[0], sizeof dir - 1); dir[sizeof dir - 1] = 0; { char *sl = strrchr(dir, '/'); if (sl) *sl = 0; else strcpy(dir, "."); }
    if (!strcmp(argv[1], "check") && argc >= 9) {
        kk = atof(argv[3]); lam = atof(argv[4]); aa = atof(argv[5]);
        setup(atoi(argv[2]), argv[8], atof(argv[6])); fb_alloc_cpu();
        seed_rng((uint64_t)strtoull(argv[7], NULL, 10));
        double sv, sda, sa, md;
        if (full(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid start\n"); return 1; }
        printf("start: R=%.6f sumV=%.6f sumR2=%.6f eps=%.4f\n", 2 * sda / sv, sv, full_sumR2, eps);
        /* 0. a few sweeps to move away from the start, with the double kernel */
        for (int i = 0; i < 5; i++) { uint64_t key = next_u64(); sweep_kd((unsigned)key, (unsigned)(key >> 32)); }
        /* 1. the kernel keeps volumes, angles, angle sums, areas and hinge volumes in step with a full recomputation */
        double *cV = malloc(sizeof(double) * NS), *cth = malloc(sizeof(double) * NS * 10), *cts = malloc(sizeof(double) * NT), *car = malloc(sizeof(double) * NT), *cvh = malloc(sizeof(double) * NT);
        memcpy(cV, V, sizeof(double) * NS); memcpy(cth, th, sizeof(double) * NS * 10); memcpy(cts, thsum, sizeof(double) * NT); memcpy(car, area, sizeof(double) * NT); memcpy(cvh, Vh, sizeof(double) * NT);
        if (full(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid after sweeps\n"); return 1; }
        double m1 = 0;
        for (int s = 0; s < NS; s++) m1 = fmax(m1, fabs(cV[s] - V[s]));
        for (long i = 0; i < (long)NS * 10; i++) m1 = fmax(m1, fabs(cth[i] - th[i]));
        for (int t = 0; t < NT; t++) { m1 = fmax(m1, fabs(cts[t] - thsum[t])); m1 = fmax(m1, fabs(car[t] - area[t])); if (aa != 0.0) m1 = fmax(m1, fabs(cvh[t] - Vh[t])); }
        printf("CHECK 1  bookkeeping after 5 sweeps (double kernel) against full recomputation: largest difference %.3e  (%s)\n", m1, m1 < 1e-9 ? "PASS" : "FAIL");
        /* 2. the action change the kernel computes equals the brute-force change of the whole action */
        double I0 = action(), m2 = 0, m2f = 0; int n2 = 0; mirror();
        for (int trial = 0; trial < 400; trial++) {
            const int e = (int)(next_u64() % (uint64_t)NE); const unsigned k0 = (unsigned)next_u64(), k1 = (unsigned)next_u64();
            const double qold = q[e], qnew = qold + eps * (2.0 * u01_kd(k0, k1, (unsigned)e, 0) - 1.0);
            if (qnew <= 0) continue;
            q[e] = qnew; const double I1 = action(); q[e] = qold;
            if (I1 != I1) { action(); continue; }                     /* invalid proposal */
            action();                                                 /* arrays back to the old state */
            /* kernel in a private copy of the state: force no acceptance side effects by restoring afterwards */
            const double keep_q = q[e];
            update_edge_kd(e, 0, eps, lam, kk, aa, k0, k1, max_eh, q, V, th, thsum, area, Vh, dacc, dvh, es_off, es_list, sedge, shinge, eh_off, eh_list, eh_cont, tedge, sV, sth, sA, sstat);
            const double dId = sV[31]; q[e] = keep_q; action();
            update_edge_kf(e, 0, (float)eps, (float)lam, (float)kk, (float)aa, k0, k1, max_eh, F.q, F.V, F.th, F.thsum, F.area, F.Vh, F.dacc, F.dvh, es_off, es_list, sedge, shinge, eh_off, eh_list, eh_cont, tedge, F.tmpV, F.tmpth, F.tmpA, F.stat);
            const double dIf = F.tmpV[31]; mirror();
            m2 = fmax(m2, fabs(dId - (I1 - I0))); m2f = fmax(m2f, fabs(dIf - dId)); n2++;
        }
        take_acc(sstat); take_acc(F.stat);
        printf("CHECK 2  action change of %d proposals: kernel (double) against brute force, largest difference %.3e  (%s)\n", n2, m2, m2 < 1e-6 ? "PASS" : "FAIL");
        printf("CHECK 3  the same proposals in single precision against double: largest difference in the action change %.3e\n", m2f);
        /* 4. single against double precision over sweeps, same random numbers; single-precision copies refreshed every sweep */
        const int NSW = 20; double *q0 = malloc(sizeof(double) * NE), *snap = malloc(sizeof(double) * NE * NSW), Rd[64], accd[64]; uint64_t keys[64];
        memcpy(q0, q, sizeof(double) * NE);
        for (int i = 0; i < NSW; i++) { keys[i] = next_u64(); sweep_kd((unsigned)keys[i], (unsigned)(keys[i] >> 32)); accd[i] = take_acc(sstat);
            if (full(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid in double chain\n"); return 1; } Rd[i] = 2 * sda / sv; memcpy(snap + (long)i * NE, q, sizeof(double) * NE); }
        memcpy(q, q0, sizeof(double) * NE); full(&sv, &sda, &sa, &md);
        printf("CHECK 4  single against double precision, same start and random numbers:\n  sweep  R(double)   R(single)   acc(d)  acc(s)  edges that took a different decision (of %d)  largest difference among the rest\n", NE);
        for (int i = 0; i < NSW; i++) {
            mirror(); sweep_kf((unsigned)keys[i], (unsigned)(keys[i] >> 32)); const double accf = take_acc(F.stat);
            for (int e = 0; e < NE; e++) q[e] = (double)F.q[e];
            const int bad = full(&sv, &sda, &sa, &md);
            long ndiff = 0; double mrest = 0; const double *sn = snap + (long)i * NE;
            for (int e = 0; e < NE; e++) { const double d = fabs(q[e] - sn[e]); if (d > 1e-4 * eps) ndiff++; else if (d > mrest) mrest = d; }
            printf("  %5d  %+.6f  %+.6f  %.4f  %.4f  %8ld   %.2e%s\n", i + 1, Rd[i], 2 * sda / sv, accd[i], accf, ndiff, mrest, bad ? "  INVALID GEOMETRY" : "");
        }
        return 0;
    }
    if (!strcmp(argv[1], "run") && argc >= 12) {
        const char *kern = argv[2];
        kk = atof(argv[4]); lam = atof(argv[5]); aa = atof(argv[6]);
        setup(atoi(argv[3]), argv[11], atof(argv[7]));
        const int nsw = atoi(argv[8]); seed_rng((uint64_t)strtoull(argv[9], NULL, 10));
        const int mode = !strcmp(kern, "kd") ? 0 : !strcmp(kern, "kf") ? 1 : 2;
        if (mode == 1) fb_alloc_cpu();
        if (mode == 2) {
#ifdef USE_METAL
            metal_setup(dir);
#else
            fprintf(stderr, "this build has no Metal\n"); return 1;
#endif
        }
        double sv, sda, sa, md;
        if (full(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid start\n"); return 1; }
        FILE *fo = fopen(argv[10], "w");
        fprintf(fo, "# kernel=%s L=%d k=%.6f lambda=%.6f a=%.6f eps=%.5f seed=%s colors=%d largest color=%d\n# cols: sweep sumV sumDA R meanq acc sumR2 sweep_sec measure_sec\n", kern, L, kk, lam, aa, eps, argv[9], ncol, nmax);
        double tsw = 0, tme = 0;
        for (int i = 0; i < nsw; i++) {
            const uint64_t key = next_u64(); double acc = 0; const double w0 = wall_sec();
            if (mode == 0) { sweep_kd((unsigned)key, (unsigned)(key >> 32)); acc = take_acc(sstat); }
            else {
                mirror(); const double w0b = wall_sec(); (void)w0b;
                if (mode == 1) sweep_kf((unsigned)key, (unsigned)(key >> 32));
#ifdef USE_METAL
                else sweep_gpu((unsigned)key, (unsigned)(key >> 32));
#endif
                acc = take_acc(F.stat);
                for (int e = 0; e < NE; e++) q[e] = (double)F.q[e];
            }
            const double w1 = wall_sec();
            if (full(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid geometry after sweep %d\n", i); return 2; }
            double sq = 0; for (int e = 0; e < NE; e++) sq += q[e];
            const double w2 = wall_sec(); tsw += w1 - w0; tme += w2 - w1;
            fprintf(fo, "%d %.10g %.10g %.8f %.8f %.5f %.10g %.4f %.4f\n", i, sv, sda, 2 * sda / sv, sq / NE, acc, full_sumR2, w1 - w0, w2 - w1);
            if ((i + 1) % 50 == 0) fflush(fo);
        }
        fclose(fo);
        printf("kernel %s L=%d: %d sweeps, %.4f s per sweep (update) + %.4f s (recompute in double); last R=%.6f\n", kern, L, nsw, tsw / nsw, tme / nsw, 2 * sda / sv);
        if (getenv("SAVE")) { FILE *fs = fopen(getenv("SAVE"), "wb"); fwrite(&L, sizeof(int), 1, fs); fwrite(&eps, sizeof(double), 1, fs); fwrite(q, sizeof(double), NE, fs); fclose(fs); }
        par_stop();
        return 0;
    }
    fprintf(stderr, "bad arguments\n"); return 1;
}
