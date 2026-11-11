/*
 * regge.c -- 4D Euclidean quantum Regge calculus Monte Carlo on a periodic
 * Freudenthal (Kuhn) triangulated hypercubic lattice Z_L^4.
 *
 * Variables: squared edge lengths q_e > 0, uniform measure prod dq_e restricted
 * to configurations in which every 4-simplex has a positive-definite Gram matrix.
 * Action (Hamber convention):  I = lambda * sum_s V_s - k * sum_h delta_h A_h + a * sum_h delta_h^2 A_h^2 / V_h
 * with V_h = sum_{s containing h} V_s / 10. The coupling a is read from the environment variable R2A
 * (unset or 0: the a = 0 engine, identical moves, random numbers and output).
 *
 * Modes:
 *   ./regge selftest L
 *   ./regge run L k lambda ntherm nmeas seed outfile [eps0] [loadcfg|-] [savecfg|-] [nscale] [nphi] [freeze]
 *   ./regge phitest L nphi nsweeps seed [cfg|-]   (scalar equipartition check on a fixed geometry)
 *   nphi = number of free massless scalars on vertices (heat bath, coupled to geometry)
 *   freeze = 0 none, 1 edges along e_0 alone fixed, 2 all edges whose step set contains axis 0 fixed
 *   (nscale = global dilation moves per sweep, default 0 = pure single-edge Metropolis)
 *
 * Build: gcc -O2 -o regge_r2 regge_r2.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>
#include <stdatomic.h>
#include <sched.h>

#define TWOPI 6.283185307179586476925286766559
#define MAXES 32          /* max simplices per edge (24 for this lattice) */

static int L, N0, NE, NT, NS;
static int PIDX[5][5];          /* local pair (i,j) -> 0..9 */
static int PA[10], PB[10];      /* pair endpoints */
static int TRI_OF_PAIR[10][3];  /* dihedral pair (p,r) -> complementary triple */
static int tricode[16][16];     /* (S1,S2) disjoint nonempty masks -> 0..49 */

static double *q;               /* [NE] squared edge lengths */
static int *svert;              /* [NS*5] chain-ordered vertex ids */
static int *sedge;              /* [NS*10] edge id of local pair k */
static int *shinge;             /* [NS*10] triangle id at which dihedral angle k sits */
static int *tedge;              /* [NT*3] edge ids of triangle */
static double *V, *th, *thsum, *area, *dacc;
static int *es_off, *es_list;   /* edge -> simplices (CSR) */
static int *eh_off, *eh_list;   /* edge -> affected hinges (CSR) */
static unsigned char *eh_cont;  /* does that hinge contain the edge? */
static int max_eh = 0;
static double *tmpA;
static double aa = 0.0;                 /* curvature-squared coupling a (env R2A) */
static double *Vh = NULL, *dvh = NULL;  /* [NT] hinge volumes V_h and pending changes */
static double run_sumR2 = 0.0, full_sumR2 = 0.0;   /* sum_h delta^2 A^2 / V_h */

/* ---------------- RNG: xoshiro256** ---------------- */
static uint64_t rs[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next_u64(void) {
    const uint64_t result = rotl(rs[1] * 5, 7) * 9;
    const uint64_t t = rs[1] << 17;
    rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3];
    rs[2] ^= t; rs[3] = rotl(rs[3], 45);
    return result;
}
static inline double urand(void) { return (double)(next_u64() >> 11) * 0x1.0p-53; }
static void seed_rng(uint64_t seed) {
    uint64_t z = seed;
    for (int i = 0; i < 4; i++) {
        z += 0x9e3779b97f4a7c15ULL;
        uint64_t x = z;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        rs[i] = x ^ (x >> 31);
    }
}

/* ---------------- geometry ---------------- */
static const double CHOL_TOL = 1e-12;   /* pivot must exceed CHOL_TOL * G_ii */

/* qq[10]: squared lengths of local pairs. Returns 0 if Gram matrix not positive definite. */
static int simplex_geom(const double *qq, double *vol, double *ang) {
    double G[4][4];
    for (int a = 1; a <= 4; a++) G[a-1][a-1] = qq[PIDX[0][a]];
    for (int a = 1; a <= 4; a++)
        for (int b = a + 1; b <= 4; b++) {
            double g = 0.5 * (qq[PIDX[0][a]] + qq[PIDX[0][b]] - qq[PIDX[a][b]]);
            G[a-1][b-1] = G[b-1][a-1] = g;
        }
    double Lm[4][4], prod = 1.0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
            double s = G[i][j];
            for (int k = 0; k < j; k++) s -= Lm[i][k] * Lm[j][k];
            if (i == j) {
                if (!(s > CHOL_TOL * G[i][i])) return 0;
                Lm[i][i] = sqrt(s);
                prod *= Lm[i][i];
            } else Lm[i][j] = s / Lm[j][j];
        }
    *vol = prod / 24.0;
    double Li[4][4];
    for (int i = 0; i < 4; i++) {
        Li[i][i] = 1.0 / Lm[i][i];
        for (int j = 0; j < i; j++) {
            double s = 0.0;
            for (int k = j; k < i; k++) s += Lm[i][k] * Li[k][j];
            Li[i][j] = -s * Li[i][i];
        }
        for (int j = i + 1; j < 4; j++) Li[i][j] = 0.0;
    }
    double Gi[4][4];
    for (int a = 0; a < 4; a++)
        for (int b = a; b < 4; b++) {
            double s = 0.0;
            for (int k = b; k < 4; k++) s += Li[k][a] * Li[k][b];
            Gi[a][b] = Gi[b][a] = s;
        }
    /* N_p.N_r for outward face normals (covectors) */
    double M[5][5], tot = 0.0;
    for (int a = 0; a < 4; a++) {
        double cs = 0.0;
        for (int c = 0; c < 4; c++) { M[a+1][c+1] = Gi[a][c]; cs += Gi[a][c]; }
        M[0][a+1] = M[a+1][0] = -cs;
        tot += cs;
    }
    M[0][0] = tot;
    double inv[5];
    for (int p = 0; p < 5; p++) inv[p] = 1.0 / sqrt(M[p][p]);
    for (int k = 0; k < 10; k++) {
        int p = PA[k], r = PB[k];
        double c = -M[p][r] * inv[p] * inv[r];
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;
        ang[k] = acos(c);
    }
    return 1;
}

static inline double tri_area(double q1, double q2, double q3) {
    double u = q1 + q2 - q3;
    double t = 4.0 * q1 * q2 - u * u;
    return t > 0.0 ? 0.25 * sqrt(t) : 0.0;
}

/* ---------------- lattice construction ---------------- */
static int vid(const int c[4]) {
    int x[4];
    for (int i = 0; i < 4; i++) x[i] = ((c[i] % L) + L) % L;
    return x[0] + L * (x[1] + L * (x[2] + L * x[3]));
}
static void vcoord(int v, int c[4]) { for (int i = 0; i < 4; i++) { c[i] = v % L; v /= L; } }

static int perms[24][4];
static void make_perms(void) {
    int n = 0;
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) for (int c = 0; c < 4; c++) for (int d = 0; d < 4; d++)
        if (a != b && a != c && a != d && b != c && b != d && c != d) {
            perms[n][0] = a; perms[n][1] = b; perms[n][2] = c; perms[n][3] = d; n++;
        }
}

static void init_tables(void) {
    int k = 0;
    for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) {
        PIDX[i][j] = PIDX[j][i] = k; PA[k] = i; PB[k] = j; k++;
    }
    for (int k2 = 0; k2 < 10; k2++) {
        int n = 0;
        for (int v = 0; v < 5; v++) if (v != PA[k2] && v != PB[k2]) TRI_OF_PAIR[k2][n++] = v;
    }
    int c = 0;
    for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) tricode[a][b] = -1;
    for (int a = 1; a < 16; a++) for (int b = 1; b < 16; b++) if ((a & b) == 0) tricode[a][b] = c++;
    if (c != 50) { fprintf(stderr, "tricode count %d != 50\n", c); exit(1); }
    make_perms();
}

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y);
}
static long count_unique(uint64_t *a, long n) {
    qsort(a, n, sizeof(uint64_t), cmp_u64);
    long u = 0;
    for (long i = 0; i < n; i++) if (i == 0 || a[i] != a[i-1]) u++;
    return u;
}
static void sort_small(int *a, int n) {
    for (int i = 1; i < n; i++) { int t = a[i], j = i - 1; while (j >= 0 && a[j] > t) { a[j+1] = a[j]; j--; } a[j+1] = t; }
}

/* Build everything. If check != 0, verify counts by brute-force unique vertex sets. */
static int build_lattice(int check) {
    N0 = L * L * L * L; NE = 15 * N0; NT = 50 * N0; NS = 24 * N0;
    q = malloc(sizeof(double) * NE);
    svert = malloc(sizeof(int) * NS * 5);
    sedge = malloc(sizeof(int) * NS * 10);
    shinge = malloc(sizeof(int) * NS * 10);
    tedge = malloc(sizeof(int) * NT * 3);
    for (int i = 0; i < NT * 3; i++) tedge[i] = -1;
    int *edge_key_u = malloc(sizeof(int) * NE * 2);
    for (int i = 0; i < NE * 2; i++) edge_key_u[i] = -1;
    int ok = 1;
    for (int v = 0; v < N0; v++) {
        int c0[4]; vcoord(v, c0);
        for (int p = 0; p < 24; p++) {
            int s = v * 24 + p;
            int cv[5][4], m[5];
            m[0] = 0;
            for (int i = 0; i < 4; i++) cv[0][i] = c0[i];
            for (int a = 1; a <= 4; a++) {
                for (int i = 0; i < 4; i++) cv[a][i] = cv[a-1][i];
                cv[a][perms[p][a-1]] += 1;
                m[a] = m[a-1] | (1 << perms[p][a-1]);
            }
            int gv[5];
            for (int a = 0; a < 5; a++) { gv[a] = vid(cv[a]); svert[s*5 + a] = gv[a]; }
            for (int k = 0; k < 10; k++) {
                int a = PA[k], b = PB[k];
                int mask = m[b] ^ m[a];
                int e = 15 * gv[a] + mask - 1;
                sedge[s*10 + k] = e;
                q[e] = (double)__builtin_popcount(mask);
                int u0 = gv[a] < gv[b] ? gv[a] : gv[b], u1 = gv[a] < gv[b] ? gv[b] : gv[a];
                if (edge_key_u[2*e] < 0) { edge_key_u[2*e] = u0; edge_key_u[2*e+1] = u1; }
                else if (edge_key_u[2*e] != u0 || edge_key_u[2*e+1] != u1) ok = 0;
            }
            for (int k = 0; k < 10; k++) {
                int a = TRI_OF_PAIR[k][0], b = TRI_OF_PAIR[k][1], cc = TRI_OF_PAIR[k][2];
                int S1 = m[b] ^ m[a], S2 = m[cc] ^ m[b];
                int t = 50 * gv[a] + tricode[S1][S2];
                shinge[s*10 + k] = t;
                int te[3] = { sedge[s*10 + PIDX[a][b]], sedge[s*10 + PIDX[a][cc]], sedge[s*10 + PIDX[b][cc]] };
                if (tedge[3*t] < 0) { tedge[3*t] = te[0]; tedge[3*t+1] = te[1]; tedge[3*t+2] = te[2]; }
                else if (tedge[3*t] != te[0] || tedge[3*t+1] != te[1] || tedge[3*t+2] != te[2]) ok = 0;
            }
        }
    }
    for (int e = 0; e < NE; e++) if (edge_key_u[2*e] < 0) ok = 0;
    for (int t = 0; t < NT; t++) if (tedge[3*t] < 0) ok = 0;
    free(edge_key_u);
    if (!ok) { fprintf(stderr, "canonical id construction inconsistent\n"); return 0; }

    if (check) {
        /* brute force: distinct vertex sets of edges, triangles, tetrahedra, 4-simplices */
        uint64_t *ke = malloc(sizeof(uint64_t) * NS * 10);
        uint64_t *kt = malloc(sizeof(uint64_t) * NS * 10);
        uint64_t *k3 = malloc(sizeof(uint64_t) * NS * 5);
        uint64_t *k4 = malloc(sizeof(uint64_t) * NS);
        for (int s = 0; s < NS; s++) {
            int *gv = &svert[s*5];
            for (int k = 0; k < 10; k++) {
                int w[2] = { gv[PA[k]], gv[PB[k]] }; sort_small(w, 2);
                ke[s*10 + k] = ((uint64_t)w[0] << 16) | (uint64_t)w[1];
                int t3[3] = { gv[TRI_OF_PAIR[k][0]], gv[TRI_OF_PAIR[k][1]], gv[TRI_OF_PAIR[k][2]] }; sort_small(t3, 3);
                kt[s*10 + k] = ((uint64_t)t3[0] << 32) | ((uint64_t)t3[1] << 16) | (uint64_t)t3[2];
            }
            for (int o = 0; o < 5; o++) {
                int w[4], n = 0;
                for (int a = 0; a < 5; a++) if (a != o) w[n++] = gv[a];
                sort_small(w, 4);
                k3[s*5 + o] = ((uint64_t)w[0] << 48) | ((uint64_t)w[1] << 32) | ((uint64_t)w[2] << 16) | (uint64_t)w[3];
            }
            int w5[5]; for (int a = 0; a < 5; a++) w5[a] = gv[a];
            sort_small(w5, 5);
            for (int a = 0; a < 4; a++) if (w5[a] == w5[a+1]) ok = 0;   /* degenerate simplex */
            /* 5 ids < 2^12 each for L <= 8 */
            k4[s] = ((uint64_t)w5[0] << 48) | ((uint64_t)w5[1] << 36) | ((uint64_t)w5[2] << 24) | ((uint64_t)w5[3] << 12) | (uint64_t)w5[4];
        }
        long ue = count_unique(ke, (long)NS * 10), ut = count_unique(kt, (long)NS * 10);
        long u3 = count_unique(k3, (long)NS * 5), u4 = count_unique(k4, NS);
        printf("COUNTS L=%d N0=%d  edges=%ld (%.4f/vertex)  triangles=%ld (%.4f/vertex)  tetrahedra=%ld (%.4f/vertex)  4-simplices=%ld (%.4f/vertex)\n",
               L, N0, ue, (double)ue / N0, ut, (double)ut / N0, u3, (double)u3 / N0, u4, (double)u4 / N0);
        printf("EULER chi = %ld\n", (long)N0 - ue + ut - u3 + u4);
        int good = (ue == 15L * N0 && ut == 50L * N0 && u3 == 60L * N0 && u4 == 24L * N0 && ok);
        printf("COUNTS_OK %d (canonical ids: %d edges, %d triangles)\n", good, NE, NT);
        free(ke); free(kt); free(k3); free(k4);
        if (!good) return 0;
    }

    /* edge -> simplices */
    es_off = calloc(NE + 1, sizeof(int));
    for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) es_off[sedge[s*10 + k] + 1]++;
    for (int e = 0; e < NE; e++) es_off[e+1] += es_off[e];
    es_list = malloc(sizeof(int) * es_off[NE]);
    int *fill = calloc(NE, sizeof(int));
    for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) { int e = sedge[s*10 + k]; es_list[es_off[e] + fill[e]++] = s; }
    free(fill);
    int maxes = 0;
    for (int e = 0; e < NE; e++) if (es_off[e+1] - es_off[e] > maxes) maxes = es_off[e+1] - es_off[e];
    if (maxes > MAXES) { fprintf(stderr, "too many simplices per edge\n"); exit(1); }

    /* edge -> affected hinges (union of triangles of its simplices) */
    int *mark = malloc(sizeof(int) * NT);
    for (int t = 0; t < NT; t++) mark[t] = -1;
    eh_off = calloc(NE + 1, sizeof(int));
    for (int e = 0; e < NE; e++) {
        int n = 0;
        for (int i = es_off[e]; i < es_off[e+1]; i++) { int s = es_list[i];
            for (int k = 0; k < 10; k++) { int t = shinge[s*10 + k]; if (mark[t] != e) { mark[t] = e; n++; } } }
        eh_off[e+1] = eh_off[e] + n;
        if (n > max_eh) max_eh = n;
    }
    eh_list = malloc(sizeof(int) * eh_off[NE]);
    eh_cont = malloc(eh_off[NE]);
    for (int t = 0; t < NT; t++) mark[t] = -1;
    for (int e = 0; e < NE; e++) {
        int n = eh_off[e];
        for (int i = es_off[e]; i < es_off[e+1]; i++) { int s = es_list[i];
            for (int k = 0; k < 10; k++) { int t = shinge[s*10 + k];
                if (mark[t] != e) { mark[t] = e; eh_list[n] = t;
                    eh_cont[n] = (tedge[3*t] == e || tedge[3*t+1] == e || tedge[3*t+2] == e); n++; } } }
    }
    free(mark);
    /* sanity: every triangle containing e must be in e's list */
    {
        long ncont = 0; for (long i = 0; i < eh_off[NE]; i++) ncont += eh_cont[i];
        if (ncont != 3L * NT) { fprintf(stderr, "hinge containment count %ld != %ld\n", ncont, 3L * NT); exit(1); }
    }
    V = malloc(sizeof(double) * NS); th = malloc(sizeof(double) * NS * 10);
    thsum = malloc(sizeof(double) * NT); area = malloc(sizeof(double) * NT);
    dacc = calloc(NT, sizeof(double)); tmpA = malloc(sizeof(double) * max_eh);
    Vh = calloc(NT, sizeof(double)); dvh = calloc(NT, sizeof(double));
    if (check) {
        int hist[33] = {0};
        for (int e = 0; e < NE; e++) hist[es_off[e+1] - es_off[e]]++;
        printf("SIMPLICES_PER_EDGE:");
        for (int i = 0; i <= 32; i++) if (hist[i]) printf(" %d:%d", i, hist[i]);
        printf("   max affected hinges per edge = %d\n", max_eh);
        int th_hist[64] = {0};
        int *cnt = calloc(NT, sizeof(int));
        for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) cnt[shinge[s*10 + k]]++;
        for (int t = 0; t < NT; t++) th_hist[cnt[t] < 63 ? cnt[t] : 63]++;
        printf("SIMPLICES_PER_TRIANGLE:");
        for (int i = 0; i < 64; i++) if (th_hist[i]) printf(" %d:%d", i, th_hist[i]);
        printf("\n");
        free(cnt);
    }
    return 1;
}

/* full recomputation of V, th, thsum, area; returns number of invalid simplices */
static int recompute_all(double *sumV, double *sumDA, double *sumA, double *maxabsdelta) {
    double qq[10];
    int bad = 0;
    for (int s = 0; s < NS; s++) {
        for (int k = 0; k < 10; k++) qq[k] = q[sedge[s*10 + k]];
        if (!simplex_geom(qq, &V[s], &th[s*10])) bad++;
    }
    memset(thsum, 0, sizeof(double) * NT);
    for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) thsum[shinge[s*10 + k]] += th[s*10 + k];
    double sv = 0, sda = 0, sa = 0, md = 0;
    for (int s = 0; s < NS; s++) sv += V[s];
    for (int t = 0; t < NT; t++) {
        area[t] = tri_area(q[tedge[3*t]], q[tedge[3*t+1]], q[tedge[3*t+2]]);
        double d = TWOPI - thsum[t];
        sda += d * area[t]; sa += area[t];
        if (fabs(d) > md) md = fabs(d);
    }
    if (aa != 0.0) {
        memset(Vh, 0, sizeof(double) * NT);
        for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) Vh[shinge[s*10 + k]] += 0.1 * V[s];
        double r2 = 0.0;
        for (int t = 0; t < NT; t++) { double d = TWOPI - thsum[t]; if (Vh[t] > 0.0) r2 += d * d * area[t] * area[t] / Vh[t]; }
        full_sumR2 = r2;
    }
    if (sumV) *sumV = sv;
    if (sumDA) *sumDA = sda;
    if (sumA) *sumA = sa;
    if (maxabsdelta) *maxabsdelta = md;
    return bad;
}

/* ---------------- Metropolis ---------------- */
static double eps = 0.5, lam = 1.0, kk = 0.0;
static long n_acc = 0, n_inval = 0, n_neg = 0, n_try = 0;
static double run_sumV = 0, run_sumDA = 0;   /* incrementally tracked totals */
static double tmpV[MAXES], tmpth[MAXES * 10];

/* ---------------- free massless scalar fields on vertices ----------------
 * S_phi = (1/2) sum_edges W_ij sum_f (phi_f,i - phi_f,j)^2,  W_ij = V_ij / q_ij,
 * V_ij = sum_{s containing ij} V_s / 10.  Exact heat bath per vertex; zero mode removed by
 * subtracting each field's mean after every heat-bath sweep (the action depends only on differences). */
static int nphi = 0;
static int seq_on = 0;
static int lines_on = 0;   /* SA74: LINES=1 appends Wilson-line data: lengths of the closed lines along axis 3 and their correlations */
static int slices_on = 0;  /* SLICES=1 appends per-slice curvature sums */   /* SA74 night: SEQ=1 sweeps the edges in index order (vertex by vertex) from a random offset */
static double *phi = NULL;      /* [nphi*N0] */
static double *Phi2 = NULL;     /* [NE] sum_f (phi_i - phi_j)^2 */
static double *Vedge = NULL;    /* [NE] */
static int *eu = NULL, *ev = NULL;
static int *vadj_off = NULL, *vadj_e = NULL, *vadj_n = NULL;
static double run_Sphi = 0.0;

static inline double gauss(void) {
    double u1 = 1.0 - urand(), u2 = urand();
    return sqrt(-2.0 * log(u1)) * cos(TWOPI * u2);
}
static void build_scalar_tables(void) {
    eu = malloc(sizeof(int) * NE); ev = malloc(sizeof(int) * NE);
    for (int e = 0; e < NE; e++) {
        int base = e / 15, mask = e % 15 + 1, c[4];
        vcoord(base, c);
        for (int i = 0; i < 4; i++) if (mask & (1 << i)) c[i] += 1;
        eu[e] = base; ev[e] = vid(c);
    }
    vadj_off = calloc(N0 + 1, sizeof(int));
    for (int e = 0; e < NE; e++) { vadj_off[eu[e] + 1]++; vadj_off[ev[e] + 1]++; }
    for (int i = 0; i < N0; i++) vadj_off[i+1] += vadj_off[i];
    vadj_e = malloc(sizeof(int) * vadj_off[N0]); vadj_n = malloc(sizeof(int) * vadj_off[N0]);
    int *fill = calloc(N0, sizeof(int));
    for (int e = 0; e < NE; e++) {
        int u = eu[e], v = ev[e];
        vadj_e[vadj_off[u] + fill[u]] = e; vadj_n[vadj_off[u] + fill[u]++] = v;
        vadj_e[vadj_off[v] + fill[v]] = e; vadj_n[vadj_off[v] + fill[v]++] = u;
    }
    free(fill);
    Phi2 = calloc(NE, sizeof(double)); Vedge = calloc(NE, sizeof(double));
    phi = calloc((size_t)(nphi > 0 ? nphi : 1) * N0, sizeof(double));
}
static void compute_Vedge(void) {
    memset(Vedge, 0, sizeof(double) * NE);
    for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) Vedge[sedge[s*10 + k]] += 0.1 * V[s];
}
static void compute_Phi2(void) {
    for (int e = 0; e < NE; e++) {
        double s = 0.0;
        for (int f = 0; f < nphi; f++) { double d = phi[f*N0 + eu[e]] - phi[f*N0 + ev[e]]; s += d * d; }
        Phi2[e] = s;
    }
}
static double compute_Sphi(void) {     /* needs Vedge and Phi2 current */
    double S = 0.0;
    for (int e = 0; e < NE; e++) S += Vedge[e] / q[e] * Phi2[e];
    return 0.5 * S;
}
static void heatbath_sweep(void) {
    compute_Vedge();
    double b[16];
    for (int i = 0; i < N0; i++) {
        double a = 0.0;
        for (int f = 0; f < nphi; f++) b[f] = 0.0;
        for (int j = vadj_off[i]; j < vadj_off[i+1]; j++) {
            const int e = vadj_e[j], n = vadj_n[j];
            const double w = Vedge[e] / q[e];
            a += w;
            for (int f = 0; f < nphi; f++) b[f] += w * phi[f*N0 + n];
        }
        const double sa = 1.0 / sqrt(a);
        for (int f = 0; f < nphi; f++) phi[f*N0 + i] = b[f] / a + gauss() * sa;
    }
    for (int f = 0; f < nphi; f++) {           /* remove zero mode */
        double m = 0.0;
        for (int i = 0; i < N0; i++) m += phi[f*N0 + i];
        m /= N0;
        for (int i = 0; i < N0; i++) phi[f*N0 + i] -= m;
    }
    compute_Phi2();
    run_Sphi = compute_Sphi();
}

static void update_edge(int e) {
    n_try++;
    double qold = q[e];
    double qnew = qold + eps * (2.0 * urand() - 1.0);
    if (qnew <= 0.0) { n_neg++; return; }
    q[e] = qnew;
    const int o0 = es_off[e], o1 = es_off[e+1];
    double dV = 0.0, qq[10];
    for (int i = o0; i < o1; i++) {
        const int s = es_list[i];
        const int *se = &sedge[s*10];
        for (int k = 0; k < 10; k++) qq[k] = q[se[k]];
        if (!simplex_geom(qq, &tmpV[i - o0], &tmpth[(i - o0) * 10])) { q[e] = qold; n_inval++; return; }
        dV += tmpV[i - o0] - V[s];
    }
    for (int i = o0; i < o1; i++) {
        const int s = es_list[i];
        for (int k = 0; k < 10; k++) dacc[shinge[s*10 + k]] += tmpth[(i - o0) * 10 + k] - th[s*10 + k];
    }
    if (aa != 0.0)
        for (int i = o0; i < o1; i++) {
            const int s = es_list[i];
            const double dv10 = 0.1 * (tmpV[i - o0] - V[s]);
            for (int k = 0; k < 10; k++) dvh[shinge[s*10 + k]] += dv10;
        }
    const int h0 = eh_off[e], h1 = eh_off[e+1];
    double dDA = 0.0, dR2 = 0.0;
    for (int j = h0; j < h1; j++) {
        const int h = eh_list[j];
        const double dold = TWOPI - thsum[h];
        const double dnew = dold - dacc[h];
        const double Aold = area[h];
        const double Anew = eh_cont[j] ? tri_area(q[tedge[3*h]], q[tedge[3*h+1]], q[tedge[3*h+2]]) : Aold;
        tmpA[j - h0] = Anew;
        dDA += dnew * Anew - dold * Aold;
        if (aa != 0.0) {
            const double vo = Vh[h], vn = vo + dvh[h];
            dR2 += dnew * dnew * Anew * Anew / vn - dold * dold * Aold * Aold / vo;
        }
    }
    double dSphi = 0.0;
    if (nphi > 0) {      /* change of the scalar action: W_ij = V_ij/q_ij changes on every edge of the touched simplices */
        double Ve_old = 0.0, dVe = 0.0;
        for (int i = o0; i < o1; i++) {
            const int s = es_list[i];
            const double dV10 = 0.1 * (tmpV[i - o0] - V[s]);
            Ve_old += 0.1 * V[s]; dVe += dV10;
            const int *se = &sedge[s*10];
            double acc = 0.0;
            for (int k = 0; k < 10; k++) { const int f = se[k]; if (f != e) acc += Phi2[f] / q[f]; }
            dSphi += dV10 * acc;
        }
        dSphi += ((Ve_old + dVe) / qnew - Ve_old / qold) * Phi2[e];
        dSphi *= 0.5;
    }
    double dI = lam * dV - kk * dDA + dSphi;
    if (aa != 0.0) dI += aa * dR2;
    const int accept = (dI <= 0.0) || (urand() < exp(-dI));
    if (accept) {
        for (int i = o0; i < o1; i++) {
            const int s = es_list[i];
            V[s] = tmpV[i - o0];
            memcpy(&th[s*10], &tmpth[(i - o0) * 10], 10 * sizeof(double));
        }
        for (int j = h0; j < h1; j++) { const int h = eh_list[j]; thsum[h] += dacc[h]; area[h] = tmpA[j - h0]; }
        run_sumV += dV; run_sumDA += dDA; run_Sphi += dSphi;
        if (aa != 0.0) { for (int j = h0; j < h1; j++) Vh[eh_list[j]] += dvh[eh_list[j]]; run_sumR2 += dR2; }
        n_acc++;
    } else q[e] = qold;
    for (int j = h0; j < h1; j++) dacc[eh_list[j]] = 0.0;
    if (aa != 0.0) for (int j = h0; j < h1; j++) dvh[eh_list[j]] = 0.0;
}


/* ---------------- PAR=n: colored parallel sweep (SA74, 7 November 2026) ----------------
 * Changing edge e alters only the hinges in eh_list[e] (the triangles of the simplices that contain e) and reads only data of those
 * simplices and hinges, so two edges whose hinge lists are disjoint can be updated at the same time, in either order. The colors are
 * (edge type tau = e mod 15, base-vertex coordinates mod p) with p the smallest divisor of L that is at least 3; setup_colors checks for
 * the lattice in use that the edges of each color have pairwise disjoint hinge lists. A sweep visits the 15 p^4 colors in turn and the
 * threads split each color's edges. Random numbers come from a counter hash of (sweep key, edge), so a run does not depend on the number
 * of threads; the full recomputation after each sweep sums in fixed chunks for the same reason. */
static int par_on = 0, nthr = 1, ncol = 0, colp = 0, blk_b = 0, nblk = 0;
static int *col_off = NULL, *col_list = NULL, *blk_off = NULL, *blk_list = NULL, *cblk_off = NULL, *cblk_list = NULL;
typedef struct { double tmpV[MAXES], tmpth[MAXES * 10]; double *tmpA; long n_acc, n_inval, n_neg, n_try; double dV, dDA, dSphi, dR2; char pad[64]; } tstate_t;
static tstate_t *TS = NULL;
static uint64_t sweep_key = 0;
static inline uint64_t mix64(uint64_t z) { z ^= z >> 30; z *= 0xbf58476d1ce4e5b9ULL; z ^= z >> 27; z *= 0x94d049bb133111ebULL; return z ^ (z >> 31); }
static inline double ctr_u01(uint64_t key, int e, int j) {
    return (double)(mix64(key ^ mix64(((uint64_t)(unsigned)e << 2) + (uint64_t)j + 0x9e3779b97f4a7c15ULL)) >> 11) * 0x1.0p-53;
}
static void update_edge_par(int e, tstate_t *ts) {
    ts->n_try++;
    double qold = q[e];
    double qnew = qold + eps * (2.0 * ctr_u01(sweep_key, e, 0) - 1.0);
    if (qnew <= 0.0) { ts->n_neg++; return; }
    q[e] = qnew;
    const int o0 = es_off[e], o1 = es_off[e+1];
    double dV = 0.0, qq[10];
    for (int i = o0; i < o1; i++) {
        const int s = es_list[i];
        const int *se = &sedge[s*10];
        for (int k = 0; k < 10; k++) qq[k] = q[se[k]];
        if (!simplex_geom(qq, &ts->tmpV[i - o0], &ts->tmpth[(i - o0) * 10])) { q[e] = qold; ts->n_inval++; return; }
        dV += ts->tmpV[i - o0] - V[s];
    }
    for (int i = o0; i < o1; i++) {
        const int s = es_list[i];
        for (int k = 0; k < 10; k++) dacc[shinge[s*10 + k]] += ts->tmpth[(i - o0) * 10 + k] - th[s*10 + k];
    }
    if (aa != 0.0)
        for (int i = o0; i < o1; i++) {
            const int s = es_list[i];
            const double dv10 = 0.1 * (ts->tmpV[i - o0] - V[s]);
            for (int k = 0; k < 10; k++) dvh[shinge[s*10 + k]] += dv10;
        }
    const int h0 = eh_off[e], h1 = eh_off[e+1];
    double dDA = 0.0, dR2 = 0.0;
    for (int j = h0; j < h1; j++) {
        const int h = eh_list[j];
        const double dold = TWOPI - thsum[h];
        const double dnew = dold - dacc[h];
        const double Aold = area[h];
        const double Anew = eh_cont[j] ? tri_area(q[tedge[3*h]], q[tedge[3*h+1]], q[tedge[3*h+2]]) : Aold;
        ts->tmpA[j - h0] = Anew;
        dDA += dnew * Anew - dold * Aold;
        if (aa != 0.0) {
            const double vo = Vh[h], vn = vo + dvh[h];
            dR2 += dnew * dnew * Anew * Anew / vn - dold * dold * Aold * Aold / vo;
        }
    }
    double dSphi = 0.0;
    if (nphi > 0) {      /* change of the scalar action: W_ij = V_ij/q_ij changes on every edge of the touched simplices */
        double Ve_old = 0.0, dVe = 0.0;
        for (int i = o0; i < o1; i++) {
            const int s = es_list[i];
            const double dV10 = 0.1 * (ts->tmpV[i - o0] - V[s]);
            Ve_old += 0.1 * V[s]; dVe += dV10;
            const int *se = &sedge[s*10];
            double acc = 0.0;
            for (int k = 0; k < 10; k++) { const int f = se[k]; if (f != e) acc += Phi2[f] / q[f]; }
            dSphi += dV10 * acc;
        }
        dSphi += ((Ve_old + dVe) / qnew - Ve_old / qold) * Phi2[e];
        dSphi *= 0.5;
    }
    double dI = lam * dV - kk * dDA + dSphi;
    if (aa != 0.0) dI += aa * dR2;
    const int accept = (dI <= 0.0) || (ctr_u01(sweep_key, e, 1) < exp(-dI));
    if (accept) {
        for (int i = o0; i < o1; i++) {
            const int s = es_list[i];
            V[s] = ts->tmpV[i - o0];
            memcpy(&th[s*10], &ts->tmpth[(i - o0) * 10], 10 * sizeof(double));
        }
        for (int j = h0; j < h1; j++) { const int h = eh_list[j]; thsum[h] += dacc[h]; area[h] = ts->tmpA[j - h0]; }
        ts->dV += dV; ts->dDA += dDA; ts->dSphi += dSphi;
        if (aa != 0.0) { for (int j = h0; j < h1; j++) Vh[eh_list[j]] += dvh[eh_list[j]]; ts->dR2 += dR2; }
        ts->n_acc++;
    } else q[e] = qold;
    for (int j = h0; j < h1; j++) dacc[eh_list[j]] = 0.0;
    if (aa != 0.0) for (int j = h0; j < h1; j++) dvh[eh_list[j]] = 0.0;
}

/* Global dilation q -> s q (all edges), ln s uniform in [-seps, seps].
 * Validity (Gram positivity, q>0) is scale invariant; V -> s^2 V, A -> s A, deficits unchanged.
 * Uniform measure prod dq gives Jacobian s^NE, so accept with min(1, s^NE exp(-dI)). Exact MH. */
static double seps = 0.01;
static long n_sacc = 0, n_stry = 0;
static void scale_move(void) {
    n_stry++;
    const double t = seps * (2.0 * urand() - 1.0);
    const double s = exp(t);
    /* W_ij = V_ij/q_ij scales as s, so S_phi -> s S_phi */
    const double dI = lam * (s * s - 1.0) * run_sumV - kk * (s - 1.0) * run_sumDA + (s - 1.0) * run_Sphi;
    const double logacc = (double)NE * t - dI;
    if (logacc >= 0.0 || urand() < exp(logacc)) {
        const double s2 = s * s;
        for (int e = 0; e < NE; e++) q[e] *= s;
        for (int i = 0; i < NS; i++) V[i] *= s2;
        for (int h = 0; h < NT; h++) area[h] *= s;
        if (aa != 0.0) for (int h = 0; h < NT; h++) Vh[h] *= s2;   /* delta^2 A^2 / V_h is scale invariant */
        run_sumV *= s2; run_sumDA *= s; run_Sphi *= s;
        n_sacc++;
    }
}

/* ---------------- self tests ---------------- */
/* generalized cross product in R^4: vector orthogonal to a,b,c */
static void cross4(const double *a, const double *b, const double *c, double *n) {
    double m[3][4];
    for (int i = 0; i < 4; i++) { m[0][i] = a[i]; m[1][i] = b[i]; m[2][i] = c[i]; }
    for (int col = 0; col < 4; col++) {
        int cc[3], k = 0;
        for (int j = 0; j < 4; j++) if (j != col) cc[k++] = j;
        double d = m[0][cc[0]] * (m[1][cc[1]] * m[2][cc[2]] - m[1][cc[2]] * m[2][cc[1]])
                 - m[0][cc[1]] * (m[1][cc[0]] * m[2][cc[2]] - m[1][cc[2]] * m[2][cc[0]])
                 + m[0][cc[2]] * (m[1][cc[0]] * m[2][cc[1]] - m[1][cc[1]] * m[2][cc[0]]);
        n[col] = (col % 2 ? -d : d);
    }
}
/* independent dihedral angles from an explicit embedding */
static void angles_by_embedding(double X[5][4], double *ang, double *vol) {
    double nrm[5][4];
    for (int p = 0; p < 5; p++) {
        int f[4], n = 0;
        for (int a = 0; a < 5; a++) if (a != p) f[n++] = a;
        double u[3][4];
        for (int j = 0; j < 3; j++) for (int i = 0; i < 4; i++) u[j][i] = X[f[j+1]][i] - X[f[0]][i];
        cross4(u[0], u[1], u[2], nrm[p]);
        double dot = 0;
        for (int i = 0; i < 4; i++) dot += nrm[p][i] * (X[p][i] - X[f[0]][i]);
        if (dot > 0) for (int i = 0; i < 4; i++) nrm[p][i] = -nrm[p][i];   /* make outward */
    }
    for (int k = 0; k < 10; k++) {
        int p = PA[k], r = PB[k];
        double d = 0, np = 0, nr = 0;
        for (int i = 0; i < 4; i++) { d += nrm[p][i] * nrm[r][i]; np += nrm[p][i] * nrm[p][i]; nr += nrm[r][i] * nrm[r][i]; }
        double c = -d / sqrt(np * nr);
        if (c > 1) c = 1;
        if (c < -1) c = -1;
        ang[k] = acos(c);
    }
    double B[4][4];
    for (int a = 0; a < 4; a++) for (int i = 0; i < 4; i++) B[a][i] = X[a+1][i] - X[0][i];
    /* det via cofactor on first row using cross4 of rows 1..3 */
    double n[4]; cross4(B[1], B[2], B[3], n);
    double det = 0; for (int i = 0; i < 4; i++) det += B[0][i] * n[i];
    *vol = fabs(det) / 24.0;
}

static void selftest(void) {
    if (!build_lattice(1)) { printf("BUILD FAILED\n"); exit(1); }
    double sv, sda, sa, md;
    int bad = recompute_all(&sv, &sda, &sa, &md);
    double vmin = 1e300, vmax = -1e300;
    for (int s = 0; s < NS; s++) { if (V[s] < vmin) vmin = V[s]; if (V[s] > vmax) vmax = V[s]; }
    printf("FLAT invalid=%d  max|delta_h|=%.3e  sumV=%.15g (L^4=%d)  V in [%.15g, %.15g] (1/24=%.15g)  sum delta*A=%.3e\n",
           bad, md, sv, N0, vmin, vmax, 1.0 / 24.0, sda);
    /* regular simplex */
    double qq[10], ang[10], vol;
    for (int k = 0; k < 10; k++) qq[k] = 1.0;
    simplex_geom(qq, &vol, ang);
    double dmax = 0;
    for (int k = 0; k < 10; k++) if (fabs(ang[k] - acos(0.25)) > dmax) dmax = fabs(ang[k] - acos(0.25));
    printf("REGULAR theta=%.15f  arccos(1/4)=%.15f  max|diff|=%.3e  V=%.15f (sqrt5/96=%.15f)\n",
           ang[0], acos(0.25), dmax, vol, sqrt(5.0) / 96.0);
    /* random simplices: Gram-dual method vs explicit embedding; Schlaefli identity */
    seed_rng(12345);
    double maxdev = 0, maxvdev = 0, maxschl = 0;
    int ntest = 0;
    while (ntest < 2000) {
        double X[5][4];
        for (int a = 0; a < 5; a++) for (int i = 0; i < 4; i++) X[a][i] = 4.0 * urand() - 2.0;
        for (int k = 0; k < 10; k++) {
            double d = 0; for (int i = 0; i < 4; i++) { double t = X[PA[k]][i] - X[PB[k]][i]; d += t * t; }
            qq[k] = d;
        }
        double a1[10], a2[10], v1, v2;
        if (!simplex_geom(qq, &v1, a1)) continue;
        angles_by_embedding(X, a2, &v2);
        for (int k = 0; k < 10; k++) if (fabs(a1[k] - a2[k]) > maxdev) maxdev = fabs(a1[k] - a2[k]);
        if (fabs(v1 - v2) / v2 > maxvdev) maxvdev = fabs(v1 - v2) / v2;
        /* Schlaefli: sum_h A_h d(theta_h)/dq_e = 0 for each edge e */
        double Ah[10];
        for (int k = 0; k < 10; k++) {
            int a = TRI_OF_PAIR[k][0], b = TRI_OF_PAIR[k][1], c = TRI_OF_PAIR[k][2];
            Ah[k] = tri_area(qq[PIDX[a][b]], qq[PIDX[a][c]], qq[PIDX[b][c]]);
        }
        for (int e = 0; e < 10; e++) {
            double hstep = 1e-6 * qq[e], qp[10], qm[10], ap[10], am[10], vv;
            memcpy(qp, qq, sizeof qp); memcpy(qm, qq, sizeof qm);
            qp[e] += hstep; qm[e] -= hstep;
            if (!simplex_geom(qp, &vv, ap) || !simplex_geom(qm, &vv, am)) continue;
            double sch = 0, scale = 0;
            for (int k = 0; k < 10; k++) { double der = (ap[k] - am[k]) / (2 * hstep); sch += Ah[k] * der; scale += fabs(Ah[k] * der); }
            if (scale > 0 && fabs(sch) / scale > maxschl) maxschl = fabs(sch) / scale;
        }
        ntest++;
    }
    printf("RANDOM %d simplices: max|theta_gram - theta_embed|=%.3e  max rel|V diff|=%.3e  max Schlaefli residual (rel)=%.3e\n",
           ntest, maxdev, maxvdev, maxschl);
    /* incremental bookkeeping check: Metropolis at k=0.05 then compare with recompute */
    kk = 0.05; lam = 1.0; eps = 0.5;
    recompute_all(&run_sumV, &run_sumDA, NULL, NULL);
    run_sumR2 = full_sumR2;
    for (int sw = 0; sw < 3; sw++) for (int i = 0; i < NE; i++) update_edge((int)(next_u64() % (uint64_t)NE));
    double *V2 = malloc(sizeof(double) * NS), *ts2 = malloc(sizeof(double) * NT), *A2 = malloc(sizeof(double) * NT);
    memcpy(V2, V, sizeof(double) * NS); memcpy(ts2, thsum, sizeof(double) * NT); memcpy(A2, area, sizeof(double) * NT);
    double incV = run_sumV, incDA = run_sumDA, incR2 = run_sumR2;
    double *Vh2 = malloc(sizeof(double) * NT); memcpy(Vh2, Vh, sizeof(double) * NT);
    recompute_all(&sv, &sda, &sa, &md);
    if (aa != 0.0) {
        double dvhm = 0; for (int t = 0; t < NT; t++) if (fabs(Vh2[t] - Vh[t]) > dvhm) dvhm = fabs(Vh2[t] - Vh[t]);
        printf("R2TERM a=%g after 3 sweeps: sumR2 full=%.10g inc-full=%.2e max|dVh|=%.2e\n", aa, full_sumR2, incR2 - full_sumR2, dvhm);
    }
    free(Vh2);
    double dv = 0, dt = 0, dA = 0;
    for (int s = 0; s < NS; s++) if (fabs(V2[s] - V[s]) > dv) dv = fabs(V2[s] - V[s]);
    for (int t = 0; t < NT; t++) { if (fabs(ts2[t] - thsum[t]) > dt) dt = fabs(ts2[t] - thsum[t]); if (fabs(A2[t] - area[t]) > dA) dA = fabs(A2[t] - area[t]); }
    printf("INCREMENTAL after 3 sweeps (k=0.05, acc=%.3f): max|dV|=%.2e max|dthsum|=%.2e max|dA|=%.2e  sumV inc-full=%.2e  sumDA inc-full=%.2e\n",
           (double)n_acc / n_try, dv, dt, dA, incV - sv, incDA - sda);
    free(V2); free(ts2); free(A2);
}

static double wall_sec(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static atomic_int bar_cnt, bar_gen, par_quit;
static void par_barrier(void) {
    const int g = atomic_load_explicit(&bar_gen, memory_order_acquire);
    if (atomic_fetch_add_explicit(&bar_cnt, 1, memory_order_acq_rel) == nthr - 1) {
        atomic_store_explicit(&bar_cnt, 0, memory_order_relaxed);
        atomic_fetch_add_explicit(&bar_gen, 1, memory_order_acq_rel);
    } else {
        int spins = 0;
        while (atomic_load_explicit(&bar_gen, memory_order_acquire) == g) if (++spins > 4096) { sched_yield(); spins = 0; }
    }
}
static int setup_blocks(void) {     /* blocks of side b = L/4 in a 2^4 checkerboard; returns 1 if every color passes the check */
    blk_b = L / 4; const int nb1 = 4; nblk = nb1 * nb1 * nb1 * nb1;
    blk_off = calloc(nblk + 1, sizeof(int)); blk_list = malloc(sizeof(int) * NE);
    int n = 0;
    for (int B = 0; B < nblk; B++) {
        const int b0 = B % 4, b1 = (B / 4) % 4, b2 = (B / 16) % 4, b3 = B / 64;
        blk_off[B] = n;
        for (int x3 = b3 * blk_b; x3 < (b3 + 1) * blk_b; x3++) for (int x2 = b2 * blk_b; x2 < (b2 + 1) * blk_b; x2++)
        for (int x1 = b1 * blk_b; x1 < (b1 + 1) * blk_b; x1++) for (int x0 = b0 * blk_b; x0 < (b0 + 1) * blk_b; x0++) {
            const int v = x0 + L * (x1 + L * (x2 + L * x3));
            for (int t = 0; t < 15; t++) blk_list[n++] = 15 * v + t;
        }
    }
    blk_off[nblk] = n;
    if (n != NE) return 0;
    ncol = 16; cblk_off = calloc(ncol + 1, sizeof(int)); cblk_list = malloc(sizeof(int) * nblk);
    for (int c = 0, m = 0; c < 16; c++) { cblk_off[c] = m;
        for (int B = 0; B < nblk; B++) { const int pc = (B % 2) | ((B / 4) % 2) << 1 | ((B / 16) % 2) << 2 | ((B / 64) % 2) << 3; if (pc == c) cblk_list[m++] = B; }
        cblk_off[c + 1] = m; }
    int *stamp = malloc(sizeof(int) * NT); for (int t = 0; t < NT; t++) stamp[t] = -1;
    int ok = 1;
    for (int c = 0; c < ncol && ok; c++) {
        for (int i = cblk_off[c]; i < cblk_off[c + 1] && ok; i++) { const int B = cblk_list[i];
            for (int j = blk_off[B]; j < blk_off[B + 1] && ok; j++) { const int e = blk_list[j];
                for (int u = eh_off[e]; u < eh_off[e + 1]; u++) { const int h = eh_list[u];
                    if (stamp[h] >= 0 && stamp[h] != c * 1000000 + B) { ok = 0; break; }    /* another block of this color */
                    stamp[h] = c * 1000000 + B; } } }
        for (int i = cblk_off[c]; i < cblk_off[c + 1]; i++) { const int B = cblk_list[i];   /* clear this color's stamps */
            for (int j = blk_off[B]; j < blk_off[B + 1]; j++) { const int e = blk_list[j]; for (int u = eh_off[e]; u < eh_off[e + 1]; u++) stamp[eh_list[u]] = -1; } }
    }
    free(stamp);
    return ok;
}
static int setup_colors(void) {
    colp = 0;
    for (int p = 3; p <= L; p++) if (L % p == 0) { colp = p; break; }
    ncol = 15 * colp * colp * colp * colp;
    col_off = calloc(ncol + 1, sizeof(int)); col_list = malloc(sizeof(int) * NE);
    int *cl = malloc(sizeof(int) * NE);
    for (int e = 0; e < NE; e++) {
        int c[4]; vcoord(e / 15, c);
        cl[e] = e % 15 + 15 * ((c[0] % colp) + colp * ((c[1] % colp) + colp * ((c[2] % colp) + colp * (c[3] % colp))));
        col_off[cl[e] + 1]++;
    }
    for (int c = 0; c < ncol; c++) col_off[c + 1] += col_off[c];
    int *fill = malloc(sizeof(int) * ncol); memcpy(fill, col_off, sizeof(int) * ncol);
    for (int e = 0; e < NE; e++) col_list[fill[cl[e]]++] = e;
    free(fill); free(cl);
    int *stamp = malloc(sizeof(int) * NT); for (int t = 0; t < NT; t++) stamp[t] = -1;
    int ok = 1;
    for (int c = 0; c < ncol && ok; c++)
        for (int i = col_off[c]; i < col_off[c + 1] && ok; i++) { const int e = col_list[i];
            for (int j = eh_off[e]; j < eh_off[e + 1]; j++) { const int h = eh_list[j]; if (stamp[h] == c) { ok = 0; break; } stamp[h] = c; } }
    free(stamp);
    return ok;
}
static void par_colors(int tid) {
    tstate_t *ts = &TS[tid];
    if (blk_b) {       /* block mode: the threads take this color's blocks in turn; each block in local order */
        for (int c = 0; c < ncol; c++) {
            for (int i = cblk_off[c] + tid; i < cblk_off[c + 1]; i += nthr) { const int B = cblk_list[i];
                for (int j = blk_off[B]; j < blk_off[B + 1]; j++) update_edge_par(blk_list[j], ts); }
            par_barrier();
        }
        return;
    }
    for (int c = 0; c < ncol; c++) {
        const int b0 = col_off[c], n = col_off[c + 1] - b0;
        const int i0 = b0 + (int)((long)n * tid / nthr), i1 = b0 + (int)((long)n * (tid + 1) / nthr);
        for (int i = i0; i < i1; i++) update_edge_par(col_list[i], ts);
        par_barrier();
    }
}
#define RCH 4096
static int *ht_off = NULL, *ht_list = NULL, nch_s = 0, nch_t = 0, par_phase = 0, *cs_bad = NULL;
static double *cs_v = NULL, *ct_da = NULL, *ct_a = NULL, *ct_md = NULL, *ct_r2 = NULL;
static void recompute_part(int tid) {
    for (int ch = tid; ch < nch_s; ch += nthr) {
        const int s0 = ch * RCH, s1 = s0 + RCH < NS ? s0 + RCH : NS;
        double qq[10], sv = 0.0; int bad = 0;
        for (int s = s0; s < s1; s++) { for (int k = 0; k < 10; k++) qq[k] = q[sedge[s*10 + k]]; if (!simplex_geom(qq, &V[s], &th[s*10])) bad++; sv += V[s]; }
        cs_v[ch] = sv; cs_bad[ch] = bad;
    }
    par_barrier();
    for (int ch = tid; ch < nch_t; ch += nthr) {
        const int t0 = ch * RCH, t1 = t0 + RCH < NT ? t0 + RCH : NT;
        double sda = 0.0, sa = 0.0, md = 0.0, r2 = 0.0;
        for (int t = t0; t < t1; t++) {
            double tsum = 0.0, vh = 0.0;
            for (int j = ht_off[t]; j < ht_off[t + 1]; j++) { const int sk = ht_list[j]; tsum += th[sk]; vh += 0.1 * V[sk / 10]; }
            thsum[t] = tsum;
            area[t] = tri_area(q[tedge[3*t]], q[tedge[3*t+1]], q[tedge[3*t+2]]);
            const double d = TWOPI - tsum;
            sda += d * area[t]; sa += area[t]; if (fabs(d) > md) md = fabs(d);
            if (aa != 0.0) { Vh[t] = vh; if (vh > 0.0) r2 += d * d * area[t] * area[t] / vh; }
        }
        ct_da[ch] = sda; ct_a[ch] = sa; ct_md[ch] = md; ct_r2[ch] = r2;
    }
    par_barrier();
}
static void *par_worker(void *arg) {
    const int tid = (int)(intptr_t)arg;
    for (;;) {
        par_barrier();
        if (atomic_load(&par_quit)) return NULL;
        if (par_phase == 1) par_colors(tid); else recompute_part(tid);
    }
}
static pthread_t *par_threads = NULL;
static void par_start(void) {
    if (L % 4 == 0 && L / 4 >= 3 && !(getenv("PARFINE") && atoi(getenv("PARFINE")))) {
        if (!setup_blocks()) { fprintf(stderr, "PAR: block check failed (b=%d)\n", blk_b); exit(1); }
    } else { blk_b = 0; if (!setup_colors()) { fprintf(stderr, "PAR: coloring check failed (p=%d)\n", colp); exit(1); } }
    TS = calloc(nthr, sizeof(tstate_t));
    for (int i = 0; i < nthr; i++) TS[i].tmpA = malloc(sizeof(double) * max_eh);
    ht_off = calloc(NT + 1, sizeof(int)); ht_list = malloc(sizeof(int) * NS * 10);
    for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) ht_off[shinge[s*10 + k] + 1]++;
    for (int t = 0; t < NT; t++) ht_off[t + 1] += ht_off[t];
    int *fill = malloc(sizeof(int) * NT); memcpy(fill, ht_off, sizeof(int) * NT);
    for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) ht_list[fill[shinge[s*10 + k]]++] = s*10 + k;   /* ascending (s, k), as the serial sum */
    free(fill);
    nch_s = (NS + RCH - 1) / RCH; nch_t = (NT + RCH - 1) / RCH;
    cs_v = malloc(sizeof(double) * nch_s); cs_bad = malloc(sizeof(int) * nch_s);
    ct_da = malloc(sizeof(double) * nch_t); ct_a = malloc(sizeof(double) * nch_t); ct_md = malloc(sizeof(double) * nch_t); ct_r2 = malloc(sizeof(double) * nch_t);
    atomic_store(&bar_cnt, 0); atomic_store(&bar_gen, 0); atomic_store(&par_quit, 0);
    par_threads = malloc(sizeof(pthread_t) * nthr);
    for (int i = 1; i < nthr; i++) pthread_create(&par_threads[i], NULL, par_worker, (void *)(intptr_t)i);
}
static void par_stop(void) {
    atomic_store(&par_quit, 1); par_barrier();
    for (int i = 1; i < nthr; i++) pthread_join(par_threads[i], NULL);
}
static void par_sweep(void) {
    sweep_key = next_u64();
    par_phase = 1; par_barrier(); par_colors(0);
    for (int i = 0; i < nthr; i++) {
        tstate_t *ts = &TS[i];
        n_acc += ts->n_acc; n_inval += ts->n_inval; n_neg += ts->n_neg; n_try += ts->n_try;
        run_sumV += ts->dV; run_sumDA += ts->dDA; run_Sphi += ts->dSphi; run_sumR2 += ts->dR2;
        ts->n_acc = ts->n_inval = ts->n_neg = ts->n_try = 0; ts->dV = ts->dDA = ts->dSphi = ts->dR2 = 0.0;
    }
}
static int recompute_all_par(double *sumV, double *sumDA, double *sumA, double *maxabsdelta) {
    par_phase = 2; par_barrier(); recompute_part(0);
    double sv = 0.0, sda = 0.0, sa = 0.0, md = 0.0, r2 = 0.0; int bad = 0;
    for (int ch = 0; ch < nch_s; ch++) { sv += cs_v[ch]; bad += cs_bad[ch]; }
    for (int ch = 0; ch < nch_t; ch++) { sda += ct_da[ch]; sa += ct_a[ch]; if (ct_md[ch] > md) md = ct_md[ch]; r2 += ct_r2[ch]; }
    if (aa != 0.0) full_sumR2 = r2;
    if (sumV) *sumV = sv;
    if (sumDA) *sumDA = sda;
    if (sumA) *sumA = sa;
    if (maxabsdelta) *maxabsdelta = md;
    return bad;
}

/* ---------------- production run ---------------- */
static void run(double k, double lambda, int ntherm, int nmeas, uint64_t seed, const char *outfile, double eps0,
                const char *loadfile, const char *savefile, int nscale, int nphi_in, int freeze) {
    if (!build_lattice(0)) exit(1);
    kk = k; lam = lambda; eps = eps0;
    nphi = nphi_in;
    if (nphi > 16) { fprintf(stderr, "nphi <= 16\n"); exit(1); }
    build_scalar_tables();
    /* time slicing (axis 0 = time): freeze=1 -> edges along e_0 alone fixed; freeze=2 -> every edge whose
     * step set contains axis 0 fixed. Frozen edges keep their initial (flat) values; only free edges are
     * proposed. The global dilation would move frozen edges, so it is switched off when freeze != 0. */
    int *free_list = malloc(sizeof(int) * NE), nfree = 0;
    for (int e = 0; e < NE; e++) {
        int mask = e % 15 + 1;
        int frz = (freeze == 1 && mask == 1) || (freeze == 2 && (mask & 1));
        if (!frz) free_list[nfree++] = e;
    }
    if (freeze != 0) nscale = 0;
    seed_rng(seed);
    int loaded = 0;
    if (loadfile && strcmp(loadfile, "-")) {     /* hot start from a saved configuration */
        FILE *fl = fopen(loadfile, "rb");
        if (!fl) { perror("load"); exit(1); }
        int Lf = 0;
        if (fread(&Lf, sizeof(int), 1, fl) != 1 || Lf != L) { fprintf(stderr, "bad config file\n"); exit(1); }
        if (fread(&eps, sizeof(double), 1, fl) != 1 || fread(q, sizeof(double), NE, fl) != (size_t)NE) { fprintf(stderr, "short config\n"); exit(1); }
        fclose(fl); loaded = 1;
    }
    if (par_on) { if (freeze != 0) { fprintf(stderr, "PAR needs freeze = 0\n"); exit(1); } par_start(); }
    double sv, sda, sa, md;
    int bad = par_on ? recompute_all_par(&sv, &sda, &sa, &md) : recompute_all(&sv, &sda, &sa, &md);
    FILE *fo = fopen(outfile, "w");
    if (!fo) { perror("fopen"); exit(1); }
    fprintf(fo, "# L=%d k=%.6f lambda=%.6f ntherm=%d nmeas=%d seed=%llu N0=%d NE=%d NT=%d NS=%d\n",
            L, k, lambda, ntherm, nmeas, (unsigned long long)seed, N0, NE, NT, NS);
    fprintf(fo, "# %s start: invalid=%d max|delta|=%.3e sumV=%.10g\n", loaded ? "loaded" : "flat", bad, md, sv);
    fprintf(fo, "# nscale=%d global dilation moves per sweep\n", nscale);
    if (par_on && blk_b) fprintf(fo, "# edge order: block-parallel sweep (PAR=%d threads; %d blocks of side %d in 16 checkerboard colors, each block in local order)\n", nthr, nblk, blk_b);
    else if (par_on) fprintf(fo, "# edge order: colored parallel sweep (PAR=%d threads; colors = edge type x base vertex mod %d, %d colors)\n", nthr, colp, ncol);
    else fprintf(fo, "# edge order: %s\n", seq_on ? "sequential from a random offset (SEQ=1)" : "random with replacement");
    if (slices_on) fprintf(fo, "# slices: 4*L columns S_d(x) after the standard columns\n");
    if (lines_on) fprintf(fo, "# lines: then 1 + (L/2 + 1) columns: mean length of the closed lines along axis 3, and P(d) = <l(x) l(x + d)> over the three spatial axes, d = 0..L/2\n");
    fprintf(fo, "# nphi=%d scalar fields (heat bath each sweep, zero mode removed by mean subtraction)\n", nphi);
    fprintf(fo, "# freeze=%d nfree=%d (0 none, 1 e0-edges fixed, 2 all edges containing axis 0 fixed; dilation off if freeze)\n", freeze, nfree);
    if (aa != 0.0) fprintf(fo, "# a=%.8g curvature-squared term a*sum_h delta^2 A^2/V_h; start sumR2=%.6e; extra last column sumR2\n", aa, full_sumR2);
    fprintf(fo, "# cols: sweep phase(0=therm,1=meas) sumV sumDA R meanq acc eps sumA maxq minV q1 q2 q3 q4 drift invalfrac negfrac sacc"
                " mean_l cv_l std_delta mean_delta frac_neg_delta cv_V Sphi_per_field meanq_free\n");
    run_sumV = sv; run_sumDA = sda; run_Sphi = 0.0; run_sumR2 = full_sumR2;
    /* histograms over measurement sweeps: deficit angles in [-4pi, 2pi], l/<l> in [0, 5] */
    enum { NBD = 600, NBL = 250 };
    static long hist_d[NBD], hist_l[NBL];
    memset(hist_d, 0, sizeof hist_d); memset(hist_l, 0, sizeof hist_l);
    const double dlo = -2.0 * TWOPI, dhi = TWOPI;
    clock_t t0 = clock(); const double w0 = wall_sec();   /* PAR: MAXSEC and the reported seconds are wall time */
    double maxdrift = 0;
    const double maxsec = getenv("MAXSEC") ? atof(getenv("MAXSEC")) : 0.0;
    for (int sw = 0; sw < ntherm + nmeas; sw++) {
        long a0 = n_acc, t0n = n_try, i0 = n_inval, g0 = n_neg, sa0 = n_sacc, st0 = n_stry;
        if (nphi > 0) heatbath_sweep();
        if (par_on) par_sweep();
        else if (seq_on) { const uint64_t off = next_u64() % (uint64_t)nfree;   /* SEQ=1: edges in index order from a random start: memory-local */
            for (int i = 0; i < nfree; i++) update_edge(free_list[(off + (uint64_t)i) % (uint64_t)nfree]); }
        else for (int i = 0; i < nfree; i++) update_edge(free_list[next_u64() % (uint64_t)nfree]);
        for (int i = 0; i < nscale; i++) scale_move();
        double sacc = nscale > 0 ? (double)(n_sacc - sa0) / (double)(n_stry - st0) : 0.0;
        double acc = (double)(n_acc - a0) / (double)(n_try - t0n);
        double ifr = (double)(n_inval - i0) / (double)(n_try - t0n);
        double nfr = (double)(n_neg - g0) / (double)(n_try - t0n);
        double incV = run_sumV, incDA = run_sumDA, incS = run_Sphi, incR2 = run_sumR2;
        bad = par_on ? recompute_all_par(&sv, &sda, &sa, &md) : recompute_all(&sv, &sda, &sa, &md);
        if (bad) { fprintf(stderr, "invalid simplex after sweep %d\n", sw); exit(2); }
        double drift = fabs(incV - sv) / sv + fabs(incDA - sda) / (sa > 0 ? sa : 1);
        double Sphi = 0.0;
        if (nphi > 0) {
            compute_Vedge();
            Sphi = compute_Sphi();
            drift += fabs(incS - Sphi) / (Sphi > 0 ? Sphi : 1);
            run_Sphi = Sphi;
        }
        if (aa != 0.0) { drift += fabs(incR2 - full_sumR2) / (full_sumR2 > 0 ? full_sumR2 : 1); run_sumR2 = full_sumR2; }
        if (drift > maxdrift) maxdrift = drift;
        run_sumV = sv; run_sumDA = sda;
        /* local roughness */
        double sl = 0, sl2 = 0;
        for (int e = 0; e < NE; e++) { double l = sqrt(q[e]); sl += l; sl2 += q[e]; }
        const double mean_l = sl / NE, cv_l = sqrt(fmax(sl2 / NE - mean_l * mean_l, 0.0)) / mean_l;
        double sdl = 0, sdl2 = 0; long nneg = 0;
        for (int t = 0; t < NT; t++) { double d = TWOPI - thsum[t]; sdl += d; sdl2 += d * d; if (d < 0) nneg++; }
        const double mean_d = sdl / NT, std_d = sqrt(fmax(sdl2 / NT - mean_d * mean_d, 0.0)), fneg = (double)nneg / NT;
        double svv = 0, svv2 = 0;
        for (int s = 0; s < NS; s++) { svv += V[s]; svv2 += V[s] * V[s]; }
        const double cv_V = sqrt(fmax(svv2 / NS - (svv / NS) * (svv / NS), 0.0)) / (svv / NS);
        if (sw >= ntherm) {
            for (int t = 0; t < NT; t++) {
                int b = (int)((TWOPI - thsum[t] - dlo) / (dhi - dlo) * NBD);
                if (b < 0) b = 0;
                if (b >= NBD) b = NBD - 1;
                hist_d[b]++;
            }
            for (int e = 0; e < NE; e++) {
                int b = (int)(sqrt(q[e]) / mean_l / 5.0 * NBL);
                if (b >= NBL) b = NBL - 1;
                hist_l[b]++;
            }
        }
        double mq = 0, maxq = 0, qc[5] = {0}; int nc[5] = {0};
        for (int e = 0; e < NE; e++) {
            mq += q[e]; if (q[e] > maxq) maxq = q[e];
            int c = __builtin_popcount((e % 15) + 1); qc[c] += q[e]; nc[c]++;
        }
        mq /= NE;
        double mqf = 0; for (int i = 0; i < nfree; i++) mqf += q[free_list[i]];
        double minV = 1e300; for (int s = 0; s < NS; s++) if (V[s] < minV) minV = V[s];
        double R = 2.0 * sda / sv;
        int phase = sw < ntherm ? 0 : 1;
        fprintf(fo, "%d %d %.10g %.10g %.10g %.10g %.5f %.6g %.10g %.6g %.6g %.6g %.6g %.6g %.6g %.3e %.5f %.5f %.3f"
                    " %.6g %.6g %.6g %.6g %.6g %.6g %.8g %.6g\n",
                sw, phase, sv, sda, R, mq, acc, eps, sa, maxq, minV,
                qc[1] / nc[1], qc[2] / nc[2], qc[3] / nc[3], qc[4] / nc[4], drift, ifr, nfr, sacc,
                mean_l, cv_l, std_d, mean_d, fneg, cv_V, nphi > 0 ? Sphi / nphi : 0.0, mqf / nfree);
        if (aa != 0.0) { fseek(fo, -1, SEEK_CUR); fprintf(fo, " %.10g\n", full_sumR2); }
        if (slices_on) {   /* S_d(x) = sum of delta_h A_h over hinges whose base vertex has coordinate x in direction d */
            static double *S = NULL; if (!S) S = malloc(sizeof(double) * 4 * L);
            memset(S, 0, sizeof(double) * 4 * L);
            for (int t = 0; t < NT; t++) { int c[4]; vcoord(t / 50, c); double w = (TWOPI - thsum[t]) * area[t];
                for (int d = 0; d < 4; d++) S[d * L + c[d]] += w; }
            fseek(fo, -1, SEEK_CUR);
            for (int i = 0; i < 4 * L; i++) fprintf(fo, " %.10g", S[i]);
            fprintf(fo, "\n");
        }
        if (lines_on) {   /* closed lines along axis 3 (edge 15 v + 7 joins v to v + e_3): length l(x) at each of the L^3 spatial sites;
                             appended: mean l, then P(d) = mean over sites and the three spatial axes of l(x) l(x + d e_a), d = 0..L/2 */
            const int L3 = L * L * L; static double *ell = NULL; if (!ell) ell = malloc(sizeof(double) * L3);
            memset(ell, 0, sizeof(double) * L3);
            for (int v = 0; v < N0; v++) ell[v % L3] += sqrt(q[15 * v + 7]);
            double mean = 0.0; for (int x = 0; x < L3; x++) mean += ell[x]; mean /= L3;
            fseek(fo, -1, SEEK_CUR); fprintf(fo, " %.10g", mean);
            for (int d = 0; d <= L / 2; d++) {
                double p = 0.0;
                for (int x = 0; x < L3; x++) { const int c0 = x % L, c1 = (x / L) % L, c2 = x / (L * L);
                    p += ell[x] * (ell[(c0 + d) % L + L * (c1 + L * c2)] + ell[c0 + L * ((c1 + d) % L + L * c2)] + ell[c0 + L * (c1 + L * ((c2 + d) % L))]); }
                fprintf(fo, " %.10g", p / (3.0 * L3));
            }
            fprintf(fo, "\n");
        }
        if (phase == 0) {   /* tune step sizes toward 50% acceptance during thermalization only */
            eps *= exp(0.8 * (acc - 0.5));
            if (eps < 1e-6) eps = 1e-6;
            if (eps > 1e6) eps = 1e6;
            if (nscale > 0) { seps *= exp(0.3 * (sacc - 0.5)); if (seps > 0.5) seps = 0.5; if (seps < 1e-5) seps = 1e-5; }
        }
        if (!isfinite(sv) || !isfinite(sda)) { fprintf(stderr, "non-finite observables\n"); break; }
        if (maxsec > 0 && (par_on ? wall_sec() - w0 : (double)(clock() - t0) / CLOCKS_PER_SEC) > maxsec) break;   /* time-boxed chunk: stop after a whole sweep */
    }
    double secs = par_on ? wall_sec() - w0 : (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(fo, "# done: cpu_seconds=%.2f sec_per_sweep=%.5f max_drift=%.3e total_acc=%.5f scale_acc=%.4f seps=%.5f\n",
            secs, secs / (ntherm + nmeas), maxdrift, (double)n_acc / n_try,
            n_stry ? (double)n_sacc / n_stry : 0.0, seps);
    fprintf(fo, "# hist_delta lo=%.10g hi=%.10g nb=%d:", dlo, dhi, (int)NBD);
    for (int b = 0; b < NBD; b++) fprintf(fo, " %ld", hist_d[b]);
    fprintf(fo, "\n# hist_lrel lo=0 hi=5 nb=%d:", (int)NBL);
    for (int b = 0; b < NBL; b++) fprintf(fo, " %ld", hist_l[b]);
    fprintf(fo, "\n");
    fclose(fo);
    if (savefile && strcmp(savefile, "-")) {
        FILE *fs = fopen(savefile, "wb");
        if (!fs) { perror("save"); exit(1); }
        fwrite(&L, sizeof(int), 1, fs); fwrite(&eps, sizeof(double), 1, fs); fwrite(q, sizeof(double), NE, fs);
        fclose(fs);
    }
    if (par_on) par_stop();
    printf("done L=%d k=%.4f secs=%.1f\n", L, k, secs);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s selftest L | run L k lambda ntherm nmeas seed outfile [eps0] [loadcfg|-] [savecfg|-] [nscale] [nphi] [freeze] | phitest L nphi nsweeps seed [cfg]\n", argv[0]); return 1; }
    init_tables();
    { const char *ra = getenv("R2A"); if (ra && *ra) aa = atof(ra); }
    { const char *sq = getenv("SEQ"); if (sq && *sq) seq_on = atoi(sq); }
    { const char *pr = getenv("PAR"); if (pr && *pr && atoi(pr) > 0) { par_on = 1; nthr = atoi(pr); } }
    { const char *ln = getenv("LINES"); if (ln && *ln) lines_on = atoi(ln); }
    { const char *sl = getenv("SLICES"); if (sl && *sl) slices_on = atoi(sl); }
    if (aa < 0.0) { fprintf(stderr, "R2A must be >= 0\n"); return 1; }
    L = atoi(argv[2]);
    if (L < 3 || L > (strcmp(argv[1], "selftest") ? 32 : 8)) { fprintf(stderr, "need 3 <= L <= 8 (selftest) or 32 (run)\n"); return 1; }
    if (!strcmp(argv[1], "selftest")) { selftest(); return 0; }
    if (!strcmp(argv[1], "run") && argc >= 9) {
        double eps0 = argc >= 10 ? atof(argv[9]) : 0.5;
        const char *lf = argc >= 11 ? argv[10] : "-", *sf = argc >= 12 ? argv[11] : "-";
        int nscale = argc >= 13 ? atoi(argv[12]) : 0;
        int np = argc >= 14 ? atoi(argv[13]) : 0;
        int frz = argc >= 15 ? atoi(argv[14]) : 0;
        run(atof(argv[3]), atof(argv[4]), atoi(argv[5]), atoi(argv[6]), (uint64_t)strtoull(argv[7], NULL, 10), argv[8], eps0, lf, sf,
            nscale, np, frz);
        return 0;
    }
    if (!strcmp(argv[1], "phitest") && argc >= 6) {
        /* equipartition check on a FIXED geometry: <S_phi>/nphi must equal (N0-1)/2 */
        nphi = atoi(argv[3]);
        int nsw = atoi(argv[4]);
        seed_rng((uint64_t)strtoull(argv[5], NULL, 10));
        if (!build_lattice(0)) return 1;
        build_scalar_tables();
        if (argc >= 7 && strcmp(argv[6], "-")) {
            FILE *fl = fopen(argv[6], "rb"); int Lf; double ed;
            if (!fl || fread(&Lf, sizeof(int), 1, fl) != 1 || Lf != L || fread(&ed, sizeof(double), 1, fl) != 1 ||
                fread(q, sizeof(double), NE, fl) != (size_t)NE) { fprintf(stderr, "bad cfg\n"); return 1; }
            fclose(fl);
        }
        double sv, sda, sa, md;
        if (recompute_all(&sv, &sda, &sa, &md)) { fprintf(stderr, "invalid geometry\n"); return 1; }
        int nb = 20, per = nsw / nb, ntherm_phi = 200;
        double bm[64] = {0};
        for (int sw = 0; sw < ntherm_phi; sw++) heatbath_sweep();
        for (int b = 0; b < nb; b++) {
            double acc = 0;
            for (int i = 0; i < per; i++) { heatbath_sweep(); acc += run_Sphi / nphi; }
            bm[b] = acc / per;
        }
        double m = 0, v = 0;
        for (int b = 0; b < nb; b++) m += bm[b];
        m /= nb;
        for (int b = 0; b < nb; b++) v += (bm[b] - m) * (bm[b] - m);
        double err = sqrt(v / (nb - 1) / nb);
        printf("PHITEST L=%d nphi=%d geometry=%s sweeps=%d: <S_phi>/nphi = %.4f +- %.4f   expected (N0-1)/2 = %.1f   (R of geometry = %.4f)\n",
               L, nphi, argc >= 7 ? argv[6] : "flat", nb * per, m, err, (N0 - 1) / 2.0, 2 * sda / sv);
        return 0;
    }
    fprintf(stderr, "bad arguments\n");
    return 1;
}
