/* rk2.h -- the Regge single-edge Metropolis update and the per-hinge recomputation with every neighbor found by formula:
 * no index tables that grow with the lattice. The lattice is the Freudenthal triangulation of Z_L^4 of regge_r2_par.c, with the same
 * numbering: edge 15*site + (step set - 1), triangle 50*site + code, simplex 24*site + permutation. Everything an edge needs lies at
 * sites whose coordinates differ from its base site by -1, 0 or +1, so one list of the 81 neighbor sites serves every lookup, and
 * two small tables (27 KB and 19 KB, the same for every site) say which neighbor and which local object.
 * State kept: q (squared edge lengths) and per hinge the angle sum, the area and the hinge volume. Old volumes and angles of the
 * simplices around an edge are recomputed when needed instead of stored.
 * Common subset of C and the Metal Shading Language. The includer defines real, RC, RKN, RK_SQRT, RK_ACOS, RK_EXP, RK_TOL, DEV, THR,
 * KTAB, rk_uint, rk_uchar. */
#ifndef RK2_LAYOUT
#define RK2_LAYOUT
#define RK2_MAXH 112
#define RK2_OFF_NS 0          /* [15] simplices around an edge of each type */
#define RK2_OFF_NH 15         /* [15] hinges affected by an edge of each type */
#define RK2_OFF_ES 30         /* [15][24][3]  neighbor code of the simplex's base site, permutation, local pair that is the edge */
#define RK2_OFF_QE 1110       /* [15][24][10][2]  the simplex's ten edges: neighbor code, edge type */
#define RK2_OFF_HJ 8310       /* [15][24][10]  the simplex's ten hinges: position in the edge's hinge list */
#define RK2_OFF_EH 11910      /* [15][112][3]  the edge's hinges: neighbor code, triangle code, contains the edge */
#define RK2_OFF_TE 16950      /* [15][112][3][2]  the three edges of each hinge: neighbor code, edge type */
#define RK2_TAB_BYTES 27030
#define RK2_MAXHT 16
#define RK2_R_NHT 0           /* [50] simplices around a triangle of each code */
#define RK2_R_ENT 50          /* [50][16][3]  neighbor code of the simplex's base site, permutation, local angle index */
#define RK2_R_QE 2450         /* [50][16][10][2]  the simplex's ten edges relative to the triangle's base site */
#define RK2_R_TE 18450        /* [50][3][2]  the triangle's three edges */
#define RK2_R_BYTES 18750
#endif

KTAB int RKN(PIDX)[25] = { 0,0,1,2,3,  0,0,4,5,6,  1,4,0,7,8,  2,5,7,0,9,  3,6,8,9,0 };
KTAB int RKN(PA)[10] = { 0,0,0,0,1,1,1,2,2,3 };
KTAB int RKN(PB)[10] = { 1,2,3,4,2,3,4,3,4,4 };

static inline rk_uint RKN(h32)(rk_uint x) { x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16; return x; }
static inline real RKN(u01)(rk_uint k0, rk_uint k1, rk_uint e, rk_uint j) {
    rk_uint h = RKN(h32)(e * 2U + j + 0x9e3779b9U);
    h = RKN(h32)(h ^ k0); h = RKN(h32)(h + k1);
    return (real)(h >> 8) * RC(5.9604644775390625e-8);
}

static inline int RKN(simplex_geom)(THR const real *qq, THR real *vol, THR real *ang) {
    real G[16], Lm[16], Li[16], Gi[16], M[25], inv[5];
    for (int a = 1; a <= 4; a++) G[(a-1)*4 + a-1] = qq[RKN(PIDX)[a]];
    for (int a = 1; a <= 4; a++)
        for (int b = a + 1; b <= 4; b++) {
            real g = RC(0.5) * (qq[RKN(PIDX)[a]] + qq[RKN(PIDX)[b]] - qq[RKN(PIDX)[a*5 + b]]);
            G[(a-1)*4 + b-1] = g; G[(b-1)*4 + a-1] = g;
        }
    real prod = RC(1.0);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
            real s = G[i*4 + j];
            for (int k = 0; k < j; k++) s -= Lm[i*4 + k] * Lm[j*4 + k];
            if (i == j) {
                if (!(s > RK_TOL * G[i*4 + i])) return 0;
                Lm[i*4 + i] = RK_SQRT(s);
                prod *= Lm[i*4 + i];
            } else Lm[i*4 + j] = s / Lm[j*4 + j];
        }
    *vol = prod / RC(24.0);
    for (int i = 0; i < 4; i++) {
        Li[i*4 + i] = RC(1.0) / Lm[i*4 + i];
        for (int j = 0; j < i; j++) {
            real s = RC(0.0);
            for (int k = j; k < i; k++) s += Lm[i*4 + k] * Li[k*4 + j];
            Li[i*4 + j] = -s * Li[i*4 + i];
        }
        for (int j = i + 1; j < 4; j++) Li[i*4 + j] = RC(0.0);
    }
    for (int a = 0; a < 4; a++)
        for (int b = a; b < 4; b++) {
            real s = RC(0.0);
            for (int k = b; k < 4; k++) s += Li[k*4 + a] * Li[k*4 + b];
            Gi[a*4 + b] = s; Gi[b*4 + a] = s;
        }
    real tot = RC(0.0);
    for (int a = 0; a < 4; a++) {
        real cs = RC(0.0);
        for (int c = 0; c < 4; c++) { M[(a+1)*5 + c+1] = Gi[a*4 + c]; cs += Gi[a*4 + c]; }
        M[a+1] = -cs; M[(a+1)*5] = -cs;
        tot += cs;
    }
    M[0] = tot;
    for (int p = 0; p < 5; p++) inv[p] = RC(1.0) / RK_SQRT(M[p*5 + p]);
    for (int k = 0; k < 10; k++) {
        int p = RKN(PA)[k], r = RKN(PB)[k];
        real c = -M[p*5 + r] * inv[p] * inv[r];
        if (c > RC(1.0)) c = RC(1.0);
        if (c < RC(-1.0)) c = RC(-1.0);
        ang[k] = RK_ACOS(c);
    }
    return 1;
}

static inline real RKN(tri_area)(real q1, real q2, real q3) {
    real u = q1 + q2 - q3;
    real t = RC(4.0) * q1 * q2 - u * u;
    return t > RC(0.0) ? RC(0.25) * RK_SQRT(t) : RC(0.0);
}

/* nb[d0 + 3 d1 + 9 d2 + 27 d3] = index of the site at coordinates c_i + d_i - 1 (periodic). */
static inline void RKN(neigh)(int v, int Ls, THR int *nb) {
    int c0 = v % Ls; int r = v / Ls; int c1 = r % Ls; r = r / Ls; int c2 = r % Ls; int c3 = r / Ls;
    int x0[3], x1[3], x2[3], x3[3];
    for (int d = 0; d < 3; d++) {
        int a = c0 + d - 1; x0[d] = a < 0 ? a + Ls : (a >= Ls ? a - Ls : a);
        a = c1 + d - 1; x1[d] = (a < 0 ? a + Ls : (a >= Ls ? a - Ls : a)) * Ls;
        a = c2 + d - 1; x2[d] = (a < 0 ? a + Ls : (a >= Ls ? a - Ls : a)) * Ls * Ls;
        a = c3 + d - 1; x3[d] = (a < 0 ? a + Ls : (a >= Ls ? a - Ls : a)) * Ls * Ls * Ls;
    }
    int n = 0;
    for (int d3 = 0; d3 < 3; d3++) for (int d2 = 0; d2 < 3; d2++) for (int d1 = 0; d1 < 3; d1++) for (int d0 = 0; d0 < 3; d0++)
        nb[n++] = x0[d0] + x1[d1] + x2[d2] + x3[d3];
}

/* One Metropolis step on the edge of type t at site v. st[0..3]: accepted, invalid geometry, negative length, tried.
 * Edges updated at the same time must have disjoint hinge lists (the colors guarantee it). Returns the action change in *dIout. */
static inline void RKN(update_edge2)(int v, int t, int Ls, real eps, real lam, real kk, real aa, rk_uint k0, rk_uint k1,
        DEV real *q, DEV real *thsum, DEV real *area, DEV real *Vh, DEV const rk_uchar *tab, DEV int *st, THR real *dIout) {
    int nb[81]; int hh[RK2_MAXH];
    real dl[RK2_MAXH], dvl[RK2_MAXH], An[RK2_MAXH], qq[10], ao[10], an[10];
    *dIout = RC(0.0);
    st[3] += 1;
    int e = 15 * v + t;
    real qold = q[e];
    real qnew = qold + eps * (RC(2.0) * RKN(u01)(k0, k1, (rk_uint)e, 0U) - RC(1.0));
    if (qnew <= RC(0.0)) { st[2] += 1; return; }
    RKN(neigh)(v, Ls, nb);
    int ns = tab[RK2_OFF_NS + t], nh = tab[RK2_OFF_NH + t];
    for (int j = 0; j < nh; j++) { dl[j] = RC(0.0); dvl[j] = RC(0.0); }
    real dV = RC(0.0);
    for (int i = 0; i < ns; i++) {
        int bq = RK2_OFF_QE + (t * 24 + i) * 20;
        for (int k = 0; k < 10; k++) qq[k] = q[nb[tab[bq + 2*k]] * 15 + tab[bq + 2*k + 1]];
        real vo, vn;
        if (!RKN(simplex_geom)(qq, &vo, ao)) { st[1] += 1; return; }
        qq[tab[RK2_OFF_ES + (t * 24 + i) * 3 + 2]] = qnew;
        if (!RKN(simplex_geom)(qq, &vn, an)) { st[1] += 1; return; }
        dV += vn - vo;
        real dv10 = RC(0.1) * (vn - vo);
        int bh = RK2_OFF_HJ + (t * 24 + i) * 10;
        for (int k = 0; k < 10; k++) { int j = tab[bh + k]; dl[j] += an[k] - ao[k]; dvl[j] += dv10; }
    }
    real dDA = RC(0.0), dR2 = RC(0.0);
    for (int j = 0; j < nh; j++) {
        int be = RK2_OFF_EH + (t * RK2_MAXH + j) * 3;
        int h = nb[tab[be]] * 50 + tab[be + 1];
        hh[j] = h;
        real dold = RC(6.283185307179586476925286766559) - thsum[h];
        real dnew = dold - dl[j];
        real Aold = area[h];
        real Anew = Aold;
        if (tab[be + 2]) {
            int bt = RK2_OFF_TE + (t * RK2_MAXH + j) * 6;
            int e1 = nb[tab[bt]] * 15 + tab[bt + 1], e2 = nb[tab[bt + 2]] * 15 + tab[bt + 3], e3 = nb[tab[bt + 4]] * 15 + tab[bt + 5];
            Anew = RKN(tri_area)(e1 == e ? qnew : q[e1], e2 == e ? qnew : q[e2], e3 == e ? qnew : q[e3]);
        }
        An[j] = Anew;
        dDA += dnew * Anew - dold * Aold;
        real vo = Vh[h], vn = vo + dvl[j];
        dR2 += dnew * dnew * Anew * Anew / vn - dold * dold * Aold * Aold / vo;
    }
    real dI = lam * dV - kk * dDA + aa * dR2;
    *dIout = dI;
    int accept = (dI <= RC(0.0)) || (RKN(u01)(k0, k1, (rk_uint)e, 1U) < RK_EXP(-dI));
    if (accept) {
        q[e] = qnew;
        for (int j = 0; j < nh; j++) { int h = hh[j]; thsum[h] += dl[j]; area[h] = An[j]; Vh[h] += dvl[j]; }
        st[0] += 1;
    }
}

/* Angle sum, area and hinge volume of the triangle with the given code at site u, from the squared edge lengths alone.
 * Returns the number of invalid simplices met. */
static inline int RKN(recompute_hinge)(int u, int code, int Ls, DEV const real *q, DEV real *thsum, DEV real *area, DEV real *Vh, DEV const rk_uchar *tab2) {
    int nb[81]; real qq[10], an[10];
    RKN(neigh)(u, Ls, nb);
    int n = tab2[RK2_R_NHT + code], bad = 0;
    real tsum = RC(0.0), vh = RC(0.0);
    for (int i = 0; i < n; i++) {
        int bq = RK2_R_QE + (code * RK2_MAXHT + i) * 20;
        for (int k = 0; k < 10; k++) qq[k] = q[nb[tab2[bq + 2*k]] * 15 + tab2[bq + 2*k + 1]];
        real vol;
        if (!RKN(simplex_geom)(qq, &vol, an)) { bad += 1; continue; }
        tsum += an[tab2[RK2_R_ENT + (code * RK2_MAXHT + i) * 3 + 2]];
        vh += RC(0.1) * vol;
    }
    int h = 50 * u + code, bt = RK2_R_TE + code * 6;
    thsum[h] = tsum; Vh[h] = vh;
    area[h] = RKN(tri_area)(q[nb[tab2[bt]] * 15 + tab2[bt + 1]], q[nb[tab2[bt + 2]] * 15 + tab2[bt + 3]], q[nb[tab2[bt + 4]] * 15 + tab2[bt + 5]]);
    return bad;
}
