/* rk.h -- the single-edge Metropolis update of regge_r2_par.c (update_edge_par, simplex_geom, tri_area) written in the common subset
 * of C and the Metal Shading Language, so that one source runs (1) on the CPU in double precision, (2) on the CPU in single precision
 * and (3) on the graphics chip in single precision. Scalar fields are not included (nphi = 0).
 * Random numbers: a 32-bit counter hash of (sweep key, edge, draw) giving 24-bit uniforms, identical in all three builds.
 * The includer defines: real, RC(x) (a real constant), RKN(name) (function name), RK_SQRT, RK_ACOS, RK_EXP, RK_TOL, DEV, THR, KTAB. */

KTAB int RKN(PIDX)[25] = { 0,0,1,2,3,  0,0,4,5,6,  1,4,0,7,8,  2,5,7,0,9,  3,6,8,9,0 };   /* local pair (i,j) -> 0..9 */
KTAB int RKN(PA)[10] = { 0,0,0,0,1,1,1,2,2,3 };
KTAB int RKN(PB)[10] = { 1,2,3,4,2,3,4,3,4,4 };

static inline rk_uint RKN(h32)(rk_uint x) { x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16; return x; }
static inline real RKN(u01)(rk_uint k0, rk_uint k1, rk_uint e, rk_uint j) {
    rk_uint h = RKN(h32)(e * 2U + j + 0x9e3779b9U);
    h = RKN(h32)(h ^ k0); h = RKN(h32)(h + k1);
    return (real)(h >> 8) * RC(5.9604644775390625e-8);        /* 2^-24 */
}

/* qq[10]: squared lengths of the local pairs. Returns 0 if the Gram matrix is not positive definite. */
static inline int RKN(simplex_geom)(THR const real *qq, THR real *vol, DEV real *ang) {
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

/* One Metropolis step on edge e. slot: index of this edge within its color (private scratch space and counters).
 * Edges updated at the same time must have disjoint hinge lists (the colors guarantee it), so nothing below is shared.
 * stat[4*slot + 0..3]: accepted, invalid geometry, negative length, tried. tmpV[slot*32 + 31] returns the action change. */
static inline void RKN(update_edge)(int e, int slot, real eps, real lam, real kk, real aa, rk_uint k0, rk_uint k1, int max_eh,
        DEV real *q, DEV real *V, DEV real *th, DEV real *thsum, DEV real *area, DEV real *Vh, DEV real *dacc, DEV real *dvh,
        DEV const int *es_off, DEV const int *es_list, DEV const int *sedge, DEV const int *shinge,
        DEV const int *eh_off, DEV const int *eh_list, DEV const rk_uchar *eh_cont, DEV const int *tedge,
        DEV real *tmpV, DEV real *tmpth, DEV real *tmpA, DEV int *stat) {
    DEV real *tV = tmpV + slot * 32; DEV real *tth = tmpth + slot * 320; DEV real *tA = tmpA + slot * max_eh;
    DEV int *st = stat + slot * 4;
    st[3] += 1;
    real qold = q[e];
    real qnew = qold + eps * (RC(2.0) * RKN(u01)(k0, k1, (rk_uint)e, 0U) - RC(1.0));
    if (qnew <= RC(0.0)) { st[2] += 1; return; }
    q[e] = qnew;
    int o0 = es_off[e], o1 = es_off[e+1];
    real dV = RC(0.0);
    real qq[10];
    for (int i = o0; i < o1; i++) {
        int s = es_list[i];
        for (int k = 0; k < 10; k++) qq[k] = q[sedge[s*10 + k]];
        real v;
        if (!RKN(simplex_geom)(qq, &v, tth + (i - o0) * 10)) { q[e] = qold; st[1] += 1; return; }
        tV[i - o0] = v;
        dV += v - V[s];
    }
    for (int i = o0; i < o1; i++) {
        int s = es_list[i];
        for (int k = 0; k < 10; k++) dacc[shinge[s*10 + k]] += tth[(i - o0) * 10 + k] - th[s*10 + k];
    }
    if (aa != RC(0.0))
        for (int i = o0; i < o1; i++) {
            int s = es_list[i];
            real dv10 = RC(0.1) * (tV[i - o0] - V[s]);
            for (int k = 0; k < 10; k++) dvh[shinge[s*10 + k]] += dv10;
        }
    int h0 = eh_off[e], h1 = eh_off[e+1];
    real dDA = RC(0.0), dR2 = RC(0.0);
    for (int j = h0; j < h1; j++) {
        int h = eh_list[j];
        real dold = RC(6.283185307179586476925286766559) - thsum[h];
        real dnew = dold - dacc[h];
        real Aold = area[h];
        real Anew = eh_cont[j] ? RKN(tri_area)(q[tedge[3*h]], q[tedge[3*h+1]], q[tedge[3*h+2]]) : Aold;
        tA[j - h0] = Anew;
        dDA += dnew * Anew - dold * Aold;
        if (aa != RC(0.0)) {
            real vo = Vh[h], vn = vo + dvh[h];
            dR2 += dnew * dnew * Anew * Anew / vn - dold * dold * Aold * Aold / vo;
        }
    }
    real dI = lam * dV - kk * dDA;
    if (aa != RC(0.0)) dI += aa * dR2;
    tV[31] = dI;
    int accept = (dI <= RC(0.0)) || (RKN(u01)(k0, k1, (rk_uint)e, 1U) < RK_EXP(-dI));
    if (accept) {
        for (int i = o0; i < o1; i++) {
            int s = es_list[i];
            V[s] = tV[i - o0];
            for (int k = 0; k < 10; k++) th[s*10 + k] = tth[(i - o0) * 10 + k];
        }
        for (int j = h0; j < h1; j++) { int h = eh_list[j]; thsum[h] += dacc[h]; area[h] = tA[j - h0]; }
        if (aa != RC(0.0)) for (int j = h0; j < h1; j++) Vh[eh_list[j]] += dvh[eh_list[j]];
        st[0] += 1;
    } else q[e] = qold;
    for (int j = h0; j < h1; j++) dacc[eh_list[j]] = RC(0.0);
    if (aa != RC(0.0)) for (int j = h0; j < h1; j++) dvh[eh_list[j]] = RC(0.0);
}
