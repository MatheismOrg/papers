"""Build regge_r2.c from the unchanged regge.c by exact, unique text replacements.
Adds Hamber's optional curvature-squared term a * sum_h delta_h^2 A_h^2 / V_h,
V_h = sum_{s containing h} V_s / 10, with a read from the environment variable R2A.
With R2A unset or 0 every code path, RNG draw and output line matches regge.c."""
from pathlib import Path
import hashlib, sys
src = Path(sys.argv[1]).read_text()
out_path = Path(sys.argv[2])
R = []
def rep(old, new):
    n = src.count(old)
    if n != 1: raise SystemExit(f"expected 1 match, got {n}: {old[:70]!r}")
    R.append((old, new))
rep(" * Action (Hamber convention, a = 0):  I = lambda * sum_s V_s - k * sum_h delta_h A_h\n",
    " * Action (Hamber convention):  I = lambda * sum_s V_s - k * sum_h delta_h A_h + a * sum_h delta_h^2 A_h^2 / V_h\n"
    " * with V_h = sum_{s containing h} V_s / 10. The coupling a is read from the environment variable R2A\n"
    " * (unset or 0: the a = 0 engine, identical moves, random numbers and output).\n")
rep(" * Build: gcc -O2 -o regge regge.c -lm\n", " * Build: gcc -O2 -o regge_r2 regge_r2.c -lm\n")
rep("static double *tmpA;\n",
    "static double *tmpA;\n"
    "static double aa = 0.0;                 /* curvature-squared coupling a (env R2A) */\n"
    "static double *Vh = NULL, *dvh = NULL;  /* [NT] hinge volumes V_h and pending changes */\n"
    "static double run_sumR2 = 0.0, full_sumR2 = 0.0;   /* sum_h delta^2 A^2 / V_h */\n")
rep("dacc = calloc(NT, sizeof(double)); tmpA = malloc(sizeof(double) * max_eh);\n",
    "dacc = calloc(NT, sizeof(double)); tmpA = malloc(sizeof(double) * max_eh);\n"
    "    Vh = calloc(NT, sizeof(double)); dvh = calloc(NT, sizeof(double));\n")
rep("""    if (sumV) *sumV = sv;
    if (sumDA) *sumDA = sda;""",
"""    if (aa != 0.0) {
        memset(Vh, 0, sizeof(double) * NT);
        for (int s = 0; s < NS; s++) for (int k = 0; k < 10; k++) Vh[shinge[s*10 + k]] += 0.1 * V[s];
        double r2 = 0.0;
        for (int t = 0; t < NT; t++) { double d = TWOPI - thsum[t]; if (Vh[t] > 0.0) r2 += d * d * area[t] * area[t] / Vh[t]; }
        full_sumR2 = r2;
    }
    if (sumV) *sumV = sv;
    if (sumDA) *sumDA = sda;""")
rep("""        for (int k = 0; k < 10; k++) dacc[shinge[s*10 + k]] += tmpth[(i - o0) * 10 + k] - th[s*10 + k];
    }
    const int h0 = eh_off[e], h1 = eh_off[e+1];
    double dDA = 0.0;""",
"""        for (int k = 0; k < 10; k++) dacc[shinge[s*10 + k]] += tmpth[(i - o0) * 10 + k] - th[s*10 + k];
    }
    if (aa != 0.0)
        for (int i = o0; i < o1; i++) {
            const int s = es_list[i];
            const double dv10 = 0.1 * (tmpV[i - o0] - V[s]);
            for (int k = 0; k < 10; k++) dvh[shinge[s*10 + k]] += dv10;
        }
    const int h0 = eh_off[e], h1 = eh_off[e+1];
    double dDA = 0.0, dR2 = 0.0;""")
rep("""        dDA += dnew * Anew - dold * Aold;
    }""",
"""        dDA += dnew * Anew - dold * Aold;
        if (aa != 0.0) {
            const double vo = Vh[h], vn = vo + dvh[h];
            dR2 += dnew * dnew * Anew * Anew / vn - dold * dold * Aold * Aold / vo;
        }
    }""")
rep("    const double dI = lam * dV - kk * dDA + dSphi;\n",
    "    double dI = lam * dV - kk * dDA + dSphi;\n"
    "    if (aa != 0.0) dI += aa * dR2;\n")
rep("""        run_sumV += dV; run_sumDA += dDA; run_Sphi += dSphi;
        n_acc++;""",
"""        run_sumV += dV; run_sumDA += dDA; run_Sphi += dSphi;
        if (aa != 0.0) { for (int j = h0; j < h1; j++) Vh[eh_list[j]] += dvh[eh_list[j]]; run_sumR2 += dR2; }
        n_acc++;""")
rep("    for (int j = h0; j < h1; j++) dacc[eh_list[j]] = 0.0;\n}\n",
    "    for (int j = h0; j < h1; j++) dacc[eh_list[j]] = 0.0;\n"
    "    if (aa != 0.0) for (int j = h0; j < h1; j++) dvh[eh_list[j]] = 0.0;\n}\n")
rep("""        for (int h = 0; h < NT; h++) area[h] *= s;
        run_sumV *= s2; run_sumDA *= s; run_Sphi *= s;""",
"""        for (int h = 0; h < NT; h++) area[h] *= s;
        if (aa != 0.0) for (int h = 0; h < NT; h++) Vh[h] *= s2;   /* delta^2 A^2 / V_h is scale invariant */
        run_sumV *= s2; run_sumDA *= s; run_Sphi *= s;""")
# self test: incremental check of V_h and the a-term
rep("""    recompute_all(&run_sumV, &run_sumDA, NULL, NULL);
    for (int sw = 0; sw < 3; sw++)""",
"""    recompute_all(&run_sumV, &run_sumDA, NULL, NULL);
    run_sumR2 = full_sumR2;
    for (int sw = 0; sw < 3; sw++)""")
rep("""    double incV = run_sumV, incDA = run_sumDA;
    recompute_all(&sv, &sda, &sa, &md);""",
"""    double incV = run_sumV, incDA = run_sumDA, incR2 = run_sumR2;
    double *Vh2 = malloc(sizeof(double) * NT); memcpy(Vh2, Vh, sizeof(double) * NT);
    recompute_all(&sv, &sda, &sa, &md);
    if (aa != 0.0) {
        double dvhm = 0; for (int t = 0; t < NT; t++) if (fabs(Vh2[t] - Vh[t]) > dvhm) dvhm = fabs(Vh2[t] - Vh[t]);
        printf("R2TERM a=%g after 3 sweeps: sumR2 full=%.10g inc-full=%.2e max|dVh|=%.2e\\n", aa, full_sumR2, incR2 - full_sumR2, dvhm);
    }
    free(Vh2);""")
# production run
rep("""    fprintf(fo, "# cols: sweep phase""",
"""    if (aa != 0.0) fprintf(fo, "# a=%.8g curvature-squared term a*sum_h delta^2 A^2/V_h; start sumR2=%.6e; extra last column sumR2\\n", aa, full_sumR2);
    fprintf(fo, "# cols: sweep phase""")
rep("    run_sumV = sv; run_sumDA = sda; run_Sphi = 0.0;\n",
    "    run_sumV = sv; run_sumDA = sda; run_Sphi = 0.0; run_sumR2 = full_sumR2;\n")
rep("        double incV = run_sumV, incDA = run_sumDA, incS = run_Sphi;\n",
    "        double incV = run_sumV, incDA = run_sumDA, incS = run_Sphi, incR2 = run_sumR2;\n")
rep("""        if (drift > maxdrift) maxdrift = drift;""",
"""        if (aa != 0.0) { drift += fabs(incR2 - full_sumR2) / (full_sumR2 > 0 ? full_sumR2 : 1); run_sumR2 = full_sumR2; }
        if (drift > maxdrift) maxdrift = drift;""")
rep("""                mean_l, cv_l, std_d, mean_d, fneg, cv_V, nphi > 0 ? Sphi / nphi : 0.0, mqf / nfree);""",
"""                mean_l, cv_l, std_d, mean_d, fneg, cv_V, nphi > 0 ? Sphi / nphi : 0.0, mqf / nfree);
        if (aa != 0.0) { fseek(fo, -1, SEEK_CUR); fprintf(fo, " %.10g\\n", full_sumR2); }""")
rep("""    init_tables();
    L = atoi(argv[2]);""",
"""    init_tables();
    { const char *ra = getenv("R2A"); if (ra && *ra) aa = atof(ra); }
    if (aa < 0.0) { fprintf(stderr, "R2A must be >= 0\\n"); return 1; }
    L = atoi(argv[2]);""")
for old, new in R: src = src.replace(old, new)
out_path.write_text(src)
print("replacements", len(R), "sha256", hashlib.sha256(src.encode()).hexdigest())
