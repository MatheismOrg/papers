#!/usr/bin/env python3
"""Driver for regge.c: self-tests, k-scan on L=4 and L=3, blocking analysis, JSON output.

Usage:
  python3 driver.py selftest            -> runs/selftest.txt
  python3 driver.py scan PLAN.json      -> runs every job in the plan (2 at a time)
  python3 driver.py analyze             -> results.json from all runs/*.dat listed in the plans
"""
import json, os, subprocess, sys, glob, time
import numpy as np
from concurrent.futures import ThreadPoolExecutor

D = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(D, "regge")
RUNS = os.path.join(D, "runs")
COLS = ["sweep", "phase", "sumV", "sumDA", "R", "meanq", "acc", "eps", "sumA", "maxq", "minV",
        "q1", "q2", "q3", "q4", "drift", "invalfrac", "negfrac", "sacc",
        "mean_l", "cv_l", "std_delta", "mean_delta", "frac_neg_delta", "cv_V", "Sphi_per_field", "meanq_free"]


def build():
    subprocess.run(["gcc", "-O2", "-Wall", "-o", EXE, os.path.join(D, "regge.c"), "-lm"], check=True)


def selftest():
    out = []
    for L in (4, 3):
        r = subprocess.run([EXE, "selftest", str(L)], capture_output=True, text=True, check=True)
        out.append(r.stdout)
    txt = "\n".join(out)
    open(os.path.join(RUNS, "selftest.txt"), "w").write(txt)
    print(txt)


def run_job(job):
    tag = job["tag"]
    rdir = os.path.join(D, job.get("dir", "runs"))
    dat = os.path.join(rdir, tag + ".dat")
    cmd = [EXE, "run", str(job["L"]), repr(job["k"]), "1.0", str(job["ntherm"]), str(job["nmeas"]),
           str(job["seed"]), dat, str(job.get("eps0", 0.5)),
           job.get("load", "-"), job.get("save", os.path.join(rdir, tag + ".cfg")), str(job.get("nscale", 0)),
           str(job.get("nphi", 0)), str(job.get("freeze", 0))]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    return tag, r.returncode, time.time() - t0, r.stderr.strip()


def scan(planfile, workers=2):
    jobs = json.load(open(planfile))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for tag, rc, dt, err in ex.map(run_job, jobs):
            print(f"{tag}: rc={rc} wall={dt:.1f}s {err}", flush=True)


# ---------------- analysis ----------------
def load(dat):
    hdr = {}
    rows = []
    tail = {}
    for line in open(dat):
        if line.startswith("# L="):
            for tok in line[2:].split():
                k, v = tok.split("=")
                hdr[k] = float(v) if "." in v else int(v)
        elif line.startswith("# done:"):
            for tok in line[7:].split():
                k, v = tok.split("=")
                tail[k] = float(v)
        elif line.startswith("#"):
            if "start:" in line:
                hdr["start"] = line.split()[1]
            if line.startswith("# nscale="):
                hdr["nscale"] = int(line.split()[1].split("=")[1])
            if line.startswith("# nphi="):
                hdr["nphi"] = int(line.split()[1].split("=")[1])
            if line.startswith("# freeze="):
                hdr["freeze"] = int(line.split()[1].split("=")[1])
                hdr["nfree"] = int(line.split()[2].split("=")[1])
            if line.startswith("# hist_"):
                name = line.split()[1]
                lo = float(line.split()[2].split("=")[1]); hi = float(line.split()[3].split("=")[1])
                counts = [int(x) for x in line.split(":", 1)[1].split()]
                hdr[name] = {"lo": lo, "hi": hi, "counts": counts}
        elif line.strip():
            rows.append([float(x) for x in line.split()])
    a = np.array(rows)
    d = {c: a[:, i] for i, c in enumerate(COLS) if i < a.shape[1]}
    if "sacc" not in d:
        d["sacc"] = np.zeros(len(a))
    return hdr, d, tail


def blocking(x, min_blocks=16):
    x = np.asarray(x, float)
    n = len(x)
    curve = []
    b = 1
    while n // b >= min_blocks:
        nb = n // b
        bm = x[: nb * b].reshape(nb, b).mean(1)
        curve.append((b, float(bm.std(ddof=1) / np.sqrt(nb))))
        b *= 2
    errs = [e for _, e in curve]
    err = max(errs)
    naive = errs[0]
    plateau = len(errs) >= 3 and errs[-1] < 1.2 * errs[-2]
    tau = 0.5 * (err / naive) ** 2 if naive > 0 else float("nan")
    return float(x.mean()), err, tau, plateau, curve


def jackknife_ratio(num, den, nblocks=16, scale=2.0):
    """scale*<num>/<den> with jackknife error over nblocks contiguous blocks."""
    n = len(num) // nblocks * nblocks
    nb = num[:n].reshape(nblocks, -1).sum(1)
    db = den[:n].reshape(nblocks, -1).sum(1)
    full = scale * nb.sum() / db.sum()
    jk = np.array([scale * (nb.sum() - nb[i]) / (db.sum() - db[i]) for i in range(nblocks)])
    err = np.sqrt((nblocks - 1) / nblocks * ((jk - jk.mean()) ** 2).sum())
    return float(full), float(err)


def analyze_run(dat):
    hdr, d, tail = load(dat)
    L, k = int(hdr["L"]), float(hdr["k"])
    N0, NE = int(hdr["N0"]), int(hdr["NE"])
    meas = d["phase"] == 1
    m = {c: d[c][meas] for c in COLS if c in d}
    nphi, freeze = int(hdr.get("nphi", 0)), int(hdr.get("freeze", 0))
    res = {"L": L, "k": k, "start": hdr.get("start", "flat"), "nscale": hdr.get("nscale", 0), "ntherm": int(hdr["ntherm"]),
           "nmeas_done": int(meas.sum()), "nmeas_planned": int(hdr["nmeas"]),
           "cpu_seconds": tail.get("cpu_seconds"), "max_drift": tail.get("max_drift"),
           "nphi": nphi, "freeze": freeze, "nfree": int(hdr.get("nfree", NE))}
    obs = {
        "R": m["R"],
        "meanq": m["meanq"],
        "acc": m["acc"],
        "sumV_per_N0": m["sumV"] / N0,
        "sumDA_per_N0": m["sumDA"] / N0,
        "sumA_per_N0": m["sumA"] / N0,
        "maxq": m["maxq"],
        "minV": m["minV"],
        "q_by_type_1": m["q1"], "q_by_type_2": m["q2"], "q_by_type_3": m["q3"], "q_by_type_4": m["q4"],
        "invalid_reject_frac": m["invalfrac"],
        "scale_acc": m["sacc"],
        # scaling Ward identity of the uniform measure: <2*lambda*sumV - k*sumDA + S_phi> = N_edges (exact);
        # not applicable when edges are frozen (global rescaling is then not a change of variables of the model)
        "virial_ratio": ((2.0 * m["sumV"] - k * m["sumDA"] + nphi * m["Sphi_per_field"]) / NE if nphi > 0
                         else (2.0 * m["sumV"] - k * m["sumDA"]) / NE),
        "R_times_meanq": m["R"] * m["meanq"],              # dimensionless (invariant under q -> s q)
        "area_weighted_mean_deficit": m["sumDA"] / m["sumA"],
    }
    if freeze:
        obs.pop("virial_ratio")
    for c in ("mean_l", "cv_l", "std_delta", "mean_delta", "frac_neg_delta", "cv_V", "meanq_free"):
        if c in m:
            obs[c] = m[c]
    if nphi > 0:
        obs["Sphi_per_field"] = m["Sphi_per_field"]
    res["eps"] = float(m["eps"][0])
    for hname in ("hist_delta", "hist_lrel"):
        if hname in hdr:
            h = hdr[hname]
            cnt = np.array(h["counts"], float)
            edges = np.linspace(h["lo"], h["hi"], len(cnt) + 1)
            cdf = np.cumsum(cnt) / cnt.sum()
            qs = {}
            for p in (0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99):
                i = int(np.searchsorted(cdf, p))
                qs[f"p{int(p * 100):02d}"] = round(float(edges[i + 1]), 4)
            res[hname + "_quantiles"] = qs
            if hname == "hist_delta":
                centers = 0.5 * (edges[1:] + edges[:-1])
                res["frac_abs_delta_gt_1rad"] = float(cnt[np.abs(centers) > 1.0].sum() / cnt.sum())
    for name, x in obs.items():
        mean, err, tau, plat, curve = blocking(x)
        res[name] = {"mean": mean, "err": err, "tau_int": tau, "plateau": plat}
        if name in ("R", "sumV_per_N0"):
            res[name]["blocking_curve"] = curve
    rr, rre = jackknife_ratio(m["sumDA"], m["sumV"])
    res["R_ratio_of_means"] = {"mean": rr, "err": rre}
    # stationarity: first vs second half of the measurement run
    h = len(m["R"]) // 2
    st = {}
    for name in ("R", "sumV_per_N0", "meanq"):
        x = obs[name]
        a1 = blocking(x[:h], 8)
        a2 = blocking(x[h:], 8)
        st[name] = {"first_half": a1[0], "second_half": a2[0],
                    "z": (a2[0] - a1[0]) / np.hypot(a1[1], a2[1]) if a1[1] + a2[1] > 0 else 0.0}
    res["stationarity"] = st
    # quarter means of the measurement window (drift diagnostic) and R histogram (bimodality check)
    nq = len(m["R"]) // 4
    res["quarters"] = {name: [float(obs[name][i * nq:(i + 1) * nq].mean()) for i in range(4)]
                       for name in ("R", "sumV_per_N0", "meanq", "maxq")}
    Rq = res["quarters"]["R"]
    zR = st["R"]["z"]
    if Rq[0] > 0 and Rq[0] > Rq[1] > Rq[2] > Rq[3]:
        status = "relaxing: started at R>0, R decreases monotonically toward the smooth state"
    elif min(Rq) > 5 and Rq[3] > Rq[1]:
        status = "degenerate runaway continues: R large and rising, edges growing (no equilibrium)"
    elif min(Rq) > 0 and not (Rq[3] > Rq[2] > Rq[1] > Rq[0]):
        status = "positive-R state throughout, no monotonic trend within run"
    elif Rq[3] > 0 and Rq[3] > Rq[2] > Rq[1] and Rq[3] > 2 * max(abs(Rq[0]), 0.5):
        status = "runaway: R rises to large positive values, still rising at end (no equilibrium)"
    elif Rq[3] > 0 or (m["R"][-500:].mean() > 0):
        status = "left smooth state during run: R crossed to positive values"
    elif abs(zR) < 3:
        status = "stationary (smooth, R<0)"
    else:
        status = f"smooth (R<0) but slow drift/fluctuation between halves (z={zR:.1f})"
    res["status"] = status
    # first 100-sweep window (counted from sweep 0, thermalization included) whose mean R is > 0
    wR = d["R"][: len(d["R"]) // 100 * 100].reshape(-1, 100).mean(1)
    res["first_positive_R_window_sweep"] = int(np.argmax(wR > 0) * 100) if (wR > 0).any() else None
    res["first_R_gt_5_window_sweep"] = int(np.argmax(wR > 5) * 100) if (wR > 5).any() else None
    if k == 0.0:
        # exact at k=0: lambda*sumV ~ Gamma(N_edges/2, 1) for the uniform measure
        x = m["sumV"] / N0
        res["gamma_check_k0"] = {"mean_sumV_per_N0": float(x.mean()), "exact_mean": NE / 2 / N0,
                                 "sd_sumV_per_N0": float(x.std()), "exact_sd": float(np.sqrt(NE / 2) / N0)}
    cnt, edges = np.histogram(m["R"], bins=24)
    res["R_histogram"] = {"counts": cnt.tolist(), "edges": [round(float(e), 4) for e in edges]}
    # coarse history over the whole run (thermalization + measurement), windows of 100 sweeps
    w = 100
    nw = len(d["R"]) // w
    res["history_100sweep_windows"] = {
        "R": [round(float(v), 5) for v in d["R"][: nw * w].reshape(nw, w).mean(1)],
        "sumV_per_N0": [round(float(v), 4) for v in (d["sumV"][: nw * w] / N0).reshape(nw, w).mean(1)],
        "meanq": [round(float(v), 4) for v in d["meanq"][: nw * w].reshape(nw, w).mean(1)],
        "maxq": [round(float(v), 3) for v in d["maxq"][: nw * w].reshape(nw, w).max(1)],
    }
    return res


def analyze(plan_files, out_json):
    results = {"description": "4D quantum Regge calculus, periodic Freudenthal-triangulated Z_L^4, "
                              "I = lambda*sum V - k*sum delta*A, lambda=1, uniform measure prod dq_e, "
                              "single-edge Metropolis (eps tuned to ~50% during thermalization, then frozen) plus, "
                              "where nscale>0, nscale exact global-dilation MH moves per sweep (acceptance "
                              "min(1, s^N_edges exp(-dI))). R = 2 sum delta*A / sum V per sweep. "
                              "Errors from blocking (max over block sizes with >=16 blocks); "
                              "'plateau'=False means the blocking error had not converged (error likely underestimated). "
                              "virial_ratio = <2 lambda sumV - k sumDA>/N_edges must equal 1 in equilibrium (exact scaling identity).",
               "selftest": open(os.path.join(RUNS, "selftest.txt")).read() if os.path.exists(os.path.join(RUNS, "selftest.txt")) else None,
               "runs": []}
    for pf in plan_files:
        for job in json.load(open(pf)):
            dat = os.path.join(RUNS, job["tag"] + ".dat")
            if os.path.exists(dat):
                r = analyze_run(dat)
                r["tag"] = job["tag"]
                results["runs"].append(r)
    json.dump(results, open(out_json, "w"), indent=1)
    return results


def table(results):
    print(f"{'tag':>22} {'L':>2} {'k':>6} {'start':>6} {'nmeas':>5} {'R':>17} {'R(ratio)':>17} {'<q>':>15} "
          f"{'acc':>6} {'V/N0':>14} {'virial':>13} {'tauR':>6} {'zR':>6}")
    for r in results["runs"]:
        print(f"{r['tag']:>22} {r['L']:>2} {r['k']:6.3f} {r['start']:>6} {r['nmeas_done']:>5} "
              f"{r['R']['mean']:8.4f}±{r['R']['err']:.4f}{'' if r['R']['plateau'] else '*'} "
              f"{r['R_ratio_of_means']['mean']:8.4f}±{r['R_ratio_of_means']['err']:.4f} "
              f"{r['meanq']['mean']:7.3f}±{r['meanq']['err']:.3f} {r['acc']['mean']:6.3f} "
              f"{r['sumV_per_N0']['mean']:7.3f}±{r['sumV_per_N0']['err']:.3f} "
              f"{r['virial_ratio']['mean']:6.4f}±{r['virial_ratio']['err']:.4f} "
              f"{r['R']['tau_int']:6.1f} {r['stationarity']['R']['z']:6.2f}")


def tile_cfg(src, dst, Ldst):
    """Periodically tile a saved configuration of Z_Lsrc^4 onto Z_Ldst^4 (Ldst a multiple of Lsrc).
    The covering map is a local isomorphism of the triangulation, so every simplex keeps its geometry."""
    raw = open(src, "rb").read()
    Ls = int(np.frombuffer(raw[:4], dtype=np.int32)[0])
    eps = np.frombuffer(raw[4:12], dtype=np.float64)[0]
    qs = np.frombuffer(raw[12:], dtype=np.float64)
    assert len(qs) == 15 * Ls ** 4 and Ldst % Ls == 0
    x = np.indices((Ldst,) * 4).reshape(4, -1)          # x[0] fastest in vid ordering below
    vd = x[0] + Ldst * (x[1] + Ldst * (x[2] + Ldst * x[3]))
    xs = x % Ls
    vs = xs[0] + Ls * (xs[1] + Ls * (xs[2] + Ls * xs[3]))
    qd = np.empty(15 * Ldst ** 4)
    for m in range(15):
        qd[15 * vd + m] = qs[15 * vs + m]
    with open(dst, "wb") as f:
        f.write(np.int32(Ldst).tobytes()); f.write(np.float64(eps).tobytes()); f.write(qd.tobytes())
    print(f"tiled {src} (L={Ls}) -> {dst} (L={Ldst}), {len(qd)} edges")


def pre_escape_segment(dat, window=200):
    """For a run that left the smooth state: blocking analysis of the measurement segment before
    the first 200-sweep window whose mean R is > 0 (the window itself and the one before it are excluded)."""
    hdr, d, _ = load(dat)
    meas = d["phase"] == 1
    R = d["R"][meas]
    nw = len(R) // window
    wm = R[: nw * window].reshape(nw, window).mean(1)
    pos = np.nonzero(wm > 0)[0]
    if len(pos) == 0:
        return None
    end = max(pos[0] - 1, 1) * window
    mean, err, tau, plat, _ = blocking(R[:end])
    sweeps = d["sweep"][meas]
    return {"k": float(hdr["k"]), "segment_meas_sweeps": int(end), "escape_near_sweep": int(sweeps[0] + pos[0] * window),
            "R_mean": mean, "R_err": err, "tau_int": tau, "plateau": plat,
            "note": "smooth-branch value conditioned on not having escaped yet (biased selection; indicative only)"}


def final():
    plans = [os.path.join(D, p) for p in ("plan_L4.json", "plan_L3.json", "plan_hyst.json", "plan_pilot.json")]
    res = analyze([p for p in plans if os.path.exists(p)], os.path.join(D, "results.json"))
    extra = {}
    for tag in ("L4_k0.070", "L4_k0.065", "L3_k0.070", "L3_k0.080", "L3_k0.060"):
        f = os.path.join(RUNS, tag + ".dat")
        if os.path.exists(f):
            seg = pre_escape_segment(f)
            if seg:
                extra[tag] = seg
    res["pre_escape_segments"] = extra
    json.dump(res, open(os.path.join(D, "results.json"), "w"), indent=1)
    table(res)
    print()
    for r in res["runs"]:
        print(f"{r['tag']:>24}: {r['status']} | R quarters {[round(x, 3) for x in r['quarters']['R']]} "
              f"| maxq quarters {[round(x, 1) for x in r['quarters']['maxq']]}")
    print(json.dumps(extra, indent=1))


def followup():
    """Analysis of the follow-up runs (runs2/) -> results_followup.json (results.json is left untouched)."""
    R2 = os.path.join(D, "runs2")
    out = {"description": "Follow-up runs: (1) local roughness in the smooth phase, (2) L=6 near the edge of the smooth "
                          "phase incl. a degenerate start, (3) N_s=4 free massless scalars acting back "
                          "(S_phi = 1/2 sum_ij W_ij (phi_i-phi_j)^2, W_ij = V_ij/q_ij, V_ij = sum_s V_s/10; exact heat bath "
                          "per vertex each sweep, zero mode removed by subtracting each field's mean after every sweep; "
                          "dS_phi included in every edge Metropolis step and in the global dilation, S_phi -> s S_phi), "
                          "(4) time slicing: freeze=1 fixes edges along e_0 alone at q=1, freeze=2 fixes all edges whose "
                          "step set contains axis 0 at their flat values; the global dilation is switched off in frozen runs "
                          "(pure single-edge Metropolis on the free edges, exact detailed balance). "
                          "Same analysis conventions as results.json (blocking errors; plateau=False means error likely underestimated).",
           "checks": {}, "runs": [], "reference_round1": {}}
    for f in ("selftest_L6.txt", "phitest.txt"):
        p = os.path.join(R2, f)
        if os.path.exists(p):
            out["checks"][f] = open(p).read()
    out["checks"]["regression"] = ("new binary with nphi=0, freeze=0 reproduces round-1 run L4_k0.000 (seed 1006) "
                                   "bit-for-bit over 60 sweeps (max |diff| = 0.0)")
    p = os.path.join(R2, "test_phi.dat")
    if os.path.exists(p):
        t = analyze_run(p)
        out["checks"]["coupled_scalar_test_L4_k0.04"] = {
            "Sphi_per_field": t["Sphi_per_field"], "exact": 127.5, "virial_ratio_with_Sphi": t["virial_ratio"],
            "max_drift_incremental_vs_full": t["max_drift"]}
    for pf in ("plan_fu_A.json", "plan_fu_B.json", "plan_fu_C.json"):
        for job in json.load(open(os.path.join(D, pf))):
            dat = os.path.join(D, job.get("dir", "runs"), job["tag"] + ".dat")
            if os.path.exists(dat) and "# done" in open(dat).read()[-20000:]:
                r = analyze_run(dat)
                r["tag"] = job["tag"]
                seg = pre_escape_segment(dat)
                if seg and r["start"] == "flat":
                    r["pre_escape_segment"] = seg
                out["runs"].append(r)
    # round-1 reference values (N_s = 0, no freezing) from results.json
    r1 = json.load(open(os.path.join(D, "results.json")))
    for r in r1["runs"]:
        if r["tag"].startswith(("L3_k", "L4_k")) and r["start"] == "flat" or r["tag"] == "L4_k0.060_fromRough":
            out["reference_round1"][r["tag"]] = {
                "R": r["R"]["mean"], "R_err": r["R"]["err"], "meanq": r["meanq"]["mean"], "acc": r["acc"]["mean"],
                "R_times_meanq": r["R"]["mean"] * r["meanq"]["mean"], "sumV_per_N0": r["sumV_per_N0"]["mean"],
                "status": r["status"], "R_quarters": r["quarters"]["R"], "nmeas": r["nmeas_done"],
                "ntherm": r["ntherm"], "maxq_quarters": r["quarters"]["maxq"]}
            dat = os.path.join(RUNS, r["tag"] + ".dat")
            if os.path.exists(dat):
                _, dd, _ = load(dat)
                wR = dd["R"][: len(dd["R"]) // 100 * 100].reshape(-1, 100).mean(1)
                out["reference_round1"][r["tag"]]["first_positive_R_window_sweep"] = int(np.argmax(wR > 0) * 100) if (wR > 0).any() else None
    out["reference_round1"]["L4_k0.070_pre_escape"] = r1.get("pre_escape_segments", {}).get("L4_k0.070")
    json.dump(out, open(os.path.join(D, "results_followup.json"), "w"), indent=1)

    def fmt(r, key, nd=3):
        x = r.get(key)
        return f"{x['mean']:.{nd}f}±{x['err']:.{nd}f}" if x else "-"
    print("tag | k | status | R | R*<q> | <q> | acc | V/N0 | virial | Sphi/Ns | cv_l | std_delta | frac_neg | Rquarters | maxq quarters")
    for r in out["runs"]:
        print(f"{r['tag']} | {r['k']} | {r['status']} | {fmt(r, 'R')} | {fmt(r, 'R_times_meanq', 2)} | {fmt(r, 'meanq')} | "
              f"{r['acc']['mean']:.3f} | {fmt(r, 'sumV_per_N0')} | {fmt(r, 'virial_ratio', 4)} | {fmt(r, 'Sphi_per_field', 2)} | "
              f"{fmt(r, 'cv_l', 4)} | {fmt(r, 'std_delta', 3)} | {fmt(r, 'frac_neg_delta', 4)} | "
              f"{[round(x, 3) for x in r['quarters']['R']]} | {[round(x, 1) for x in r['quarters']['maxq']]}"
              + (f" | meanq_free {fmt(r, 'meanq_free')}" if r['freeze'] else "")
              + (f" | pre-escape R {r['pre_escape_segment']['R_mean']:.3f}±{r['pre_escape_segment']['R_err']:.3f} "
                 f"escape~{r['pre_escape_segment']['escape_near_sweep']}" if r.get('pre_escape_segment') else ""))
    print()
    for r in out["runs"]:
        if r["tag"].startswith("R4_"):
            print(r["tag"], "delta quantiles", r.get("hist_delta_quantiles"), "| l/<l> quantiles", r.get("hist_lrel_quantiles"),
                  "| frac|delta|>1rad", round(r.get("frac_abs_delta_gt_1rad", float('nan')), 4),
                  "| area-weighted mean deficit", fmt(r, "area_weighted_mean_deficit", 4), "| cv_V", fmt(r, "cv_V", 3))


if __name__ == "__main__":
    os.makedirs(RUNS, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else "help"
    if cmd == "final":
        final(); sys.exit(0)
    if cmd == "followup":
        followup(); sys.exit(0)
    if cmd == "build":
        build()
    elif cmd == "selftest":
        build(); selftest()
    elif cmd == "scan":
        scan(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    elif cmd == "tile":
        tile_cfg(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    elif cmd == "analyze":
        res = analyze(sys.argv[2:-1], sys.argv[-1])
        table(res)
    elif cmd == "one":
        r = analyze_run(sys.argv[2])
        r.pop("history_100sweep_windows")
        print(json.dumps(r, indent=1))
    else:
        print(__doc__)
