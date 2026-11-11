
# ------------------------------------------------------------------ SA74 T31 (Nov 9, Gary 08:57: "update and package all for a42 and code to verify claims in a zip")
# Plans A2 (16^4: k = 0.244 from the 0.248 and 0.252 configurations, second chains at 0.248 and 0.252 from the other coupling's configuration,
# 0.250; a = 0.02 on 10^4 and 12^4 at k = 0.45, 0.46), A3 (two 16^4 forks at 0.244), A4 (a = 0.02 at 0.47 on 10^4 and 12^4) and A5 (20^4 at
# 0.244, 0.248, 0.252 with the multi-core engine) finished on the Mac on Nov 7 by 19:33. Numbers from regge/r2/night: multihist_T30.json and
# multihist_skip50_T30.json (4^4 to 12^4, unchanged since T30), a1_results_T31.json (16^4 pooled per coupling by lineage, 20^4 single chains;
# fixed_k_power_lowbias; start_check), au_results_T31.json; regge/r2/fss/fss_results_T30.json (autocorrelation on 4^4 to 10^4, as at T30).
_A31 = json.load(open(_ND + 'a1_results_T31.json')); _U31 = json.load(open(_ND + 'au_results_T31.json'))
_P31 = _A31['fixed_k_power']; _P31b = _A31['fixed_k_power_lowbias']
for _k in ('0.2440', '0.2480', '0.2520'):
    assert _P31[_k]['L'] == [4, 5, 6, 8, 10, 12, 16, 20], (_k, _P31[_k]['L'])
    assert _A31['fixed_k'][_k]['16']['lineages'] == 2 and _A31['fixed_k'][_k]['20']['lineages'] == 1
_c31 = [r for r in _A31['chains'] if r['L'] >= 16 and r['finished']]
assert len(_c31) == 12 and all(r['chi'] is not None for r in _c31)
_fk31 = lambda k: [v['chi'] for v in _A31['fixed_k'][k].values()]
# autocorrelation near the crossing: medians per grid (4^4 to 12^4 as at T30; 16^4 and 20^4 from all their finished chains)
_tau31 = {L: v for L, v in _tauT.items() if L <= 12}
for _L in (16, 20): _tau31[_L] = [r['tau'] for r in _c31 if r['L'] == _L and 0.2425 <= r['k'] <= 0.2575]
_tm31 = [float(np.median(v)) for v in _tau31.values()]
# zero crossing of the average curvature on 16^4 and 20^4 from the chain averages per coupling (20^4: a short extrapolation above 0.252)
def _cross31(L):
    by = {}
    for r in _c31:
        if r['L'] == L: by.setdefault(r['k'], []).append(r['R'])
    p = sorted((k, float(np.mean(v))) for k, v in by.items())
    for (k1, r1), (k2, r2) in zip(p, p[1:]):
        if r1 < 0 <= r2: return k1 + (k2 - k1) * (-r1) / (r2 - r1)
    (k1, r1), (k2, r2) = p[-2], p[-1]; c = k1 + (k2 - k1) * (-r1) / (r2 - r1)
    assert r2 < 0 and c - k2 < 0.0005, (L, p, c); return c
_cr31 = [v['crossing'] for v in _MT['per_L'].values()] + [_cross31(16), _cross31(20)]
assert '%.4f' % min(_cr31) == T30_N['cr_lo'] and '%.4f' % max(_cr31) == T30_N['cr_hi']   # the 16^4 and 20^4 crossings fall inside the earlier range
# start check on 16^4: chains from another coupling's final configuration against flat starts
_s31 = _A31['start_check']; _z31 = [_s31['A1_L16_k0.248 - A2_L16_k0.248_r2']['z'], _s31['A1_L16_k0.252 - A2_L16_k0.252_r2']['z']]
assert max(abs(z) for z in _z31) < 1.2, _z31
# Hamber's (1993) fitted curve at k = 0.244: its rise from L = 4 to 20 as a power, and the gap to the fitted power here
_h31 = lambda L: (0.16 + 1.11 / L) * L ** 0.64
_hp31 = math.log(_h31(20) / _h31(4)) / math.log(5)
_g31a = (_hp31 - _P31['0.2440']['power']) / _P31['0.2440']['se']; _g31b = (_hp31 - _P31b['0.2440']['power']) / _P31b['0.2440']['se']
assert 0.35 < _h31(20) / _h31(4) - 1 < 0.42 and 2.0 < _g31b < _g31a < 3.5, (_h31(20) / _h31(4), _g31a, _g31b)
T31_N = dict(crit='%.1f' % 5 ** 0.94,
             f244_lo='%.1f' % min(_fk31('0.2440')), f244_hi='%.1f' % max(_fk31('0.2440')), p244=_p9(_P31['0.2440']),
             f248_lo='%.1f' % min(_fk31('0.2480')), f248_hi='%.1f' % max(_fk31('0.2480')), p248=_p9(_P31['0.2480']),
             f252_lo='%.1f' % min(_fk31('0.2520')), f252_hi='%.1f' % max(_fk31('0.2520')), p252=_p9(_P31['0.2520']),
             chi_lo=T30_N['chi_lo'], chi_hi=T30_N['chi_hi'], pw=T30_N['pw'], pw2=T30_N['pw2'],
             tau_lo='%d' % (10 * math.floor(min(_tm31) / 10)), tau_hi='%d' % (10 * math.ceil(max(_tm31) / 10)),
             cr_lo='%.4f' % min(_cr31), cr_hi='%.4f' % max(_cr31), hp='%.2f' % _hp31, g1='%.1f' % _g31a, g2='%.1f' % _g31b,
             c14a=T29_N['c14a'], c14b=T29_N['c14b'])
T31_OLD5 = md[md.index("Larger grids test whether the transition is critical."):]; T31_OLD5 = T31_OLD5[:T31_OLD5.index('\n\n')]
T31_P5 = ("Larger grids test whether the transition is critical. At a critical point the curvature fluctuation $χ = 4\\,\\mathrm{var}(Σ_h δ_h A_h)/\\langle Σ_s V_s \\rangle$ "
          "grows with the grid as $L^{2/ν - 4}$ (Hamber, 2015), as $L^{0.94}$ for the ν = 0.405 that δ = 0.62 implies, a factor of <<crit>> from 4⁴ to 20⁴, "
          "and the chains slow down as the grid grows. Chains on grids of 4⁴ to 20⁴ sites at k = 0.244 to 0.255 show neither. Combined across couplings "
          "by reweighting (Ferrenberg and Swendsen, 1989) on grids up to 12⁴, with two chains from independent starts per coupling on 16⁴ (and two forks "
          "at k = 0.244) and one chain on 20⁴, the fluctuation at k = 0.244, Hamber's $k_c$, is <<f244_lo>> to <<f244_hi>> on every grid, a fitted power "
          "of <<p244>> in L; at k = 0.248 it is <<f248_lo>> to <<f248_hi>>, a power of <<p248>>, and at k = 0.252 it is <<f252_lo>> to <<f252_hi>>, "
          "a power of <<p252>>. Its largest value near the crossing, in a range that includes Hamber's $k_c$, is <<chi_lo>> to <<chi_hi>> on grids up to 12⁴, "
          "a power of <<pw>>, or <<pw2>> from the second halves of the chains. Autocorrelation times measured within chains near the crossing stay between "
          "<<tau_lo>> and <<tau_hi>> sweeps, with no rise from 4⁴ to 20⁴. The average curvature crosses zero between k = <<cr_lo>> and <<cr_hi>> on grids "
          "of 4⁴ to 20⁴, a little above Hamber's fitted $k_c$, and its slope, which the fluctuation sets, does not steepen with size. Hamber (1993) fitted "
          "his fluctuation at k = 0.244 on grids up to 16⁴ as $L^{0.64(11)}$ with a large correction for small grids; that fitted curve rises by about two "
          "fifths from L = 4 to 20, a power of <<hp>> over that range, <<g1>> standard errors above the power found here at his $k_c$, or <<g2>> after "
          "correcting the 16⁴ and 20⁴ values for the low bias of a variance from few independent sweeps, so these grids lean against his fitted growth "
          "without ruling it out. Grids below 4⁴, which Hamber (1993) left out because the term couples next-nearest neighbors, are left out here too. "
          "On 12⁴, 16⁴ and 20⁴, where shorter chains still drift, each chain runs at least 1,500 sweeps at its coupling before 1,500 to 4,000 measured "
          "ones. Even so, on 12⁴ the chains continued from earlier runs sit 0.001 to 0.0025 below fresh flat starts in the average curvature at all six "
          "couplings, one to two standard errors each, a small memory of the start; their fluctuations show no such pattern. On 16⁴, chains started from "
          "the final configuration of the coupling above or below agree with flat starts within about one standard error. The two 14⁴ chains, of 1,600 "
          "measured sweeps each, are too short for the fits, and their fluctuations, <<c14a>> at k = 0.2475 and <<c14b>> at 0.2525, show no growth either.")
for _k, _v in T31_N.items(): T31_P5 = T31_P5.replace('<<' + _k + '>>', _v)
assert '<<' not in T31_P5
rep(T31_OLD5, T31_P5, label='Section 8.10, larger grids (T31: 16⁴ at Hamber\'s k_c and start check; 20⁴)')
# a = 0.02 (plans A2 and A4 add 10^4 and 12^4): the crossing stops moving; T29's "drifting down with size" rested on 4^4 to 8^4 and is corrected
_au31 = {(r['L'], round(r['k'], 3)): r for r in _U31['chains']}
_x8b = _au31[(8, 0.46)]; _cr8b = 0.46 - _x8b['R'] / _x8b['slope']; _cu = {int(L): v for L, v in _U31['crossing'].items()}
assert _cu[8] is None and abs(_cr8b - _cr8) < 1e-12 and _cu[4] > 0.475 and all(0.466 < _cu[L] < 0.469 for L in (6, 10, 12))
_g31 = max(_au31[(L, k)]['chi'] / _au31[(4, k)]['chi'] for (L, k) in _au31 if L > 4 and (4, k) in _au31)
_nu02 = (1 + 0.628) / 4; _tr = 3 ** (2 / _nu02 - 4)
assert 1.0 < _g31 < 1.2 and 2.5 < _tr < 3.0, (_g31, _tr)
T31_AU_NEW = ("the fluctuation is a tenth as large and grows by at most a fifth on grids from 4⁴ to 12⁴ (at k = 0.46 it is %.2f, %.2f, %.2f, %.2f and %.2f "
              "on 4⁴, 6⁴, 8⁴, 10⁴ and 12⁴), where critical growth with his δ = 0.628 would nearly triple it, and the curvature crosses zero at k = %.3f on 4⁴, "
              "%.3f on 6⁴, about %.3f on 8⁴ (from one coupling) and %.3f and %.3f on 10⁴ and 12⁴, with no drift toward his $k_c$ beyond 6⁴."
              % tuple([_au31[(L, 0.46)]['chi'] for L in (4, 6, 8, 10, 12)] + [_cu[4], _cu[6], _cr8b, _cu[10], _cu[12]]))
rep(T29_AU_NEW, T31_AU_NEW, label='Section 8.10, a = 0.02 (T31: 10⁴ and 12⁴; the crossing stops moving, correcting "drifting down with size")')
rep("per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T30.",
    "per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T31.", label='Section 8.10, supplement (T31)')
rep("grids up to 16⁴ show no critical growth yet, which leaves the critical exponent,",
    "grids up to 20⁴ show no critical growth yet, which leaves the critical exponent,", label='Section 8.10, costs (T31)')
rep("no critical growth up to 16⁴ here (Section 8.10) |", "no critical growth up to 20⁴ here (Section 8.10) |", label='Table 12, gravity row (T31)')
rep("On grids up to 16⁴ that transition shows no critical growth,", "On grids up to 20⁴ that transition shows no critical growth,", label='Conclusion (T31)')
