
# ------------------------------------------------------------------ SA74 T32 (Nov 10; Gary 19:03 forwarded the portfolio session's review of the T31 package)
# Editorial, no result changed: (1) the one authorial "we" goes; (2) the abstract states the larger-grid result and stays at 250 words;
# (3) every paragraph over 250 words, and Section 8.10's 248-word one, is broken at its topic changes; in Section 8.10's larger-grids
# paragraph the reweighting clause becomes a sentence of its own so that the new paragraph opens on its numbers; (4) the autocorrelation
# sentence says "median", which is what the 120 to 240 sweeps are (single chains on 12^4 run from 76 to 354); (5) the a = 0.02 values carry
# their standard errors (au_results_T31.json; crossings: the two bracketing curvature errors carried through the interpolation).
rep("None of these, as far as we can find, reports a short rule", "None of these appears to report a short rule", label='Introduction: the authorial "we" removed (T32)')
rep("transitive percolation reaches four on the chain at a million elements while its midpoint drifts to 5.8. From the top,",
    "transitive percolation's chain reaches four at a million elements while its midpoint drifts to 5.8. From the top,", label='Abstract: two words saved (T32)')
rep("All four groups in a new 24-chain Regge audit fail the mixing screen; a curvature-squared term, which quantum matter requires, makes all four pass.",
    "All four groups in a Regge audit fail the mixing screen; a curvature-squared term, which quantum matter requires, makes all pass, without critical growth through 20⁴.",
    label='Abstract: the larger-grid result (T32)')
_abs32 = md[md.index('## Abstract'):md.index('Keywords:')].replace('## Abstract', '').strip()
assert len(_abs32.split()) == 250, len(_abs32.split())
# Section 8.10, larger grids: the paragraph on the numbers opens with them; the reweighting clause follows as its own sentence
_rw = ("by reweighting (Ferrenberg and Swendsen, 1989) on grids up to 12⁴, with two chains from independent starts per coupling on 16⁴ (and two forks at k = 0.244) and one chain on 20⁴")
rep("Combined across couplings " + _rw + ", the fluctuation at k = 0.244, Hamber's $k_c$, is ", "At k = 0.244, Hamber's $k_c$, the fluctuation is ", label='Section 8.10, larger grids: the numbers first (T32)')
_m = re.search(r"a power of -?\d\.\d\d\(\d+\)\. Its largest value near the crossing,", md); assert _m and md.count(_m.group(0)) == 1
rep(_m.group(0), _m.group(0).replace(" Its largest value near the crossing,", " These values combine the chains across couplings " + _rw + ". The fluctuation's largest value near the crossing,"),
    label='Section 8.10, larger grids: the reweighting sentence (T32)')
rep("Autocorrelation times measured within chains near the crossing stay between ", "Median autocorrelation times of chains near the crossing stay between ", label='Section 8.10: "median" (T32)')
_mt = re.search(r"Median autocorrelation times of chains near the crossing stay between (\d+) and (\d+) sweeps, with no rise from 4⁴ to 20⁴\.", md); assert _mt
rep(_mt.group(0), _mt.group(0).replace(" sweeps, with no rise", " sweeps on every grid, with no rise"), label='Section 8.10: per grid (T32)')
assert all(int(_mt.group(1)) <= m <= int(_mt.group(2)) for m in _tm31)
# a = 0.02: standard errors
_e = lambda x, s, n=2: ('%.' + str(n) + 'f(%d)') % (x, max(1, round(s * 10 ** n)))
_au32 = {(r['L'], round(r['k'], 3)): r for r in _U31['chains']}
def _cross32(L):
    q = sorted((r['k'], r['R'], r['R_se']) for (l, k), r in _au32.items() if l == L)
    (k1, r1, s1), (k2, r2, s2) = next((a, b) for a, b in zip(q, q[1:]) if a[1] < 0 <= b[1])
    return k1 + (k2 - k1) * (-r1) / (r2 - r1), abs(k2 - k1) / (r2 - r1) ** 2 * math.hypot(r2 * s1, r1 * s2)
_c32 = {L: _cross32(L) for L in (4, 6, 10, 12)}
assert all(abs(_c32[L][0] - _cu[L]) < 1e-9 for L in _c32)
_old_chi = "(at k = 0.46 it is %.2f, %.2f, %.2f, %.2f and %.2f on 4⁴, 6⁴, 8⁴, 10⁴ and 12⁴)" % tuple(_au32[(L, 0.46)]['chi'] for L in (4, 6, 8, 10, 12))
_new_chi = "(at k = 0.46 it is %s, %s, %s, %s and %s on 4⁴, 6⁴, 8⁴, 10⁴ and 12⁴)" % tuple(_e(_au32[(L, 0.46)]['chi'], _au32[(L, 0.46)]['chi_se']) for L in (4, 6, 8, 10, 12))
rep(_old_chi, _new_chi, label='Section 8.10, a = 0.02: standard errors of the fluctuation (T32)')
_old_cr = "the curvature crosses zero at k = %.3f on 4⁴, %.3f on 6⁴, about %.3f on 8⁴ (from one coupling) and %.3f and %.3f on 10⁴ and 12⁴," % (_cu[4], _cu[6], _cr8b, _cu[10], _cu[12])
_new_cr = "the curvature crosses zero at k = %s on 4⁴, %s on 6⁴, about %.3f on 8⁴ (from one coupling) and %s and %s on 10⁴ and 12⁴," % (
    _e(*_c32[4], 3), _e(*_c32[6], 3), _cr8b, _e(*_c32[10], 3), _e(*_c32[12], 3))
rep(_old_cr, _new_cr, label='Section 8.10, a = 0.02: standard errors of the crossings (T32)')
# paragraph breaks at topic changes (no word changed)
T32_BREAKS = [
    ("yet doubles the source of gravity.", "Section 8.7 records a proposed specification"),
    ("The target is the rule both routes reach.", "The paper's results are the census and its failures"),
    ("needs about $10^8$.", "Dimension can depend on scale:"),
    ("could read four at large scales.", "The test that counts is the trend with size:"),
    ("while its midpoint moves away (Section 5.1).", "A rule that stays three-dimensional at large scales,"),
    ("Section 8 runs the first tests.", "The hard parts of this route"),
    ("finite sizes it tests.", "The first cap is in hand for part of the known laws."),
    ("times per halving.", "The geometry and the gauge fields are supplied, not grown,"),
    ("What matter does at the center must come from the rule.", "Discrete quantum geometry shows what such dynamics can do."),
    ("instead of reaching a singularity.", "A bounce can survive the quantum matter it holds."),
    ("The fields are free and massless.", "The coupling to curvature squared, which the quantum fields leave free, decides the outcome."),
    ("it does not test four-dimensional SU(3) dynamics.", "The weak force acts on left-handed fields only."),
    ("that Lee and Yang (1956) proposed.", "The anomalies must cancel,"),
    ("has not been performed here.", "Joining general relativity with quantum theory"),
    ("Jurkiewicz and Loll, 2008).", "On the supplied lattice itself,"),
    ("show neither.", "At k = 0.244, Hamber's $k_c$, the fluctuation is"),
    ("from the second halves of the chains.", "Median autocorrelation times of chains near the crossing"),
    ("does not steepen with size.", "Hamber (1993) fitted his fluctuation at k = 0.244"),
    ("are left out here too.", "On 12⁴, 16⁴ and 20⁴, where shorter chains still drift,"),
    ("as a local defect's share of an average should.", "Rough starts carry many such clusters."),
    ("with no sign of a second bulk phase.", "Toward smaller a, mixing worsens:"),
    ("would lead one to expect.", "At a = 0.02, where Hamber (1993) found"),
    ("with no drift toward his $k_c$ beyond 6⁴.", "The term buys an ensemble that mixes;"),
    ("differ in exponent.", "The term also bounds the action below."),
    ("the bound stops short of $k_c$.", "The term has costs."),
    ("the opposite of an expanding universe.", "A rule that rewrites a few relations at a time bounds the links of every update,"),
    ("runs its chains too long.", "The search goes on with longer rules,"),
    ("remain separate requirements.", "The independent Regge audit fails the prespecified mixing screen"),
    ("the continuum limit the gravity sector needs remains open.", "Harder tests sharpen the search:"),
]
for _a, _b in T32_BREAKS:
    rep(_a + ' ' + _b, _a + '\n\n' + _b, label='Paragraph break (T32): ' + _b[:40])
_long32 = [len(l.split()) for l in md.split('\n') if not l.startswith(('|', '#')) and len(l.split()) > 250]
assert not _long32, _long32
