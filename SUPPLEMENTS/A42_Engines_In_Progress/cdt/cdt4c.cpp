// SA74 cdt4.cpp: causal dynamical triangulations in 3+1 dimensions, written from scratch (November 9, 2026), after cdt3.cpp.
// Spatial slices are 3-spheres, time is periodic with T slices, the ensemble is strict simplicial manifolds.
// Action (Ambjorn, Jurkiewicz, Loll): S = -(k0 + 6 D) N0 + k4 N4 + D (2 N41 + N32) + eps (N41 - V)^2, where N41 counts (4,1) and (1,4)
// simplices, N32 counts (3,2) and (2,3), N4 = N41 + N32.
// Simplex types in the slab [tt, tt+1]: 0 = (4,1) (v[0..3] a spatial tetrahedron at tt, v[4] the apex at tt+1), 1 = (1,4) (v[0..3] at tt+1,
// apex v[4] at tt), 2 = (3,2), 3 = (2,3). n[i] is the neighbor across the face opposite v[i].
// Moves and their selection (detailed balance; dS is the change of S):
//   (2,8): uniform (4,1) simplex, a vertex goes into its spatial tetrahedron; accept min(1, Nup/(N0+1) e^{-dS}), Nup = number of (4,1).
//   (8,2): uniform vertex, eligible if its star is four (4,1) and four (1,4) and its slice has at least 6 tetrahedra; min(1, N0/(Nup-3) e^{-dS}).
//   (4,6): uniform (4,1) and one of its 4 spatial triangles (two spatial tetrahedra become three around a new spatial edge);
//          min(1, Nup/(Nup+1) e^{-dS}).   (6,4): uniform (4,1) and one of its 6 spatial edges, which must have exactly three spatial
//          tetrahedra; min(1, Nup/(Nup-1) e^{-dS}).
//   (2,4): uniform simplex and one of its 5 faces, both simplices in one slab, the two opposite vertices at different times (a new
//          timelike edge); min(1, N4/(N4+2) e^{-dS}).   (4,2): uniform simplex and one of its 10 edges, timelike with exactly four
//          simplices around it; min(1, N4/(N4-2) e^{-dS}).
//   (3,3): uniform simplex and one of its 10 triangles, timelike with exactly three simplices around it, the new triangle timelike; accept 1.
// Usage: cdt4 T V k0 Delta k4 eps therm meas seed out.dat [check_every]   (V = target N41; sweep = V attempted moves; k4 tuned in thermalization)
// cdt4c (Nov 10, plan C): cdt4 plus the option KAVG (see main); runs without KAVG are unchanged
// Exact test: env CANON=Nmax (N4 capped), MOVES=bitmask of the 7 move kinds, THIN=moves between samples.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <string>
#include <unordered_map>
#include <algorithm>
using namespace std;
struct Sx { int v[5], n[5]; int type, tt; };
static vector<Sx> sx; static vector<int> freeS; static vector<char> liveS;
static vector<int> vtime, vsx, freeV; static vector<char> liveV;
static vector<int> bagU, posU, bagA, posA, bagV, posV;
static vector<int> sliceN; static int T; static long N0 = 0, N4 = 0;
static uint64_t s_[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnd() { uint64_t r = rotl(s_[1] * 5, 7) * 9, t = s_[1] << 17; s_[2] ^= s_[0]; s_[3] ^= s_[1]; s_[1] ^= s_[2]; s_[0] ^= s_[3]; s_[2] ^= t; s_[3] = rotl(s_[3], 45); return r; }
static inline double u01() { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static bool scripted = false; static vector<int> scr, radix; static size_t sp = 0; static double lastA = 1;
static inline int rint_(int n) { if (scripted) { if (sp >= scr.size()) scr.push_back(0); if (sp >= radix.size()) radix.push_back(n); else radix[sp] = n; return scr[sp++]; } return (int)(u01() * n); }
static void seedrng(uint64_t x) { for (int i = 0; i < 4; i++) { x += 0x9e3779b97f4a7c15ULL; uint64_t z = x; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; s_[i] = z ^ (z >> 31); } }

static int newV(int t) { int v; if (!freeV.empty()) { v = freeV.back(); freeV.pop_back(); } else { v = vtime.size(); vtime.push_back(0); vsx.push_back(-1); liveV.push_back(0); posV.push_back(-1); }
    vtime[v] = t; liveV[v] = 1; posV[v] = bagV.size(); bagV.push_back(v); N0++; return v; }
static void delV(int v) { liveV[v] = 0; int p = posV[v], w = bagV.back(); bagV[p] = w; posV[w] = p; bagV.pop_back(); freeV.push_back(v); N0--; }
static long att[7]; static int curmove = -1;
static int newS(const int *w, int tt) {   // the type follows from the vertex times; pyramids get their apex last
    int lo = 0, a[5]; for (int i = 0; i < 5; i++) { a[i] = w[i]; if (vtime[w[i]] == tt) lo++; else if (vtime[w[i]] != (tt + 1) % T) { fprintf(stderr, "newS: vertex outside the slab\n"); exit(3); } }
    int type; if (lo == 4) { type = 0; for (int i = 0; i < 4; i++) if (vtime[a[i]] != tt) swap(a[i], a[4]); }
    else if (lo == 1) { type = 1; for (int i = 0; i < 4; i++) if (vtime[a[i]] == tt) swap(a[i], a[4]); }
    else if (lo == 3) type = 2; else if (lo == 2) type = 3; else { fprintf(stderr, "move %d ", curmove); fprintf(stderr, "newS: bad type lo=%d tt=%d times %d %d %d %d %d; attempts %ld %ld %ld %ld %ld %ld %ld\n", lo, tt, vtime[a[0]], vtime[a[1]], vtime[a[2]], vtime[a[3]], vtime[a[4]], att[0], att[1], att[2], att[3], att[4], att[5], att[6]); exit(3); }
    int t; if (!freeS.empty()) { t = freeS.back(); freeS.pop_back(); } else { t = sx.size(); sx.push_back(Sx()); liveS.push_back(0); posU.push_back(-1); posA.push_back(-1); }
    Sx &x = sx[t]; for (int i = 0; i < 5; i++) { x.v[i] = a[i]; x.n[i] = -1; vsx[a[i]] = t; } x.type = type; x.tt = tt; liveS[t] = 1; N4++;
    posA[t] = bagA.size(); bagA.push_back(t); if (type == 0) { posU[t] = bagU.size(); bagU.push_back(t); sliceN[tt]++; } return t; }
static void delS(int t) { Sx &x = sx[t]; liveS[t] = 0; N4--; { int p = posA[t], w = bagA.back(); bagA[p] = w; posA[w] = p; bagA.pop_back(); }
    if (x.type == 0) { int p = posU[t], w = bagU.back(); bagU[p] = w; posU[w] = p; bagU.pop_back(); sliceN[x.tt]--; } freeS.push_back(t); }
static inline bool has(const Sx &x, int v) { return x.v[0] == v || x.v[1] == v || x.v[2] == v || x.v[3] == v || x.v[4] == v; }
static inline int idx(const Sx &x, int v) { for (int i = 0; i < 5; i++) if (x.v[i] == v) return i; return -1; }
static inline int common(const Sx &x, const Sx &y) { int c = 0; for (int i = 0; i < 5; i++) c += has(y, x.v[i]); return c; }
static void link(int a, int b) { Sx &x = sx[a], &y = sx[b]; int i = -1, j = -1; for (int k = 0; k < 5; k++) { if (!has(y, x.v[k])) i = k; if (!has(x, y.v[k])) j = k; } x.n[i] = b; y.n[j] = a; }
static void relink(const vector<int> &nw, const vector<int> &outer) {
    for (size_t i = 0; i < nw.size(); i++) { for (size_t j = i + 1; j < nw.size(); j++) if (common(sx[nw[i]], sx[nw[j]]) == 4) link(nw[i], nw[j]);
        for (int o : outer) if (common(sx[nw[i]], sx[o]) == 4) link(nw[i], o); } }
static void outers(const vector<int> &old, vector<int> &out) { out.clear(); for (int t : old) for (int i = 0; i < 5; i++) { int u = sx[t].n[i];
        if (find(old.begin(), old.end(), u) == old.end() && find(out.begin(), out.end(), u) == out.end()) out.push_back(u); } }
static vector<int> stampS; static int stamp = 0; static vector<int> starbuf;
static int star(int v, int cap) {
    if ((int)stampS.size() < (int)sx.size()) stampS.resize(sx.size() * 2 + 16, 0);
    stamp++; starbuf.clear(); starbuf.push_back(vsx[v]); stampS[vsx[v]] = stamp;
    for (size_t h = 0; h < starbuf.size(); h++) { const Sx &x = sx[starbuf[h]]; for (int i = 0; i < 5; i++) if (x.v[i] != v) { int u = x.n[i]; if (stampS[u] != stamp) { stampS[u] = stamp; starbuf.push_back(u); if ((int)starbuf.size() > cap) return starbuf.size(); } } }
    return starbuf.size(); }
static bool edgeExists(int a, int b) { star(a, 1 << 30); for (int t : starbuf) if (has(sx[t], b)) return true; return false; }
static bool triExists(int a, int b, int c) { star(a, 1 << 30); for (int t : starbuf) if (has(sx[t], b) && has(sx[t], c)) return true; return false; }
static bool tetExists(int a, int b, int c, int d) { star(a, 1 << 30); for (int t : starbuf) if (has(sx[t], b) && has(sx[t], c) && has(sx[t], d)) return true; return false; }

static double k0, Dl, k4, eps; static long Vtar, NMAX = 0;
static inline long N41() { return 2 * (long)bagU.size(); }
static inline double fix(long d) { return eps * ((double)(N41() + d - Vtar) * (N41() + d - Vtar) - (double)(N41() - Vtar) * (N41() - Vtar)); }
static long acc[7];
static inline bool take(double a) { if (scripted) { lastA = a < 1 ? a : 1; return true; } return a >= 1 || u01() < a; }
static void build(const vector<int> &old, const vector<array<int, 5>> &nwv, int tt_up, int tt_alt, const vector<int> &alt) {   // delete old, create new, relink
    vector<int> out; outers(old, out); for (int t : old) delS(t); vector<int> nw; for (size_t i = 0; i < nwv.size(); i++) nw.push_back(newS(nwv[i].data(), alt[i] ? tt_alt : tt_up)); relink(nw, out); }

static bool move28() { att[0]++; if (NMAX && N4 + 6 > NMAX) return false; int U = bagU[rint_(bagU.size())], D = sx[U].n[4];
    if (!take((double)bagU.size() / (N0 + 1) * exp(k0 - 6 * k4 - 6 * Dl - fix(6)))) return false;
    Sx u = sx[U], d = sx[D]; int v = newV(vtime[u.v[0]]); vector<array<int, 5>> nw; vector<int> alt;
    for (int i = 0; i < 4; i++) { array<int, 5> a = {u.v[0], u.v[1], u.v[2], u.v[3], u.v[4]}; a[i] = v; nw.push_back(a); alt.push_back(0); a[4] = d.v[4]; nw.push_back(a); alt.push_back(1); }
    build({U, D}, nw, u.tt, d.tt, alt); acc[0]++; return true; }
static bool move82() { att[1]++; int v = bagV[rint_(bagV.size())]; if (sliceN[vtime[v]] < 6) return false; if (star(v, 8) != 8) return false;
    int n0 = 0, n1 = 0, U = -1, D = -1; for (int t : starbuf) { const Sx &x = sx[t]; if (x.type > 1 || idx(x, v) == 4) return false; if (x.type == 0) { n0++; U = t; } else { n1++; D = t; } }
    if (n0 != 4 || n1 != 4) return false;
    if (!take((double)N0 / (bagU.size() - 3) * exp(-k0 + 6 * k4 + 6 * Dl - fix(-6)))) return false;
    vector<int> old(starbuf.begin(), starbuf.end()); set<int> nb; for (int t : old) for (int i = 0; i < 4; i++) if (sx[t].v[i] != v) nb.insert(sx[t].v[i]);
    if (nb.size() != 4) { fprintf(stderr, "move82: link is not a tetrahedron\n"); exit(3); }
    vector<int> b(nb.begin(), nb.end()); int du = sx[U].v[4], dd = sx[D].v[4], tu = sx[U].tt, td = sx[D].tt;
    build(old, {{b[0], b[1], b[2], b[3], du}, {b[0], b[1], b[2], b[3], dd}}, tu, td, {0, 1}); delV(v); acc[1]++; return true; }
static bool move46() { att[2]++; if (NMAX && N4 + 2 > NMAX) return false; int U1 = bagU[rint_(bagU.size())], i = rint_(4); Sx u1 = sx[U1]; int U2 = u1.n[i]; if (sx[U2].type != 0) return false;
    int D1 = u1.n[4], D2 = sx[U2].n[4], x = u1.v[i]; if (sx[D1].n[idx(sx[D1], x)] != D2) return false;
    int f[3], m = 0; for (int k = 0; k < 4; k++) if (k != i) f[m++] = u1.v[k]; int y = -1; for (int k = 0; k < 4; k++) if (sx[U2].v[k] != f[0] && sx[U2].v[k] != f[1] && sx[U2].v[k] != f[2]) y = sx[U2].v[k];
    if (edgeExists(x, y)) return false; if (!take((double)bagU.size() / (bagU.size() + 1) * exp(-2 * k4 - 4 * Dl - fix(2)))) return false;
    int e = u1.v[4], g = sx[D1].v[4], tu = u1.tt, td = sx[D1].tt; vector<array<int, 5>> nw; vector<int> alt;
    for (int k = 0; k < 3; k++) { nw.push_back({x, y, f[k], f[(k + 1) % 3], e}); alt.push_back(0); nw.push_back({x, y, f[k], f[(k + 1) % 3], g}); alt.push_back(1); }
    build({U1, U2, D1, D2}, nw, tu, td, alt); acc[2]++; return true; }
static bool move64() { att[3]++; int U = bagU[rint_(bagU.size())]; static const int pr[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}; int r = rint_(6); Sx u = sx[U];
    int x = u.v[pr[r][0]], y = u.v[pr[r][1]], ab[2], m = 0; for (int k = 0; k < 4; k++) if (k != pr[r][0] && k != pr[r][1]) ab[m++] = u.v[k]; int a = ab[0], b = ab[1];
    int U2 = u.n[idx(u, a)]; if (sx[U2].type != 0) return false; int c = -1; for (int k = 0; k < 4; k++) { int w = sx[U2].v[k]; if (w != x && w != y && w != b) c = w; }
    int U3 = sx[U2].n[idx(sx[U2], b)]; if (sx[U3].type != 0 || !has(sx[U3], a) || U3 == U) return false; if (idx(sx[U3], a) == 4) return false;
    int D1 = u.n[4], D2 = sx[U2].n[4], D3 = sx[U3].n[4]; if (sx[D1].n[idx(sx[D1], a)] != D2 || sx[D2].n[idx(sx[D2], b)] != D3) return false;
    if (triExists(a, b, c)) return false; if (!take((double)bagU.size() / (bagU.size() - 1) * exp(2 * k4 + 4 * Dl - fix(-2)))) return false;
    int e = u.v[4], g = sx[D1].v[4], tu = u.tt, td = sx[D1].tt;
    build({U, U2, U3, D1, D2, D3}, {{a, b, c, x, e}, {a, b, c, y, e}, {a, b, c, x, g}, {a, b, c, y, g}}, tu, td, {0, 0, 1, 1}); acc[3]++; return true; }
static bool move24() { att[4]++; if (NMAX && N4 + 2 > NMAX) return false; int S1 = bagA[rint_(bagA.size())], i = rint_(5); Sx s = sx[S1]; int S2 = s.n[i]; if (sx[S2].tt != s.tt) return false;
    int p = s.v[i], q = -1; for (int k = 0; k < 5; k++) if (!has(s, sx[S2].v[k])) q = sx[S2].v[k]; if (vtime[p] == vtime[q]) return false; if (edgeExists(p, q)) return false;
    if (!take((double)N4 / (N4 + 2) * exp(-2 * k4 - 2 * Dl))) return false;
    int f[4], m = 0; for (int k = 0; k < 5; k++) if (k != i) f[m++] = s.v[k]; vector<array<int, 5>> nw; for (int k = 0; k < 4; k++) nw.push_back({p, q, f[k], f[(k + 1) % 4], f[(k + 2) % 4]});
    build({S1, S2}, nw, s.tt, s.tt, {0, 0, 0, 0}); acc[4]++; return true; }
static bool move42() { att[5]++; int S1 = bagA[rint_(bagA.size())]; Sx s = sx[S1]; int i = rint_(5), j = rint_(4); if (j >= i) j++; int p = s.v[i], q = s.v[j]; if (vtime[p] == vtime[q]) return false;
    vector<int> ring = {S1}; for (size_t h = 0; h < ring.size(); h++) { const Sx &x = sx[ring[h]]; for (int k = 0; k < 5; k++) if (x.v[k] != p && x.v[k] != q) { int u = x.n[k]; if (find(ring.begin(), ring.end(), u) == ring.end()) { ring.push_back(u); if (ring.size() > 4) return false; } } }
    if (ring.size() != 4) return false; set<int> lk; for (int t : ring) for (int k = 0; k < 5; k++) if (sx[t].v[k] != p && sx[t].v[k] != q) lk.insert(sx[t].v[k]); if (lk.size() != 4) return false;
    vector<int> f(lk.begin(), lk.end()); if (vtime[f[0]] == vtime[f[1]] && vtime[f[1]] == vtime[f[2]] && vtime[f[2]] == vtime[f[3]]) return false;   // the new shared tetrahedron must be timelike
    if (tetExists(f[0], f[1], f[2], f[3])) return false; if (!take((double)N4 / (N4 - 2) * exp(2 * k4 + 2 * Dl))) return false;
    build(ring, {{f[0], f[1], f[2], f[3], p}, {f[0], f[1], f[2], f[3], q}}, s.tt, s.tt, {0, 0}); acc[5]++; return true; }
static bool move33() { att[6]++; int S1 = bagA[rint_(bagA.size())]; Sx s = sx[S1]; int i = rint_(5), j = rint_(4); if (j >= i) j++; int x = s.v[i], y = s.v[j], t3[3], m = 0; for (int k = 0; k < 5; k++) if (k != i && k != j) t3[m++] = s.v[k];
    if (vtime[t3[0]] == vtime[t3[1]] && vtime[t3[1]] == vtime[t3[2]]) return false;
    int S2 = s.n[i]; int z = -1; for (int k = 0; k < 5; k++) if (!has(s, sx[S2].v[k])) z = sx[S2].v[k]; int S3 = sx[S2].n[idx(sx[S2], y)]; if (S3 == S1 || !has(sx[S3], x)) return false;
    if (vtime[x] == vtime[y] && vtime[y] == vtime[z]) return false; if (triExists(x, y, z)) return false;
    build({S1, S2, S3}, {{x, y, z, t3[0], t3[1]}, {x, y, z, t3[1], t3[2]}, {x, y, z, t3[0], t3[2]}}, s.tt, s.tt, {0, 0, 0}); acc[6]++; return true; }

static void domove(int r) { curmove = r; switch (r) { case 0: move28(); break; case 1: move82(); break; case 2: move46(); break; case 3: move64(); break; case 4: move24(); break; case 5: move42(); break; default: move33(); } }

static void check(const char *where) {
#define FAIL(msg) { fprintf(stderr, "CHECK FAILED (%s): %s\n", where, msg); exit(2); }
    map<array<int, 4>, int> faces; long n4 = 0, nu = 0; vector<long> sl(T, 0), inc(vtime.size(), 0); vector<map<array<int, 3>, int>> stri(T);
    for (size_t t = 0; t < sx.size(); t++) if (liveS[t]) { const Sx &x = sx[t]; n4++; int lo = 0, hi = 0;
        for (int i = 0; i < 5; i++) { if (!liveV[x.v[i]]) FAIL("dead vertex"); inc[x.v[i]]++; for (int j = i + 1; j < 5; j++) if (x.v[i] == x.v[j]) FAIL("repeated vertex");
            if (vtime[x.v[i]] == x.tt) lo++; else if (vtime[x.v[i]] == (x.tt + 1) % T) hi++; else FAIL("vertex time");
            int u = x.n[i]; if (u < 0 || !liveS[u]) FAIL("dead neighbor"); if (has(sx[u], x.v[i])) FAIL("neighbor has opposite vertex"); if (common(x, sx[u]) != 4) FAIL("neighbor shares != 4");
            int back = 0; for (int k = 0; k < 5; k++) if (sx[u].n[k] == (int)t) back++; if (back != 1) FAIL("back pointer");
            array<int, 4> f; int m = 0; for (int k = 0; k < 5; k++) if (k != i) f[m++] = x.v[k]; sort(f.begin(), f.end()); faces[f]++; }
        int ty = lo == 4 ? 0 : lo == 1 ? 1 : lo == 3 ? 2 : lo == 2 ? 3 : -1; if (ty != x.type) FAIL("type");
        if (x.type == 0) { if (vtime[x.v[4]] != (x.tt + 1) % T) FAIL("apex up"); nu++; sl[x.tt]++; if (sx[x.n[4]].type != 1) FAIL("base partner");
            for (int i = 0; i < 4; i++) { array<int, 3> f; int m = 0; for (int k = 0; k < 4; k++) if (k != i) f[m++] = x.v[k]; sort(f.begin(), f.end()); stri[x.tt][f]++; } }
        if (x.type == 1 && vtime[x.v[4]] != x.tt) FAIL("apex down"); }
    for (auto &f : faces) if (f.second != 2) FAIL("tetrahedron not shared by exactly two simplices");
    for (int t = 0; t < T; t++) { for (auto &f : stri[t]) if (f.second != 2) FAIL("spatial triangle not shared by exactly two spatial tetrahedra"); if (sl[t] != sliceN[t] || sl[t] < 5) FAIL("slice count"); }
    long n0 = 0; for (size_t v = 0; v < vtime.size(); v++) if (liveV[v]) { n0++; if (!liveS[vsx[v]] || !has(sx[vsx[v]], v)) FAIL("vsx"); if (star(v, 1 << 30) != inc[v]) FAIL("vertex star not connected"); }
    if (n4 != N4 || nu != (long)bagU.size() || n0 != N0 || n0 != (long)bagV.size() || n4 != (long)bagA.size()) FAIL("counts");
    long chi = 0; { set<array<int, 2>> E; set<array<int, 3>> F; for (size_t t = 0; t < sx.size(); t++) if (liveS[t]) { const Sx &x = sx[t]; for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) { E.insert({min(x.v[i], x.v[j]), max(x.v[i], x.v[j])});
                    for (int k = j + 1; k < 5; k++) { array<int, 3> f = {x.v[i], x.v[j], x.v[k]}; sort(f.begin(), f.end()); F.insert(f); } } }
        chi = n0 - (long)E.size() + (long)F.size() - (long)faces.size() + n4; }
    if (chi != 0) FAIL("Euler characteristic of S3 x S1 is not 0");
}
static void canonical(string &best, int &aut) {
    vector<int> live; for (size_t t = 0; t < sx.size(); t++) if (liveS[t]) live.push_back(t);
    int pm[5] = {0, 1, 2, 3, 4}; vector<array<int, 5>> perms; do { perms.push_back({pm[0], pm[1], pm[2], pm[3], pm[4]}); } while (next_permutation(pm, pm + 5));
    vector<int> lab(vtime.size()), seen(sx.size()), q; string code; best.clear(); aut = 0; bool first = true;
    for (int t0 : live) for (auto &p : perms) { fill(lab.begin(), lab.end(), -1); fill(seen.begin(), seen.end(), 0); q.clear(); int nl = 0; code.clear(); bool worse = false;
        for (int k = 0; k < 5; k++) { lab[sx[t0].v[p[k]]] = nl++; code.push_back((char)vtime[sx[t0].v[p[k]]]); } q.push_back(t0); seen[t0] = 1;
        for (size_t h = 0; h < q.size() && !worse; h++) { const Sx &x = sx[q[h]]; int o[5] = {0, 1, 2, 3, 4}; sort(o, o + 5, [&](int a, int b) { return lab[x.v[a]] < lab[x.v[b]]; });
            for (int k = 0; k < 5; k++) code.push_back((char)lab[x.v[o[k]]]);
            for (int k = 0; k < 5; k++) { int u = x.n[o[k]]; if (!seen[u]) { seen[u] = 1; q.push_back(u); for (int m = 0; m < 5; m++) if (lab[sx[u].v[m]] < 0) { lab[sx[u].v[m]] = nl++; code.push_back((char)(100 + vtime[sx[u].v[m]])); } } }
            if (!first && code.compare(0, code.size(), best, 0, code.size()) > 0) worse = true; }
        if (worse) continue; if (first || code < best) { best = code; aut = 1; first = false; } else if (code == best) aut++; } }

static void save_config(const char *fn) {   // vertices renumbered 0..N0-1; one line per simplex: five vertices and the slab time
    vector<int> id(vtime.size(), -1); int n = 0; for (size_t v = 0; v < vtime.size(); v++) if (liveV[v]) id[v] = n++;
    string tmp = string(fn) + ".tmp"; FILE *f = fopen(tmp.c_str(), "w"); fprintf(f, "cdt4config %d %ld %ld\n", T, N0, N4);
    for (size_t v = 0; v < vtime.size(); v++) if (liveV[v]) fprintf(f, "%d\n", vtime[v]);
    for (size_t t = 0; t < sx.size(); t++) if (liveS[t]) fprintf(f, "%d %d %d %d %d %d\n", id[sx[t].v[0]], id[sx[t].v[1]], id[sx[t].v[2]], id[sx[t].v[3]], id[sx[t].v[4]], sx[t].tt);
    fclose(f); rename(tmp.c_str(), fn); }
static void neighbors_from_faces(const vector<int> &all) {
    map<array<int, 4>, pair<int, int>> fm; for (int t : all) for (int i = 0; i < 5; i++) { array<int, 4> f; int m = 0; for (int k = 0; k < 5; k++) if (k != i) f[m++] = sx[t].v[k]; sort(f.begin(), f.end());
        auto it = fm.find(f); if (it == fm.end()) fm[f] = {t, i}; else { sx[t].n[i] = it->second.first; sx[it->second.first].n[it->second.second] = t; } } }
static bool load_config(const char *fn) { FILE *f = fopen(fn, "r"); if (!f) return false; int t_, tt; long n0, n4; if (fscanf(f, "cdt4config %d %ld %ld", &t_, &n0, &n4) != 3 || t_ != T) { fprintf(stderr, "load: bad header or T\n"); exit(3); }
    for (long v = 0; v < n0; v++) { int tm; if (fscanf(f, "%d", &tm) != 1) exit(3); newV(tm); } vector<int> all;
    for (long i = 0; i < n4; i++) { int w[5]; if (fscanf(f, "%d %d %d %d %d %d", &w[0], &w[1], &w[2], &w[3], &w[4], &tt) != 6) exit(3); all.push_back(newS(w, tt)); }
    fclose(f); neighbors_from_faces(all); return true; }
int main(int argc, char **argv) {
    if (argc < 11) { fprintf(stderr, "usage: cdt4 T V k0 Delta k4 eps therm meas seed out.dat [check_every]\n"); return 1; }
    T = atoi(argv[1]); Vtar = atol(argv[2]); k0 = atof(argv[3]); Dl = atof(argv[4]); k4 = atof(argv[5]); eps = atof(argv[6]); long therm = atol(argv[7]), meas = atol(argv[8]); seedrng(atoll(argv[9]));
    int chk = argc > 11 ? atoi(argv[11]) : 0; sliceN.assign(T, 0);
    // start: every slice the boundary of a 4-simplex (5 tetrahedra); each prism over a tetrahedron i<j<k<l is cut as
    // (bi bj bk bl | tl), (bi bj bk | tk tl), (bi bj | tj tk tl), (bi | ti tj tk tl)
    bool loaded = getenv("LOAD") && load_config(getenv("LOAD"));
    if (!loaded) {
    vector<array<int, 5>> sv(T); for (int t = 0; t < T; t++) for (int i = 0; i < 5; i++) sv[t][i] = newV(t);
    vector<int> all;
    for (int t = 0; t < T; t++) { auto &b = sv[t]; auto &u = sv[(t + 1) % T]; for (int o = 0; o < 5; o++) { int c[4], m = 0; for (int k = 0; k < 5; k++) if (k != o) c[m++] = k; int i = c[0], j = c[1], k = c[2], l = c[3];
            int w1[5] = {b[i], b[j], b[k], b[l], u[l]}, w2[5] = {b[i], b[j], b[k], u[k], u[l]}, w3[5] = {b[i], b[j], u[j], u[k], u[l]}, w4[5] = {b[i], u[i], u[j], u[k], u[l]};
            all.push_back(newS(w1, t)); all.push_back(newS(w2, t)); all.push_back(newS(w3, t)); all.push_back(newS(w4, t)); } }
    neighbors_from_faces(all); }
    check("start");
    if (getenv("DBTEST")) {   // exact local test of detailed balance: for sampled states x, every possible selection of every move is applied;
        // for every neighbor state y, every selection from y is applied; pi(x) P(x->y) must equal pi(y) P(y->x), pi = e^{-S}.
        struct St { vector<Sx> sx; vector<int> freeS; vector<char> liveS; vector<int> vtime, vsx, freeV; vector<char> liveV; vector<int> bagU, posU, bagA, posA, bagV, posV, sliceN; long N0, N4; };
        auto save = [&]() { return St{sx, freeS, liveS, vtime, vsx, freeV, liveV, bagU, posU, bagA, posA, bagV, posV, sliceN, N0, N4}; };
        auto load = [&](const St &z) { sx = z.sx; freeS = z.freeS; liveS = z.liveS; vtime = z.vtime; vsx = z.vsx; freeV = z.freeV; liveV = z.liveV; bagU = z.bagU; posU = z.posU; bagA = z.bagA; posA = z.posA; bagV = z.bagV; posV = z.posV; sliceN = z.sliceN; N0 = z.N0; N4 = z.N4; };
        auto hsh = [&]() { uint64_t h = 0; for (size_t t = 0; t < sx.size(); t++) if (liveS[t]) { int w[5]; for (int i = 0; i < 5; i++) w[i] = sx[t].v[i]; sort(w, w + 5); uint64_t z = 1469598103934665603ULL; for (int i = 0; i < 5; i++) { z ^= (uint64_t)w[i] + 0x9e3779b97f4a7c15ULL; z *= 1099511628211ULL; z ^= z >> 29; } h += z * z + z; } return h; };
        auto action = [&]() { double n41 = N41(); return -(k0 + 6 * Dl) * N0 + k4 * N4 + Dl * (2 * n41 + (N4 - n41)) + eps * (n41 - Vtar) * (n41 - Vtar); };
        // all transitions out of the current state: target hash -> (probability, a saved copy)
        auto expand = [&](map<uint64_t, double> &P, map<uint64_t, St> *keep, long *nsel) { St base = save(); scripted = true;
            for (int mv = 0; mv < 7; mv++) { scr.clear(); radix.clear();
                while (true) { sp = 0; lastA = 1; bool ok = false; switch (mv) { case 0: ok = move28(); break; case 1: ok = move82(); break; case 2: ok = move46(); break; case 3: ok = move64(); break; case 4: ok = move24(); break; case 5: ok = move42(); break; default: ok = move33(); }
                    scr.resize(sp); radix.resize(sp); if (nsel) (*nsel)++;
                    if (ok) { double g = 1.0 / 7; for (size_t i = 0; i < sp; i++) g /= radix[i]; uint64_t h = hsh(); P[h] += g * lastA; if (keep && !keep->count(h)) (*keep)[h] = save(); load(base); }
                    int i = (int)sp - 1; while (i >= 0 && scr[i] + 1 >= radix[i]) i--; if (i < 0) break; scr[i]++; scr.resize(i + 1); radix.resize(i + 1); } }
            scripted = false; load(base); };
        int nst = atoi(getenv("DBTEST")); double worst = 0; long pairs = 0, nsel = 0, oneway = 0;
        for (int it = 0; it < nst; it++) { for (long m = 0; m < therm; m++) domove(rint_(7)); check("dbtest state");
            St X = save(); uint64_t hx = hsh(); double Sx_ = action(); map<uint64_t, double> P; map<uint64_t, St> keep; expand(P, &keep, &nsel);
            for (auto &e : P) { load(keep[e.first]); double Sy = action(); map<uint64_t, double> Q; expand(Q, nullptr, &nsel); double back = Q.count(hx) ? Q[hx] : 0; if (back == 0) oneway++;
                double l = exp(-Sx_) * e.second, r = exp(-Sy) * back, d = fabs(l - r) / max(l, r); if (d > worst) worst = d; pairs++; }
            load(X); fprintf(stderr, "state %d: N0=%ld N4=%ld N41=%ld, %zu neighbor states, worst relative violation so far %.3g, one-way transitions %ld\n", it, N0, N4, N41(), P.size(), worst, oneway); }
        printf("DBTEST: %d states, %ld state pairs, %ld selections tried, worst relative violation of detailed balance %.3g, one-way transitions %ld\n", nst, pairs, nsel, worst, oneway); return worst < 1e-9 && oneway == 0 ? 0 : 4; }
    if (getenv("CANON")) { NMAX = atol(getenv("CANON")); eps = 0; int thin = getenv("THIN") ? atoi(getenv("THIN")) : 100, mask = getenv("MOVES") ? atoi(getenv("MOVES")) : 127;
        unordered_map<string, array<long, 4>> tally; string c; int aut;
        for (long i = 0; i < therm; i++) { for (int m = 0; m < thin; m++) { int r = rint_(7); if (mask >> r & 1) domove(r); } canonical(c, aut); auto &e = tally[c]; e[0]++; e[1] = N0; e[2] = N4; e[3] = aut; if (e[0] == 1) tally[c][3] = aut; }
        check("canon end"); FILE *fc = fopen(argv[10], "w"); fprintf(fc, "# count N0 N4 aut N41\n"); for (auto &e : tally) { int n41 = 0; (void)n41; fprintf(fc, "%ld %ld %ld %ld\n", e.second[0], e.second[1], e.second[2], e.second[3]); } fclose(fc);
        fprintf(stderr, "accepted: %ld %ld %ld %ld %ld %ld %ld\n", acc[0], acc[1], acc[2], acc[3], acc[4], acc[5], acc[6]); return 0; }
    FILE *fo = fopen(argv[10], "w"); fprintf(fo, "# cdt4 T=%d V(N41)=%ld k0=%g Delta=%g k4(start)=%g eps=%g therm=%ld meas=%ld seed=%s\n# sweep phase N0 N4 N41 k4 | N3(t), spatial tetrahedra per slice, t = 0..T-1\n", T, Vtar, k0, Dl, k4, eps, therm, meas, argv[9]);
    long outevery = getenv("OUTEVERY") ? atol(getenv("OUTEVERY")) : 1;
    long ramp = getenv("RAMP") ? atol(getenv("RAMP")) : 0, Vfinal = Vtar, V0 = getenv("RAMPFROM") ? atol(getenv("RAMPFROM")) : 2000;   // RAMP: the target volume grows linearly from RAMPFROM to V over the first RAMP sweeps
    const char *savefn = getenv("SAVE"); long saveevery = getenv("SAVEEVERY") ? atol(getenv("SAVEEVERY")) : 0;
    bool kavg = getenv("KAVG") != NULL; long kfrom = ramp > therm / 2 ? ramp : therm / 2, kn = 0; double ksum = 0;   // KAVG: at the end of thermalization k4 is set to its mean over the thermalization sweeps from max(RAMP, therm/2) on (without it k4 stays at its last tuned value)
    for (long sw = 0; sw < therm + meas; sw++) { double sumN = 0; long cnt = 0;
        if (kavg && sw == therm && kn > 0) k4 = ksum / kn;
        Vtar = (ramp && sw < ramp && V0 < Vfinal) ? V0 + (Vfinal - V0) * sw / ramp : Vfinal;
        for (long m = 0; m < Vfinal; m++) { domove(rint_(7)); if ((m & 63) == 0) { sumN += N41(); cnt++; } }
        if (sw < therm) { double d = sumN / cnt - Vtar; if (fabs(d) > 0.02 * Vtar + 5) k4 += (d > 0 ? 1 : -1) * 0.001; else k4 += 2 * eps * d * 0.2; if (kavg && sw >= kfrom) { ksum += k4; kn++; } }
        if (chk && sw % chk == 0) check("sweep");
        if (savefn && saveevery && sw % saveevery == 0 && sw > 0) { check("save"); if (getenv("SAVESERIES") && sw >= therm) { char b[600]; snprintf(b, 600, "%s.%ld", savefn, sw); save_config(b); } save_config(savefn); }
        if (sw % outevery) continue;
        fprintf(fo, "%ld %d %ld %ld %ld %.5f |", sw, sw >= therm, N0, N4, N41(), k4); for (int t = 0; t < T; t++) fprintf(fo, " %d", sliceN[t]); fprintf(fo, "\n"); if (sw % 20 == 0) fflush(fo); }
    check("end"); if (savefn) save_config(savefn);
    fprintf(fo, "# done; accepted (2,8) %ld (8,2) %ld (4,6) %ld (6,4) %ld (2,4) %ld (4,2) %ld (3,3) %ld of attempts %ld %ld %ld %ld %ld %ld %ld\n", acc[0], acc[1], acc[2], acc[3], acc[4], acc[5], acc[6], att[0], att[1], att[2], att[3], att[4], att[5], att[6]);
    fclose(fo); return 0; }
