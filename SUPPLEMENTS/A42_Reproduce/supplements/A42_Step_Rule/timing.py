"""Time the compact rule against the supplied numpy step and the spectral reference (massive Dirac, curved static 1+1).
usage: python3 timing.py [sizes...]   default sizes 256 1024 4096; the spectral reference runs only up to 1024 sites."""
import os, sys, time, importlib.util
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from A42_Step_Rule import step
SUP = os.path.join(H, '..', 'A42_Supplied_Physics', 'massive', 'massive_curved.py')
s = importlib.util.spec_from_file_location('mc', SUP); mc = importlib.util.module_from_spec(s); s.loader.exec_module(mc)
for n in [int(a) for a in sys.argv[1:]] or (256, 1024, 4096):
    t = time.time(); ref, rec = mc.local_circuit(n, .25, 2., .9); tn = time.time() - t
    dx, x, N, Nmid, y, Y = mc.geometry(n, .25); dt = rec['dt']; steps = rec['steps']
    psi = mc.initial(n, .25, .9); F = [list(psi[:, 0]), list(psi[:, 1])]; V = [list(Nmid), list(-Nmid)]; M = list(.5 * dt * 2. * N)
    t = time.time()
    for k in range(steps): step(F, V, None, M, dt, dx, k)
    tp = time.time() - t
    if n <= 1024:
        t = time.time(); sp, _ = mc.spectral(n, .25, 2., .9); ts = time.time() - t
    else: ts = float('nan')
    print(f'n={n:5d} steps={steps:5d} numpy local {tn:7.3f}s  pure-python rule {tp:7.2f}s  spectral reference {ts:7.3f}s  site-updates/s (numpy) {n*steps/tn:,.0f}')
