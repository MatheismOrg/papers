"""T39 (PROTOCOL_T39.txt): T37's MOGF0IES under the alias MOGF0IESY, for chunk_run_loco.py (CLASH out of every training set)."""
import harness_T37 as S
H = S.H
H.MODELS['MOGF0IESY'] = H.MODELS['MOGF0IES']
_g, _e = H.gpred, H.esd
def gpred(model, th, a): return _g('MOGF0IES' if model == 'MOGF0IESY' else model, th, a)
def esd(model, th, st): return _e('MOGF0IES' if model == 'MOGF0IESY' else model, th, st)
H.gpred, H.esd = gpred, esd
