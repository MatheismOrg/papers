"""T41 (PROTOCOL_T41.txt): harness_T37 with the halo mass function chosen by the environment variable A19_HMF
(tinker08, despali16, watson13 or bocquet16; mdef 200c). Everything else as T37."""
import os
import harness_T37 as S
from colossus.lss import mass_function
MOD = os.environ.get('A19_HMF', 'tinker08')
def hmf(z):
    k = round(float(z), 6)
    if k not in S._HMF: S._HMF[k] = mass_function.massFunction(S.MG * S.HH, k, mdef='200c', model=MOD, q_out='dndlnM')
    return S._HMF[k]
S.hmf = hmf
H = S.H
