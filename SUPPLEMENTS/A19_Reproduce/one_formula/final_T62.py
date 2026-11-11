"""T62: the final formula with its core (the model of PROTOCOL_T61, MOGF0IESBC, adopted on 10 November 2026 after all its
preset marks were met), with t and kappa at the values fitted on all 149 benchmark disks. It exposes the interface that
the earlier check scripts use for the formula without the core (harness_T48.P: .H, .B, .S, .P, .MODEL), so that the same
scripts run on the final formula with one import changed."""
import types
import harness_T61 as T61
T61.U['t'], T61.U['k'] = 0.1501, 0.6229
P = types.SimpleNamespace(**{k: getattr(T61.P, k) for k in dir(T61.P) if not k.startswith('__')})
P.MODEL = 'MOGF0IESBC'; P.B = T61.B; P.H = T61.H; P.S = T61.S; T = T61
