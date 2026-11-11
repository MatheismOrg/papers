"""Build report, paper section and figure directly from the audit results."""
from pathlib import Path
import json,hashlib,platform,lzma,ast
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parent;root=D.parent
r=json.loads((D/'results.json').read_text());groups=r['groups'];runs=r['runs']
# Reproduce the previously reported description length, without running unrelated expensive tests.
base=root/'T18_ORIGINAL/SUPPLEMENTS'
a=ast.parse((base/'A42_SA74/A42_SA74_Supplied_Lattice.py').read_text());texts={}
for n in a.body:
 if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ['LATTICE','DECLARED']:
  texts[n.targets[0].id]=ast.literal_eval(n.value)
b=(base/'A38_Symbolic_Benchmark_Spec.txt').read_bytes()+texts['LATTICE']+texts['DECLARED']
bits=8*len(lzma.compress(b,format=lzma.FORMAT_RAW,filters=[{'id':lzma.FILTER_LZMA2,'preset':9|lzma.PRESET_EXTREME}]))
(D/'compression_check.json').write_text(json.dumps(dict(description_bits=bits,raw_bytes=len(b),lzma_version='Python stdlib; raw LZMA2 preset 9 extreme',is_complete_simulator=False),indent=2))
(D/'requirements.txt').write_text(f'numpy=={np.__version__}\nmatplotlib=={matplotlib.__version__}\n')
fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True)
for g,ax in zip(groups,axes.flat):
 for row in [x for x in runs if x['L']==g['L'] and x['k']==g['k']]:
  dat=np.loadtxt(D/'runs'/f'{row["tag"]}.dat')
  # Nonoverlapping 100-sweep means only for display; diagnostics use every sweep.
  n=len(dat)//100;t=dat[:n*100,0].reshape(n,100).mean(1);v=dat[:n*100,4].reshape(n,100).mean(1)
  ax.plot(t,v,color='#1f77b4' if row['start']=='flat' else '#d95f02',alpha=.75,lw=1,label=row['start'] if row['seed']%100 in (0,1) else None)
 ax.axvspan(0,3000,color='gray',alpha=.12);ax.axhline(0,color='black',lw=.5)
 ax.set_yscale('symlog',linthresh=2);ax.set_title(f'L={g["L"]}, k={g["k"]}; split R-hat={g["split_rhat"]:.3f}')
 ax.set_ylabel('Curvature R (symlog scale)');ax.grid(alpha=.2)
for ax in axes[-1]:ax.set_xlabel('Monte Carlo sweep')
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([0],[0],color='#1f77b4',label='Flat starts (3 seeds)'),Line2D([0],[0],color='#d95f02',label='Rough starts (3 seeds)')],loc='upper center',ncol=2)
# The PDF caption supplies the title; keep the legend unobstructed.
fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(D/'regge_audit.png',dpi=160,bbox_inches='tight');plt.close(fig)
lines=['### 8.9 Independent execution and equilibration audit','',
'A new audit reran the supplied Regge C engine without changing its action, measure or update rules. The protocol was fixed before the main runs. At each of L = 3 and 4 and k = 0.05 and 0.065, three seeded chains began flat and three began from independently prepared rough configurations. Every rough preparation used the same prescribed 3,000-sweep evolution at k = 0.1; no configuration was selected for its observed outcome. Each of the 24 main chains used 3,000 equilibration sweeps followed by 12,000 measurement sweeps. All per-sweep data, starting preparations and ending configurations are retained. These are new executions and independently written diagnostics, not a second implementation of the gravitational dynamics.','',
'The geometry self-tests reproduced, including the flat-start deficit angles and incremental/full recomputation check. The statistical audit used ordinary split R-hat across all six chains at each parameter pair, with a prespecified screening threshold of 1.05, and batch-mean, drift and start-dependence diagnostics. Passing this screen would not prove equilibrium; failing it prevents treating the pooled chains as a settled ensemble. Curvature means below are averages over each start class and are descriptive when the chains have not mixed.','',
'| L | k | Flat-start mean R | Rough-start mean R | Split R-hat | Screen |',
'|---|---|---|---|---|---|']
for g in groups:
 lines.append(f'| {g["L"]} | {g["k"]} | {g["by_start"]["flat"]["mean_R"]:.4f} | {g["by_start"]["rough"]["mean_R"]:.4f} | {g["split_rhat"]:.3f} | {"passes" if g["passes_rhat_screen"] else "fails"} |')
fails=sum(not g['passes_rhat_screen'] for g in groups)
fail_label='All four' if fails==4 else str(fails)+' of four'
scale_error=max(abs(x['scaling_identity_ratio']-1) for x in runs)*100
lines+=['',f'{fail_label} parameter groups fail the prespecified R-hat screen. The table and trajectories distinguish initial-state dependence from an established equilibrium phase. Negative mean curvature after a flat start is insufficient to demonstrate that the ensemble has equilibrated. These finite-run tests neither locate a critical point nor establish the nonexistence of an equilibrium limit; they restrict what can be claimed from the sampled histories.', '', f'All 24 chains satisfy the scaling identity in their measurement averages to within {np.ceil(scale_error*100)/100:.2f} percent, despite the failed mixing screen. This separates a successful overall-scale consistency check from equilibration of the geometry. It provides a concrete example of an implementation check that cannot establish the sampled phase.', '',
'![Regge audit trajectories. Each curve averages 100 sweeps for display; all diagnostics use the unaveraged measurement histories. Gray shading marks equilibration sweeps.](regge_audit.png)','',
'The next tests are normalizability of the specified ensemble and convergence from both start classes at longer durations, followed by finite-size scaling and curvature correlations only where convergence is supported. A Newtonian potential and a Lorentzian matter-coupled theory require separate tests. The original histories remain unavailable; this audit supplies its own complete histories rather than claiming to reconstruct them.']
(D/'paper_audit_section.md').write_text('\n'.join(lines)+'\n')
report=['A42 T19 REPRODUCIBILITY AUDIT — 31 October 2026, Los Angeles','',f'24 main chains; 360,000 main sweeps, including 288,000 measurement sweeps. Six preparation chains add 18,000 sweeps. Timing pilot adds 200, excluded from inference.',f'R-hat failures: {fails}/4 groups.',f'Recomputed specification compression: {bits} bits. This is not a complete executable ToE.','', '\n'.join(lines),'', 'PER-CHAIN DIAGNOSTICS (SE is maximum over the retained batch sizes; no calibrated significance claim):']
for x in runs:
 report.append(f'{x["tag"]}: R={x["mean_R"]:.8g}, batch-SE={x["batch_se"]:.5g}, ESS={x["ess_R"]:.3g}, half-drift z={x["drift_z_descriptive"]:.4g}, scale-identity ratio={x["scaling_identity_ratio"]:.7g}, mean deficit spread={x["mean_deficit_std"]:.5g}')
report+=['','Interpretation: the scale identity is a useful implementation diagnostic, but satisfying it does not imply mixing of geometrical shapes. A global scaling move can equilibrate an overall size while shape variables remain trapped. This is a possible explanation to test, not a claim of a new physical law.','The supplied paper, code, inputs and transcripts are retained under T18_ORIGINAL for provenance; regenerable bulk outputs are omitted from this lean distribution. Statements in those historical files are not all endorsed by the revised paper. The T19 paper supersedes the T18 paper for this package.','The title and specific overclaims were revised; exact old/new replacements are in PAPER/CHANGES.json and PAPER/CHANGES.md. No fundamental physical update law was changed by this audit.']
report += ['', 'Data recovery: one output artifact was truncated despite a successful process return. The exact seeded job was replayed; all 3,440,640 bytes of the partial history matched the replay prefix. The complete 15,000-sweep replay replaces only that artifact and is not an independent chain. The partial output and replay are retained in the full T19 archive; they are omitted from this lean distribution.']
(root/'AUDIT_REPORT.txt').write_text('\n'.join(report)+'\n')
print('Generated report/section/figure; failed groups',fails,'bits',bits)
