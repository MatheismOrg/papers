from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent/'results'
s=json.loads((P/'external_validation_results.json').read_text());a0=s['frozen_a0_m_s2']
rar=lambda g:g/(-np.expm1(-np.sqrt(g/a0)))
e=list(csv.DictReader((P/'etg_external_pairs.csv').open()))
l=list(csv.DictReader((P/'lensing_external_pairs.csv').open()))
gr=list(csv.DictReader((P/'group_velocity_dispersion_test.csv').open()))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,3,figsize=(13.8,4.7),layout='constrained')
ax=axs[0]
for tr,label,col,mark in [('inner_JAM','Inner stellar JAM','#b97842','s'),('outer_HI','Outer HI','#24617e','o')]:
 r=[q for q in e if q['tracer']==tr];x=np.array([float(q['gbar_m_s2']) for q in r]);y=np.array([float(q['gobs_m_s2']) for q in r]);ax.scatter(x,y,c=col,marker=mark,label=label,s=30)
xx=np.geomspace(2e-12,2e-8,300);ax.plot(xx,rar(xx),color='black',label='Frozen SPARC RAR');ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Baryonic acceleration (m/s²)');ax.set_ylabel('Observed acceleration estimate (m/s²)');ax.set_title('16 early types: 32 published pairs');ax.legend(fontsize=8)
ax=axs[1]
for sample,label,col,mark in [('GAMA_isolated','GAMA isolated stack','#24617e','o'),('KiDS_isolated','KiDS isolated stack','#b97842','s')]:
 r=[q for q in l if q['system']==sample];x=np.array([float(q['gbar_m_s2']) for q in r]);y=np.array([float(q['gobs_m_s2']) for q in r]);dy=np.array([float(q['gobs_error_m_s2']) for q in r]);ax.errorbar(x,y,dy,fmt=mark,ms=3.5,color=col,label=label,alpha=.8)
xx=np.geomspace(7e-16,8e-12,300);ax.plot(xx,rar(xx),color='black');ax.axvline(1e-13,color='#777777',ls=':',lw=1);ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Baryonic acceleration (m/s²)');ax.set_title('Correlated lensing bins; SIS conversion');ax.legend(fontsize=8)
ax.text(.02,.02,'KiDS isolation uncertain left of dotted line',transform=ax.transAxes,fontsize=8)
ax=axs[2]
for excluded,label,col,mark in [('False','55 source-selected rows','#24617e','o'),('True','11 source-excluded rows','#b97842','x')]:
 r=[q for q in gr if q['source_excluded']==excluded];x=[float(q['sigma_predicted_frozen_a0_km_s']) for q in r];y=[float(q['sigma_observed_km_s']) for q in r];ax.scatter(x,y,label=label,c=col,marker=mark,s=25)
ax.plot([20,450],[20,450],color='black');ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(20,450);ax.set_ylim(20,450);ax.set_xlabel('Predicted dispersion (km/s)');ax.set_ylabel('Observed dispersion (km/s)');ax.set_title('Groups: separate global virial prediction');ax.legend(fontsize=8)
for ax in axs:ax.grid(alpha=.15)
fig.suptitle('Frozen galaxy calibration transferred to different observables',fontsize=13)
fig.savefig(P/'external_observables.png',dpi=180,bbox_inches='tight')
