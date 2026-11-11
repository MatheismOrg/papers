"""Print 500-sweep segment means of a Regge history (R, deficit spread, max q, and sumR2 when present)."""
import sys, numpy as np
for fn in sys.argv[1:]:
    try: a=np.loadtxt(fn)
    except Exception as e: print(fn, 'unreadable', e); continue
    if a.ndim<2 or len(a)==0: print(fn,'empty'); continue
    hdr=[l.strip() for l in open(fn) if l.startswith('# a=') or l.startswith('# loaded') or l.startswith('# flat')]
    print('==',fn,'rows',len(a),'|',' | '.join(hdr)[:200])
    for lo in range(0,len(a),500):
        s=a[lo:lo+500]
        extra=f" sumR2/NT={s[:,-1].mean()/(50*int(round((a.shape[0] and 1)*0+1))):.4g}" if a.shape[1]>27 else ""
        print(f"  {lo:5d}-{lo+len(s)-1:5d} R={s[:,4].mean():9.4f} std_d={s[:,21].mean():.4f} maxq={s[:,9].mean():9.3f} meanq={s[:,5].mean():.4f} acc={s[:,6].mean():.3f}"+(f" sumR2={s[:,-1].mean():.5g}" if a.shape[1]>27 else ""))
