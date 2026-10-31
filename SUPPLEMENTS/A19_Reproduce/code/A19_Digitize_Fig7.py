"""Digitize Brouwer et al. 2021 Fig. 7 (isolated red and blue KiDS-bright galaxies, stellar mass histograms)
from a screenshot of the A&A figure (aa40108-20-fig7.jpg). Axis calibration from the tick marks.
Output: counts per 0.05-dex interval (mean of the step curve over the interval) for red, blue and min(red, blue)."""
import sys, json, numpy as np
from PIL import Image
im = np.array(Image.open(sys.argv[1]).convert('RGB')).astype(int)
R, G, B = im[..., 0], im[..., 1], im[..., 2]
red = (R - B > 60) & (R - G > 60)
blue = (B - R > 45) & (R < 160)
X0, PXD = 135.5, (764.5 - 135.5) / 3.0      # x = 135.5 at log M = 8.0; 3 dex over 629 px
Y0, PXC = 892.5, (892.5 - 434.0) / 6000.0  # y = 892.5 at 0; 6000 at 434
def top(mask, x):
    col = mask[380:893, x]; ys = np.where(col)[0]
    if x < 200: ys = ys[ys + 380 > 520]      # skip the legend box
    return (Y0 - (ys.min() + 380)) / PXC if len(ys) else 0.0
xs = np.arange(137, 764)
lm = 8.0 + (xs - X0) / PXD
r = np.array([top(red, x) for x in xs]); b = np.array([top(blue, x) for x in xs])
r[r < 0] = 0; b[b < 0] = 0
m = np.minimum(r, b)
edges = np.round(np.arange(8.0, 11.0001, 0.05), 3)
out = []
for lo, hi in zip(edges[:-1], edges[1:]):
    k = (lm >= lo) & (lm < hi)
    out.append(dict(lo=float(lo), hi=float(hi), red=round(float(r[k].mean()), 1), blue=round(float(b[k].mean()), 1), matched=round(float(m[k].mean()), 1)))
json.dump(dict(source='Brouwer et al. 2021, A&A 650, A113, Fig. 7 (screenshot of aa40108-20-fig7.jpg)', calibration=dict(x0=X0, px_per_dex=PXD, y0=Y0, px_per_count=PXC), bins=out), open(sys.argv[2], 'w'), indent=1)
for o in out:
    if o['red'] > 20 or o['blue'] > 20: print(o)
