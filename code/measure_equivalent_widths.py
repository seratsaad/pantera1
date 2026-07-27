"""
Step 26 — Measure equivalent widths (mA) of the Fe I/II lines in a normalised
spectrum by Gaussian fitting with a local linear continuum.
  Usage: python 26_measure_ews.py <norm_npz> <out_csv>
Outputs: wl_air, species, Xexc, loggf, ew_mA  (only well-measured lines)
"""
import sys, numpy as np, pandas as pd
from scipy.optimize import curve_fit

ROOT = '/Users/saad.104/Downloads/pantera'
npzf, outcsv = sys.argv[1], sys.argv[2]
lines = pd.read_csv(f'{ROOT}/results/ew_fe_lines.csv')

d = np.load(npzf); w = d['wav']; f = d['flux']
o = np.argsort(w); w, f = w[o], f[o]
ok = np.isfinite(f); w, f = w[ok], f[ok]

def gauss(x, a, mu, sig):
    return -a*np.exp(-(x-mu)**2/(2*sig**2))

rows = []
for _, ln in lines.iterrows():
    wc = ln['wl_air']
    # continuum from the WINGS (upper envelope), core fit separately
    wing = (np.abs(w-wc) > 0.45) & (np.abs(w-wc) < 1.6)
    core = np.abs(w-wc) < 0.40
    if wing.sum() < 8 or core.sum() < 5: continue
    # robust linear continuum through the high points of the wings
    xw, yw = w[wing], f[wing]
    hi = yw > np.percentile(yw, 60)
    try:
        c1, c0 = np.polyfit(xw[hi]-wc, yw[hi], 1)
    except Exception:
        continue
    cont = lambda xx: c0 + c1*(xx-wc)
    xc, yc = w[core], f[core]
    norm = yc / cont(xc)                     # continuum-normalised core
    depth0 = max(1 - np.min(norm), 0.02)
    try:
        p, _ = curve_fit(lambda x,a,mu,sig: 1+gauss(x,a,mu,sig), xc, norm,
                         p0=[depth0, wc, 0.09],
                         bounds=([0.0, wc-0.12, 0.025],[1.1, wc+0.12, 0.35]),
                         maxfev=6000)
    except Exception:
        continue
    a, mu, sig = p
    ew = 1000.0 * a * sig * np.sqrt(2*np.pi)   # mA (core already normalised)
    model = 1 + gauss(xc, *p)
    rms = np.sqrt(np.nanmean((norm-model)**2))
    fwhm_kms = 2.355*sig/wc*3e5
    if not (5 < ew < 200): continue
    if not (2.0 < fwhm_kms < 30): continue
    if rms > 0.10: continue
    if a < 0.025: continue
    rows.append(dict(wl_air=round(wc,4), species=ln['species'], Xexc=ln['Xexc'],
                     loggf=ln['loggf'], ew_mA=round(ew,1)))

df = pd.DataFrame(rows)
df.to_csv(outcsv, index=False)
n1 = (df.species=='Fe I').sum(); n2 = (df.species=='Fe II').sum()
print(f'{npzf.split("/")[-1]}: measured {len(df)} lines (Fe I={n1}, Fe II={n2}) -> {outcsv}')
