"""Negative control: overdensity radius already encodes the mass target.

Usage: python dm_audit_r500.py --fit-root tmp/dm_fit --output audit.json
Requires numpy and astropy. This is an identity check, not a DM fit.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fit-root', type=Path, required=True)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    # The current reproducibility workspace kept optional astropy here.
    if (args.fit_root / 'deps').is_dir():
        sys.path.insert(0, str(args.fit_root / 'deps'))
    from astropy.io import fits

    with (args.fit_root / 'results' / 'cluster_test.csv').open() as f:
        rows = list(csv.DictReader(f))
    radius = np.array([float(row['r500_kpc']) for row in rows])
    observed = np.array([float(row['gobs']) for row in rows])
    redshift = []
    for row in rows:
        name = row['cluster'] + '_fgas_profile.fits'
        with fits.open(args.fit_root / 'data' / 'xcop' / name) as hdus:
            redshift.append(hdus['FGAS'].header['REDSHIFT'])
    redshift = np.array(redshift)
    e2 = 0.3 * (1 + redshift) ** 3 + 0.7
    ratio = observed / radius
    # Geometric least squares, slope fixed to one, intercept from six others.
    pred_radius = np.array([
        radius[i] * np.exp(np.mean(np.log(np.delete(ratio, i))))
        for i in range(len(rows))
    ])
    pred_radius_redshift = np.array([
        radius[i] * e2[i]
        * np.exp(np.mean(np.log(np.delete(ratio / e2, i))))
        for i in range(len(rows))
    ])
    kpc_m = 3.085677581491367e19
    h0 = 70 * 1000 / (1000 * kpc_m)
    pred_identity = 250 * h0 ** 2 * e2 * radius * kpc_m

    def metrics(pred):
        return {
            'median_absolute_acceleration_fraction': float(
                np.median(abs(pred / observed - 1))),
            'maximum_absolute_acceleration_fraction': float(
                np.max(abs(pred / observed - 1))),
            'rms_log10_acceleration': float(
                np.sqrt(np.mean(np.log10(pred / observed) ** 2))),
            'predictions_m_s2': pred.tolist(),
        }

    result = {
        'purpose': 'Negative control; target construction, not a physical fit.',
        'identity': 'g500=250 H(z)^2 R500',
        'cosmology': 'H0=70 km/s/Mpc, Omega_m=0.3, Omega_Lambda=0.7',
        'cluster_names': [row['cluster'] for row in rows],
        'redshifts': redshift.tolist(),
        'r500_kpc': radius.tolist(),
        'gobs_m_s2': observed.tolist(),
        'leave_one_out_radius_only': metrics(pred_radius),
        'leave_one_out_radius_redshift': metrics(pred_radius_redshift),
        'no_fit_overdensity_identity': metrics(pred_identity),
    }
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
