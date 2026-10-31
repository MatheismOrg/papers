FIXED-PHYSICAL-RADIUS CLUSTER STRESS CHECK

Purpose: determine whether the prior R500 radius-augmented fit remains accurate
away from the radius whose definition algebraically links radius to target mass.

Run:
  python stress_test.py --data-root /path/to/dm_fit
Requires numpy, scipy, and astropy. Inputs are the original dm_fit data and result
files. All outputs are written beside this script. No network access is needed.

Seven previously examined X-COP clusters are sampled at 500, 750, and 1000 kpc.
All 21 locations lie within the supplied gas, stellar, and hydrostatic mass grids.
No extrapolation is used. These are 21 correlated, model-derived profile samples,
not 21 independent raw acceleration observations. The seven clusters are not new
external validation data.

Primary target: forward hydrostatic mass profile M_FORW. NFW and Einasto target
profiles are sensitivity cases. The forward profile is also an HSE reconstruction,
not a gravity-independent measurement.

Baryons: spherical gas mass plus 0.75 times projected enclosed stellar mass.
The 0.75 conversion was motivated at R500 in the source study; it is an explicitly
approximate treatment at the fixed physical radii here. We also compute no-star
and full-projected-star cases. These are sensitivity cases, not confidence bounds
on all stellar-systematic uncertainty. Full projected stellar mass is 3.7-9.8%
of gas mass across these locations, so its influence is modest.

All fitted cluster tests leave a whole cluster out, retaining three radii from
each of the other six. Every system has total training weight one. The joint
quadratic also trains on all 149 SPARC galaxies; this test evaluates cluster
prediction only and must not be reported as held-out galaxy validation.

The frozen previous R500 polynomial already saw all seven clusters during its
original fit at R500. Its transfer to these fixed radii is a within-cluster radial
stress test, not a held-out-cluster validation.

Primary-case median absolute acceleration errors:
  Galaxy-calibrated RAR, frozen:                  56.16%
  Previous R500 radius-quadratic, frozen:         39.02%
  Cluster-calibrated one-a0 RAR, leave-one-out:    15.98%
  Joint galaxy/cluster quadratic, leave-one-out:  21.02%
  Cluster radius-only power law, leave-one-out:   17.96%

The full-sample cluster RAR scale is 7.1591e-10 m/s^2, versus the prior galaxy
scale of 1.1603e-10 m/s^2. This is a class-specific calibration, not one universal
coefficient set. The no-star/full-projected-star cases change the cluster RAR
error to 16.99%/15.57%. With NFW/Einasto targets and the middle stellar case, it
is 13.61%/17.29%. The model ranking is consequently sensitive to mass inference.

The small earlier R500 radius-fit error does not persist at these physical radii.
This supports treating that earlier score as a limited, target-dependent result
rather than evidence of a universal acceleration formula. A gas-informed RAR can
still be useful as a piecemeal description; the radius-only null is competitive.

All model metrics in stress_test_results.json were independently recomputed from
fixed_radius_predictions.csv. Formula parameters, source hashes, per-cluster
metrics, and all nine target/stellar sensitivity cases are supplied in the JSON.
