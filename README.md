# PANTERA I

Data and code for the first paper of PANTERA, a high-resolution spectroscopic survey
of the most metal-rich stars in the solar neighborhood.

Paper Link: https://arxiv.org/abs/2607.27328

The sample is 56 stars selected from the Gaia DR3 XP metallicities of Andrae et al.
(2023) and observed with PEPSI on the Large Binocular Telescope, the Levy spectrograph
on the Automated Planet Finder, and HIRES on Keck I. Iron abundances come from an
equivalent-width analysis with the temperature fixed from photometry and the surface
gravity from the Gaia parallax.

## Tables

| File | Rows | Contents |
| --- | --- | --- |
| `data/table_targets.csv` | 56 | Positions, G magnitude, XP metallicity, and the instrument used |
| `data/table_results.csv` | 62 | Adopted parameters and iron abundances, one row per spectrum |
| `data/table_equivalent_widths.csv` | 4700 | Line-by-line equivalent widths with the atomic data |
| `data/table_kinematics.csv` | 56 | Radial velocities, space motions, actions, and orbital quantities |
| `data/feh_ewcaps.csv` | 248 | Iron abundances at four equivalent-width caps, one row per spectrum and cap |
| `data/feh_ewcaps_star.csv` | 224 | The same, averaged over the spectra of each star |

Stars are identified by their Gaia DR3 source identifier in every table. Six stars were
observed with more than one spectrograph, which is why the results table has 62 rows for
56 stars.

Column names follow the paper. `teff_K` and `logg_dex` are the adopted temperature and
surface gravity, `vmic_kms` is the microturbulence from the reduced-equivalent-width
balance, `FeH_dex` is the abundance from the Fe I lines with its uncertainty
`eFeH_dex`, `FeII_dex` is the Fe II abundance, `dion_dex` is the difference between the
two, and `nFeI` is the number of Fe I lines. `MH_XP_dex` is the Gaia XP value the star
was selected on.

## The equivalent-width cap

The abundances in the paper use the Fe I lines with equivalent widths below 120 mA.
Lines stronger than that are affected by blends in their wings, and because the strong
lines are also the ones that set the microturbulence, including them biases the result
low. For the benchmark giant mu Leo, whose reference value is +0.25, truncating at 100,
120, and 140 mA returns +0.250, +0.231, and +0.172, and with no cap the same procedure
returns +0.079. Relative to an uncapped solution the adopted cap raises the giants by
about 0.14 dex, the subgiants by 0.07 dex, and the dwarfs by 0.05 dex.

`data/feh_ewcaps.csv` gives the abundance, the microturbulence, and the number of
retained lines for each spectrum at each of the four caps, so the effect of this choice
can be checked directly. The `w120` rows are the values adopted in the paper and match
`data/table_results.csv`. The `railed` column flags the spectra where the abundance
against reduced-equivalent-width slope does not cross zero inside the search range, so
the microturbulence sits at a boundary rather than at a balanced solution.

## Code

| File | Purpose |
| --- | --- |
| `code/measure_equivalent_widths.py` | Measures equivalent widths by fitting a Gaussian to the line core over a local linear continuum |
| `code/solve_abundances.jl` | Solves the microturbulence by reduced-equivalent-width balance and returns the iron abundance |
| `code/solve_abundances_ewcaps.jl` | Rebuilds `feh_ewcaps.csv` from the released tables |

The abundance step uses Korg (Wheeler et al. 2023) with MARCS model atmospheres. To
rebuild the cap table from the released data alone, run

```
julia code/solve_abundances_ewcaps.jl
```

Paths at the top of that script point at the paper directory and need editing to match
where the tables sit.

## Reduced spectra

The reduced spectra are not in this repository because of their size. They are available
from the corresponding author on request.
