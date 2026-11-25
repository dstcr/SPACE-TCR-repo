# Provenance and licensing of `SPACE2TCR`

`SPACE2TCR/` is a TCR adaptation of **SPACE2** (Structural Profiling of
Antibodies to Cluster by Epitope 2; Spoendlin et al., 2023), distributed under
its BSD 3-Clause licence. `_upstream/LICENSE` is the upstream licence text,
© 2022 Brennan Abanades Kenyon, from <https://github.com/fspoendlin/SPACE2>; the
fork is based on upstream `version='1.0.1'`.

`_upstream/setup.py` keeps the upstream metadata:

```
name='SPACE2TCR'            version='1.0.1'
description='Structural clustering of antibody CDR loops'
license='BSD 3-clause license'
author=['Fabian Spoendlin', 'Brennan Abanades']
maintainer='Fabian Spoendlin'   <fabian.spoendlin@stats.ox.ac.uk>
```

Upstream: <https://github.com/fspoendlin/SPACE2> and
<https://github.com/oxpig/SPACE2>.

`_upstream/README.md` is the upstream SPACE2 README. It describes antibodies
(chains `H` and `L`) and is kept for provenance; for this pipeline, follow
[`INSTALL.md`](../../INSTALL.md) instead.

## Changes for TCRs

The region definitions keep upstream's residue numbering and rename the regions
from heavy/light to beta/alpha:

| upstream (antibody) | fork (TCR) | residues |
|---|---|---|
| `CDRH1` `CDRH2` `CDRH3` | `CDRB1` `CDRB2` `CDRB3` | 17, 12, 13 |
| `CDRL1` `CDRL2` `CDRL3` | `CDRA1` `CDRA2` `CDRA3` | 17, 15, 13 |
| `fwH` | `fwB` | 74 |
| `fwL` | `fwA` | 67 |

IMGT unique numbering is defined across both IG and TR V-domains, so the CDR
boundaries apply to TCRs unchanged. As upstream, the second chain is offset by
+128 so both chains share one residue-number array:

```
reg_def["CDRB1"]  24-40        reg_def["CDRA1"]  152-168
reg_def["CDRB2"]  55-66        reg_def["CDRA2"]  183-197
reg_def["CDRB3"]  105-117      reg_def["CDRA3"]  233-245
reg_def["fwB"]    (74 positions)   reg_def["fwA"]    (67 positions)

reg_def["CDR_all"] = CDR_B_all + CDR_A_all     6 parts, 87 residues
reg_def["fw_all"]  = [fwB, fwA]                141 anchor residues
```

`get_antibody()` applies the offset while parsing and requires chains labelled
`A` and `B`. Residues numbered above 128 within a chain are not used.

Two upstream names remain from the antibody original: the chain check's error
message (`util.py`, line 70) mentions heavy and light chains, and
`fw = ["fwH", "fwL"]` (`util.py`, line 9) is unused.

## Where the rest of the method lives

Two parts of SPACE-TCR are implemented in the pipeline scripts rather than in
this package:

- **The dihedral gate**: φ/ψ extraction is in `scripts/extract_dihedrals.py`,
  and the circular-distance filter in `scripts/cluster_one_bin.py` (with
  `scripts/dihedral_util.py`).
- **Clustering**: the RMSD-threshold graph and its connected components, and
  the `--max-exact-bin` switch to greedy representatives for very large bins,
  are in `scripts/cluster_one_bin.py`. This package's `greedy_clustering.py` and
  `exhaustive_clustering.py` differ from upstream only in import names and in
  `length_tolerance` taking an array.

## Not included

- `__pycache__/`: the compiled and numba cache, regenerated on first run.
- `build/` and `SPACE2TCR.egg-info/`: build output. `requires.txt` is kept as
  the record of the declared dependencies.

## Runtime dependencies

From `_upstream/requires.txt`:

```
numba>=0.51.0
joblib>=1.0.0
scipy
pandas
scikit-learn>=1.4.0
```

`get_residues` and `compare_CDRs_for_cluster` are compiled with
`@nb.njit(cache=True, fastmath=True)`, so the cache depends on the numba and
Python versions.

## Checking the install

`tests/verify_space2tcr.py` checks the names the pipeline imports and the region
definitions, and computes RMSDs, including that the distance is unchanged by
rigid-body rotation and translation:

```bash
PYTHONPATH=. python3 tests/verify_space2tcr.py .
```
