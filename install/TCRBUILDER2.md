# Installing TCRBuilder2 (Stage B)

TCRBuilder2 turns the stitched paired sequences from Stage A into PDB structures
for Stage C. It ships inside **ImmuneBuilder** (Abanades et al., 2023).

These steps were tested in October 2026.

---

## Dependencies ImmuneBuilder does not declare

ImmuneBuilder's package metadata lists:

```
numpy, scipy>=1.6, einops>=0.3, torch>=1.8, requests
```

but importing it also needs **`pdbfixer`** and **`anarci`**. After
`pip install ImmuneBuilder` alone:

```python
>>> from ImmuneBuilder import TCRBuilder2
ModuleNotFoundError: No module named 'pdbfixer'
# fix that, and then:
ModuleNotFoundError: No module named 'anarci'
```

`pdbfixer` installs from PyPI. ANARCI does not, which is why the pipeline uses
conda environments.

---

## Recommended: conda

ANARCI needs a set of HMMs built from IMGT germline data. Bioconda ships them
prebuilt; installing from source builds them for you, which needs network access
to `imgt.org` at install time (see below).

```bash
conda create -n myenv2TCR python=3.10 -y
conda activate myenv2TCR

# ANARCI and HMMER together, with the HMMs already built
conda install -c bioconda -c conda-forge anarci hmmer -y

# the structure refinement stack
conda install -c conda-forge openmm pdbfixer -y

# TCRBuilder2 itself
pip install ImmuneBuilder
```

Then set the path in `config.sh`:

```bash
PY_TCRBUILDER="/path/to/envs/myenv2TCR/bin/python"
```

On ARC, load HMMER as a module instead of installing it — `run_B_slurm.sh`
already does `module load HMMER/3.3.2-gompi-2022a`.

## Verify

```bash
$PY_TCRBUILDER -c "from ImmuneBuilder import TCRBuilder2; print('ok')"
$PY_TCRBUILDER -c "import anarci, pdbfixer, openmm; print('deps ok')"
which hmmscan ANARCI
```

All four should succeed; the first exercises the whole import chain.

---

## Pip-only, if you cannot use conda

`pdbfixer` and `openmm` are on PyPI (tested with openmm 8.6.1):

```bash
pip install ImmuneBuilder pdbfixer openmm
```

ANARCI is not on PyPI. Install it from source:

```bash
sudo apt-get install hmmer          # or: brew install hmmer
pip install biopython
git clone https://github.com/oxpig/ANARCI
cd ANARCI && python setup.py install
```

`setup.py` downloads germline data from IMGT/GENE-DB and builds the HMMs during
installation, so it needs access to `imgt.org`. Without it, the install stops
with:

```
HMM file ./HMMs/ALL.hmm not found (nor an .h3m binary of it)
error: [Errno 2] No such file or directory: 'curated_alignments/germlines.py'
```

In that case use the conda route.

Requirements per ANARCI's own `INSTALL`: **HMMER3 ≥ 3.1** and **biopython**.

---

## Model weights: download once, before the array job

TCRBuilder2 fetches four model weight files from Zenodo on first use, into
`weights_dir`, which defaults to a `trained_model/` folder *inside the installed
package*. Two sets exist:

| weights | Zenodo record | selected by |
|---|---|---|
| TCRBuilder2**+** | `10892159` | `use_TCRBuilder2_PLUS_weights=True` *(the default)* |
| TCRBuilder2 | `7258553` | `use_TCRBuilder2_PLUS_weights=False` |

The paper uses the default, the TCRBuilder2+ weights (Quast et al., 2025): in
TCRBuilder2+ output the PDB
B-factor column holds an ensemble-disagreement estimate (§2.1).

Download the weights before launching Stage B. `run_B_slurm.sh` is a 100-task
array, and without the weights every task would download the same four files
into the same directory at once. Download them once from a login node:

```bash
$PY_TCRBUILDER -c "from ImmuneBuilder import TCRBuilder2; TCRBuilder2()"
```

That downloads the weights and exits. To keep them outside site-packages, pass
`weights_dir` explicitly.

Every run also prints this warning, which can be ignored:

```
FutureWarning: You are using `torch.load` with `weights_only=False` ...
```

Pinning the `torch` version keeps this behaviour fixed across releases.

---

## Input and output contract

**Input.** A FASTA per receptor with exactly two records, **named `B` and `A`**
(beta and alpha) — `scripts/stage1_stitchr.py` writes these. The CLI takes
sequences directly too:

```bash
ImmuneBuilder2 --help      # installed entry point
python -c "from ImmuneBuilder import TCRBuilder2; \
           TCRBuilder2().predict({'B':'<beta seq>','A':'<alpha seq>'})"
```

**Output.** `to_pdb(..., chain_ids="BA")` — so the PDB chains are labelled
**`B`** and **`A`**, IMGT-numbered (`numbering_scheme='imgt'` by default).

This is what SPACE2TCR expects: `get_antibody()` checks that each chain is `A`
or `B` and offsets the alpha chain by +128, so no renaming or renumbering is
needed between Stage B and Stage C. (The check's error message still refers to
heavy and light chains, from the antibody original; it means A and B. See
`SPACE2TCR/_upstream/NOTICE.md`.)

---

## Performance

About one minute per paired receptor on CPU and roughly five seconds on GPU,
against ~30 minutes for AlphaFold-Multimer.

For the 4,779-receptor cohort that is ~80 CPU-hours, hence the 100-task array
with a 48-hour limit. Budget disk too: one PDB per receptor, a few
GB for this cohort and tens of GB for a large one.

---

## Troubleshooting

| Symptom | Cause |
|---|---|
| `No module named 'pdbfixer'` | `pip install pdbfixer openmm`, or use conda |
| `No module named 'anarci'` | ANARCI is not on PyPI; install with conda or from source |
| `HMM file ./HMMs/ALL.hmm not found` during ANARCI install | the IMGT download failed; use the conda route |
| `hmmscan: command not found` | HMMER missing; `module load HMMER/3.3.2-gompi-2022a` on ARC |
| Stage B array tasks all stall at the start | the tasks are all downloading the weights; download them first |
| `ERROR: <model> not downloaded or corrupted` | partial Zenodo download; delete the file in `weights_dir` and re-fetch |
| Stage C rejects the PDBs | chains are not `A`/`B`, or not IMGT-numbered — check Stage A's FASTA record names |

---

## References

- Abanades, Wong, Boyles, Georges, Bujotzek & Deane (2023). *ImmuneBuilder:
  Deep-Learning models for predicting the structures of immune proteins.*
  Communications Biology 6, 575.
- Quast, Abanades, Guloglu, Karuppiah, Harper, Raybould & Deane (2025).
  *T-cell receptor structures and predictive models reveal comparable alpha and
  beta chain structural diversity despite differing genetic complexity.*
  Communications Biology 8, 362. <https://doi.org/10.1038/s42003-025-07708-6>
  (TCRBuilder2+, the weights used in the paper.)
- <https://github.com/brennanaba/ImmuneBuilder>
- <https://github.com/oxpig/ANARCI>
