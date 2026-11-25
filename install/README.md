# install/

Four installers. [`../install.sh`](../install.sh) picks the right ones and
checks the result:

```bash
./install.sh                 # reproduce the paper (R only)
./install.sh --own-data      # + structural clustering
./install.sh --everything    # + the ML distillation
./install.sh --check         # check, change nothing
```

Each can also be used on its own: they are independent of one another, can be
re-run safely, take `--check`, and exit non-zero on failure.

| What | Run | Time | Needs |
|---|---|---|---|
| **R** — analysis and figures | `Rscript install/r_packages.R [--figures]` | ~3 min (+10 with `--figures`) | R ≥ 4.1 |
| **SPACE2-TCR** — Stages C & D | `bash install/space2tcr.sh` | ~2 min | Python 3.10–3.13 |
| **Distillation model** (`ml/`) | `bash install/distill.sh --venv .venv-ml` | ~5 min | Python 3.10–3.13 |
| **TCRBuilder2** — Stage B | see [`TCRBUILDER2.md`](TCRBUILDER2.md) | ~20 min | conda |

Tested versions and the full laptop-and-cluster walkthrough are in
[`../INSTALL.md`](../INSTALL.md).
