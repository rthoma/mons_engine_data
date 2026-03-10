# MONS v1.0.0 — Data Repository (internal)

Internal companion repository for MONS v1.0.0 that hosts large data files.

## Contents

- `data/` — large datasets and supporting files
- `MiniWECC/` — MiniWECC dataset/materials (as provided)

## Where to put this repo's folders

Place the contents of this data repo into your MONS code repo so the final layout looks like this:

```text
mons_v1p0p0/
├─ src/
│  └─ TestCases/
│     └─ MiniWECC/          (from this repo: MiniWECC/)
├─ tests/
│  └─ data/                 (from this repo: data/)
├─ run_mons_tests.m
└─ mons_v1p0p0_doc.pptx
```

## Notes
- This repo is for data only; the MONS engine/source code lives in the separate code repository.
- POC: Ryan Elliott, rtellio@sandia.gov
