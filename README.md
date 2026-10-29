# secretscores

Calculate polygenic scores from PLINK binary genotypes (`.bed/.bim/.fam`) and
[PGS Catalog](https://www.pgscatalog.org/) scoring files.

Version 1.0.0 — (c) 2026 Wei-Min Chen. See the license notice at the top of each source file.

## Build

You need a C++ compiler and zlib (`zlib1g-dev` on Debian/Ubuntu). OpenMP is optional but
gives multithreading.

```bash
make
```

This creates the `secretscores` executable. Without OpenMP, use `make CXXFLAGS=-O2`.

## Quick start

The `example/` folder has a script that writes a tiny dataset (4 people, 5 SNPs) and a model.

```bash
cd example
python3 make_example.py
../secretscores -b ex --readmodel model.txt --prefix out
```

`-b ex` is the PLINK file prefix: `ex.bed`, `ex.bim` and `ex.fam` must exist. Output:

| File | Contents |
|---|---|
| `outgrs.txt` | Single model only: `FID IID InfoSNP Score`. `InfoSNP` is the fraction of the model's SNPs with a genotype for that person. |
| `out_scores.txt` | Tab-delimited `IID` plus one column per model, people in `.fam` order. |
| `out_scores.coverage` | One line per model: SNPs in the file, SNPs found, SNPs used, coverage, and why SNPs were dropped. |

Expected `out_scores.txt`:

```
IID     model
ind1    0.33
ind2    0.21
ind3    0.25
ind4    0.38
```

The score is the sum over SNPs of `effect_weight × copies of the effect allele`. A missing
genotype contributes 0 (which is why `ind4` has `InfoSNP` 0.8).

## Model files

**Text model:** a PGS Catalog *harmonized* scoring file (plain or `.gz`). The header needs the
columns `effect_allele`, `effect_weight`, `hm_chr`, `hm_pos`, and either `other_allele` or
`hm_inferOtherAllele`. Lines starting with `#` are ignored.

**How SNPs are matched:** by chromosome and position, then by alleles. If the effect/other
alleles equal the data's A1/A2 the action is `OK`; if reversed it is `SWITCHED`. Alleles are
compared as given, with no strand flipping, so SNPs with mismatched alleles are dropped (see
the `AlleleMismatch` column of the coverage file). Only autosomes and XY are scored.

**Binary model:** `--makebin` converts text models to compact `.bin` files (no genotypes
needed). Use a `.bin` file with `--readmodel` exactly like a text model. Only
single-letter A/C/G/T alleles can be stored.

```bash
./secretscores --makebin --readmodel model.txt     # writes model.bin
./secretscores -b ex --readmodel model.bin
```

## Scoring many models

Give several files separated by commas:

```bash
./secretscores -b ex --readmodel a.txt,b.txt.gz,c.bin
```

or list them in a file, one model per line, with an optional second column naming the score
(lines starting with `#` are skipped):

```
# models.list
PGS000001.txt.gz   height
PGS000002.txt.gz   bmi
```

```bash
./secretscores -b ex --readlist models.list --prefix mine
```

Without a name, the score is named after the file (no directory, `.gz` or extension). Score
names must be unique.

## Options

| Option | Meaning |
|---|---|
| `-b <prefix>` | PLINK binary file prefix (`.bed/.bim/.fam`, `.bim`/`.fam` may be `.gz`) |
| `--readmodel <files>` | Model file(s), comma-separated |
| `--readlist <file>` | File listing model files (and optional names) |
| `--makebin` | Convert text models to `.bin`; no genotypes are read |
| `--cpus <n>` | Threads (default 0 = automatic, half the available cores) |
| `--prefix <name>` | Output file prefix (default `secretscores`) |

## Tips

- **Speed:** individual-major `.bed` files are faster than the default SNP-major ones. Convert once with
  `plink2 --bfile ex --export ind-major-bed --out ex_indmajor` and use `-b ex_indmajor`.
- **Memory:** about half of the available RAM is used to hold genotypes; many models are
  scored in batches.
- **Check coverage:** look at `*_scores.coverage` before trusting a score. Low coverage
  usually means a genome-build mismatch between the model and your data.
