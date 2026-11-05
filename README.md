# secretscores

Calculate polygenic scores from PLINK binary genotypes (`.bed/.bim/.fam`) and
[PGS Catalog](https://www.pgscatalog.org/) scoring files.

Version 1.0.2 — (c) 2026 Wei-Min Chen. See the license notice at the top of each source file.

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
| `out_scores.coverage` | One line per model: SNPs in the file, SNPs found, SNPs used, coverage, why SNPs were dropped, the percentage of genotypes imputed, and the trait. |

Expected `out_scores.txt`:

```
IID     model
ind1    0.33
ind2    0.21
ind3    0.25
ind4    0.38
```

The score is the sum over SNPs of `effect_weight × copies of the effect allele`. A missing
genotype contributes 0 unless allele frequencies are available (see *Missing genotypes*), which is
why `ind4` has `InfoSNP` 0.8.

## Model files

**Text model:** a PGS Catalog *harmonized* scoring file (plain or `.gz`). The header needs the
columns `effect_allele`, `effect_weight`, `hm_chr`, `hm_pos`, and either `other_allele` or
`hm_inferOtherAllele`. The optional column `allelefrequency_effect` (the frequency of the effect
allele) is used for mean imputation. Lines starting with `#` are ignored.

**How SNPs are matched:** by chromosome and position, then by alleles. If the effect/other
alleles equal the data's A1/A2 the action is `OK`; if reversed it is `SWITCHED`. Alleles are
compared as given, with no strand flipping, so SNPs with mismatched alleles are dropped (see
the `AlleleMismatch` column of the coverage file). Autosomes and XY are scored; the X chromosome
is scored only with `--allowX 1` or `--allowX 2`.

**Binary model:** `--makebin` converts text models to compact `.bin` files (no genotypes
needed). Use a `.bin` file with `--readmodel` exactly like a text model. Only
single-letter A/C/G/T alleles can be stored. If the text model has the column
`allelefrequency_effect`, the `.bin` file (version 4) stores the frequency of the effect allele of
each SNP as a 4-byte float, right after the weights; otherwise a version 3 file is written.
`--viewbin` prints a `.bin` model as text.

```bash
./secretscores --makebin --readmodel model.txt     # writes model.bin
./secretscores -b ex --readmodel model.bin
```

## Missing genotypes

By default a missing genotype is replaced by its mean: 2 × the frequency of the effect allele
(the mean genotype dosage), times the weight. The frequency comes from

1. the `--readfreq` file, if given (a PLINK `.frq` file with columns `SNP A1 MAF`, or a PLINK 2
   `.afreq` file with columns `ID ALT ALT_FREQS`, plain or `.gz`; SNPs are matched by ID and
   allele). A SNP that the file has no frequency for is still scored, without imputation, and a
   warning is printed;
2. otherwise the `allelefrequency_effect` column of the model.

A SNP with no frequency contributes 0 for a missing genotype. Use `--no-mean-imputation` to
count every missing genotype as 0 (for imputed data, which has no missing genotypes). The
`PercentImputed` column of the coverage file is the percentage of the genotypes (SNPs used ×
individuals) that were imputed.

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
| **Input Files** | |
| `--readmodel <files>` | Model file(s), comma-separated |
| `--readlist <file>` | File listing model files (and optional names) |
| `--readfreq <file>` | Allele frequencies of the genotype data (PLINK `.frq` or PLINK 2 `.afreq`), used for mean imputation instead of the model's `allelefrequency_effect` |
| **Score Action** | |
| `--score` | Compute the scores (the default action) |
| `--makebin` | Convert text models to `.bin`; no genotypes are read |
| `--viewbin` | Print a `.bin` model as text |
| **Score Parameter** | |
| `--no-mean-imputation` | Count a missing genotype as 0, do not use allele frequencies |
| `--allowX <0/1/2>` | X chromosome: 0 = not scored (default), 1 = males 0/1 and females 0/1/2, 2 = males and females 0/1/2 |
| **Other** | |
| `--cpus <n>` | Threads (default 0 = automatic, half the available cores) |
| `--prefix <name>` | Output file prefix (default `secretscores`) |

## Tips

- **Speed:** individual-major `.bed` files are faster than the default SNP-major ones. Convert once with
  `plink2 --bfile ex --export ind-major-bed --out ex_indmajor` and use `-b ex_indmajor`.
- **Memory:** about half of the available RAM is used to hold genotypes; many models are
  scored in batches of at most 1000 models and fewer than 9 million different SNPs.
- **Check coverage:** look at `*_scores.coverage` before trusting a score. Low coverage
  usually means a genome-build mismatch between the model and your data.
