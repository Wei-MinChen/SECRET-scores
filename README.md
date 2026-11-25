# secretscores

Calculate polygenic scores from PLINK binary genotypes (`.bed/.bim/.fam`), a VCF file (`--vcf`) or
PLINK 2 files (`--pfile`: `.pgen/.pvar/.psam`) and [PGS Catalog](https://www.pgscatalog.org/) scoring files.

Version 1.0.4 — (c) 2026 Wei-Min Chen. See the license notice at the top of each source file.

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
../secretscores --bfile ex --readmodel model.txt --prefix out
```

`--bfile ex` is the PLINK file prefix: `ex.bed`, `ex.bim` and `ex.fam` must exist. Output:

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

## Genotype input

Give exactly one of:

| Option | Input |
|---|---|
| `--bfile <prefix>` | PLINK 1 binary files `.bed/.bim/.fam` (`<prefix>` or `<prefix>.bed`) |
| `--vcf <file>` | A VCF, plain or gzipped (`.vcf`, `.vcf.gz`). BCF is not read. |
| `--pfile <prefix>` | PLINK 2 files `<prefix>.pgen`, `.pvar` (or `.pvar.gz`) and `.psam` |

Only hard calls are used (the `GT` field of a VCF, the genotype calls of a `.pgen`); dosages are
ignored. A record with several ALT alleles gives one SNP for each ALT allele (that allele against REF),
so a model SNP is matched to the right allele by position and alleles as usual. A VCF has no sex: its
samples are of unknown sex (scored as females at the X chromosome with `--allowX`), and their FID is 0
(`--keep`/`--remove` lines can hold the IID alone). A `.psam` gives the sex when it has a `SEX` column.
A haploid genotype is taken as homozygous, a half call (`./1`) as missing. A VCF is read through for
every batch of models; for many models or large data, a `.pgen` is much faster.

## ABO blood type

`--inferABO` infers the ABO type of every sample from a few SNPs of the genotype data (`--bfile`, `--vcf` or `--pfile`;
no model is needed) and writes `<prefix>_ABO.txt` (`IID`, `ABO`, `O_copies`, `A_copies`, `Flag`).

```bash
./secretscores --inferABO --vcf ex.vcf.gz --prefix out
```

The SNPs are found by their variant ID (rsID). A SNP not found by ID is looked for by chromosome 9 and position, in GRCh38 or GRCh37 (the build with more SNPs found), and it has to have the expected alleles (or their complement); the list printed says how each SNP was found. Positions are in the table for rs8176719, rs8176746, rs505922, rs657152, rs8176749, rs8176720, rs8176741, rs8176722 (both builds) and rs687289 (GRCh38 only); the other tag SNPs are found by ID only. O versus not O comes from
rs8176719 (the functional deletion; the shorter allele is O), else from the tag SNPs rs505922, rs657152, rs687289, rs612169,
rs529565 or rs514659; A versus B comes from rs8176746 (G = A, T = B), else from rs8176749, rs8176672, rs8176741, rs8176722 or
rs8176720. Alleles may be on either strand. A missing genotype at a SNP is replaced by the next SNP of the list. Types: A, B, AB, O;
`nonO` when the sample is not O and no A/B SNP is available; `NA` when no O genotype is available. A sample with one O allele that is
heterozygous at rs8176746 is called B (the O allele is assumed to be on an A background) and flagged `assumed`.
The tag SNPs are much less accurate than the functional ones, most of all in African ancestry. The Rh type is not inferred.
SNPs and rules: the ABO review at <https://www.medrxiv.org/content/10.1101/2024.02.21.24302372v1.full>.

## APOE genotype

`--inferAPOE` infers the APOE genotype (e2/e3/e4) of every sample from rs429358 and rs7412 (`--bfile`, `--vcf` or `--pfile`;
no model is needed) and writes `<prefix>_APOE.txt` (`IID`, `APOE`, `rs429358_C`, `rs7412_T`, `E4_copies`, `Flag`).
The SNPs are found by ID, else by chromosome 19 position (GRCh38 or GRCh37). Samples heterozygous at both SNPs are called e2/e4
and flagged `ambiguous` (e1/e3 is the other reading). If only one of the two SNPs is in the data, `--apoe-partial` reports e2 (from rs7412) or e4 (from rs429358) alone, flagged `partial`; without it a missing SNP is an error. `--inferABO` and `--inferAPOE` are run one at a time. The SNPs, positions
and rules of both options are in `ALGORITHMS.md`.

## Sex

`--infersex` infers the sex of every sample from the X chromosome (the F coefficient of PLINK `--check-sex`, outside the
pseudoautosomal regions) and the Y chromosome (the rate of genotypes), for `--bfile`, `--vcf` or `--pfile`, and writes
`<prefix>_sex.txt`. Where the data record a sex (.fam, .psam), `Status` says OK or PROBLEM. At most 50000 SNPs of the X chromosome and 50000 of the Y chromosome
are used. See `ALGORITHMS.md`.

```bash
./secretscores --infersex --bfile ex --prefix out
```

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
single-letter A/C/G/T alleles can be stored. A SNP without an other allele is stored with the
effect allele only (as SECRETS does; it cannot be matched to genotypes and is not scored). Binary
models of the earlier versions 0, 1 and 2 are no longer read (make them again with `--makebin`), and
SECRETS GWAS files are refused. If the text model has the column
`allelefrequency_effect`, the `.bin` file (version 4) stores the frequency of the effect allele of
each SNP as a 4-byte float, right after the weights; otherwise a version 3 file is written.
A text model with no valid SNP (for example, one whose `other_allele` column is empty throughout)
gets no `.bin` file; it is skipped with a note and listed in `<prefix>_scores.skipped`.
`--viewbin` prints a `.bin` model as text.

```bash
./secretscores --makebin --readmodel model.txt     # writes model.bin
./secretscores --bfile ex --readmodel model.bin
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

## Selecting samples

As in PLINK, `--keep` and `--remove` take a file with the family ID (`FID`) in the first column and
the individual ID (`IID`) in the second, separated by spaces or tabs. A header line starting with `#` is skipped, columns
after the second are ignored, and empty lines are skipped. The file may be gzipped. Samples are
matched to the `.fam` file by FID and IID together.

```bash
./secretscores --bfile ex --readmodel model.txt --keep keep.txt
./secretscores --bfile ex --readmodel model.txt --remove remove.txt
```

`--keep` keeps only the listed samples and `--remove` drops the listed ones; with both, `--keep` is
applied first and then `--remove`. The scores are written for the remaining samples, in `.fam` order.
The IDs in the files that are not in the `.fam` file are counted in a note, and it is an error if no
sample remains. The `.bed` file is still read as a whole, but only the genotypes of the remaining
samples are held in memory.

## Scoring many models

Give several files separated by commas:

```bash
./secretscores --bfile ex --readmodel a.txt,b.txt.gz,c.bin
```

or list them in a file, one model per line, with an optional second column naming the score
(lines starting with `#` are skipped):

```
# models.list
PGS000001.txt.gz   height
PGS000002.txt.gz   bmi
```

```bash
./secretscores --bfile ex --readlist models.list --prefix mine
```

Without a name, the score is named after the file (no directory, `.gz` or extension). Score
names must be unique.

## Options

| Option | Meaning |
|---|---|
| **Input Files** | |
| `--bfile <prefix>` | PLINK binary file prefix (`.bed/.bim/.fam`, `.bim`/`.fam` may be `.gz`); `ex` and `ex.bed` both work |
| `--vcf <file>` | VCF genotypes (hard calls), plain or `.gz`; instead of `--bfile` |
| `--pfile <prefix>` | PLINK 2 `.pgen/.pvar/.psam`; instead of `--bfile` |
| `--readmodel <files>` | Model file(s), comma-separated |
| `--readlist <file>` | File listing model files (and optional names) |
| `--readfreq <file>` | Allele frequencies of the genotype data (PLINK `.frq` or PLINK 2 `.afreq`), used for mean imputation instead of the model's `allelefrequency_effect` |
| **Sample Filter** | |
| `--keep <file>` | Keep only the samples listed in the file (FID and IID, as in PLINK) |
| `--remove <file>` | Remove the samples listed in the file (FID and IID, as in PLINK) |
| **Score Action** | |
| `--score` | Compute the scores (the default action) |
| `--makebin` | Convert text models to `.bin`; no genotypes are read |
| `--viewbin` | Print a `.bin` model as text |
| **Inferences** | |
| `--inferABO` | Infer the ABO blood type of every sample from the genotypes (no model needed); see below |
| `--inferAPOE` | Infer the APOE genotype (e2/e3/e4) of every sample from the genotypes (no model needed); see below |
| `--apoe-partial` | With `--inferAPOE`, report e2 or e4 alone when only one of the two SNPs is in the data |
| `--infersex` | Infer the sex of every sample from the X and Y chromosomes (no model needed); see Sex |
| **Score Parameter** | |
| `--no-mean-imputation` | Count a missing genotype as 0, do not use allele frequencies |
| `--allowX <0/1/2>` | X chromosome: 0 = not scored (default), 1 = males 0/1 and females 0/1/2, 2 = males and females 0/1/2 |
| **Other** | |
| `--cpus <n>` | Threads (default 0 = automatic, half the available cores) |
| `--prefix <name>` | Output file prefix (default `secretscores`) |

## Tips

- **Speed:** individual-major `.bed` files are faster than the default SNP-major ones. Convert once with
  `plink2 --bfile ex --export ind-major-bed --out ex_indmajor` and use `--bfile ex_indmajor`.
- **Memory:** about half of the available RAM is used to hold genotypes; many models are
  scored in batches of 100 to 1000 models (a batch stops once its models touch 9 million different SNPs).
- **Check coverage:** look at `*_scores.coverage` before trusting a score. Low coverage
  usually means a genome-build mismatch between the model and your data.
