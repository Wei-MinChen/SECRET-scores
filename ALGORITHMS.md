# Algorithms of `--inferABO`, `--inferAPOE` and `--infersex`

Both options read hard genotype calls from `--bfile` (PLINK bed), `--vcf` or `--pfile`, need no score model, and write one
line per sample to `<prefix>_ABO.txt` or `<prefix>_APOE.txt`. They are run one at a time. The code is in `src/InferABO.cpp`
and `src/InferAPOE.cpp` (the SNP tables are at the top of each file).

Common steps:

1. **Find each SNP of the table in the data.** First by variant ID (rsID); a SNP of the right ID but other alleles is not accepted.
   If it is not found by ID, by chromosome and position, in GRCh38 or GRCh37: every SNP not found by ID is tried in both builds
   and the build that finds more SNPs is used for all of them. A SNP found by position needs the same alleles as the SNP
   of the table, or their complement (the other strand). The run prints how each SNP was found.
2. **Count the copies of one allele** of each SNP found, per sample (0, 1, 2; missing if the genotype is missing).
   An allele given on the other strand is translated by complement. A/T and C/G SNPs are not used, since the strand cannot be told.
3. **Combine the counts by a fixed rule** to the type of the sample.

The positions and alleles in the tables are from dbSNP as shown by genebe.net, SNPedia and the medRxiv study below; they
were not checked against a genome browser. A position that is not in the tables is not searched (only the ID is).

## ABO blood type

Source of the SNPs: Disparities in ABO blood type determination across diverse ancestries (medRxiv 2024.02.21.24302372).

ABO is decided by the O versus not O question first, then by A versus B. Thirteen SNPs are checked, but only one SNP of
each question is used for a sample: the first of the list below that has a genotype for that sample, so a missing genotype at
the first SNP is replaced by the second, and so on. The functional SNPs are rs8176719 and rs8176746; the others are tags
(SNPs correlated with them), used only when the functional SNP is absent from the data or missing for a sample.

| Question | Order | SNP | Counted allele | Other allele | Positions in the table (GRCh38 / GRCh37) |
|---|---|---|---|---|---|
| O / not O | 1 | rs8176719 (deletion of a G, exon 6) | T = O (the shorter allele; TC = not O) | TC | 133257521 / 136132908 |
| | 2 | rs505922 | T = O | C | 133273813 / 136149229 |
| | 3 | rs657152 | C = O | A | 133263862 / 136139265 |
| | 4 | rs687289 | G = O | A | 133261703 / not in the table |
| | 5 | rs612169 | A = O | G | by ID only |
| | 6 | rs529565 | T = O | C | by ID only |
| | 7 | rs514659 | A = O | C | by ID only |
| A / B | 1 | rs8176746 (exon 7) | G = A | T = B | 133255935 / 136131322 |
| | 2 | rs8176749 | C = A | T = B | 133255801 / 136131188 |
| | 3 | rs8176672 | C = A | T = B | by ID only |
| | 4 | rs8176741 | G = A | A = B | 133256074 / 136131461 |
| | 5 | rs8176722 | C = A | A = B | 133257367 / 136132754 |
| | 6 | rs8176720 | T = A | C = B | 133257486 / 136132873 |

Not used: rs8176704 and rs8176693 (the most discordant with serology in the study) and rs8176645 (A/T, strand cannot be told).
The indel rs8176719 is accepted when one allele is a single letter and the other is longer (T against TC), or D against I.

The type of a sample, with `O_copies` (nO) the copies of the O allele at the SNP used for O, and `A_copies` (a) the copies of the
A-associated allele at the SNP used for A/B:

| nO | a | ABO | Reason |
|---|---|---|---|
| 2 | any | O | two O chromosomes; the A/B SNP is ignored |
| 0 | 2 | A | |
| 0 | 1 | AB | |
| 0 | 0 | B | |
| 1 | 2 | A | AO |
| 1 | 0 | B | BO |
| 1 | 1 | B, flag `assumed` | one chromosome is O, the other has T at the A/B SNP |
| 1 or 0 | no SNP or genotype | nonO | |
| no genotype | | NA | |

`A_copies` is a genotype at one SNP, not a count of A alleles. The O chromosomes carry the G of rs8176746, as the O allele descends
from A, so an O sample usually has `A_copies` = 2. In the row nO = 1 and a = 1, the O chromosome accounts for one G, so the
other chromosome has the T (B). The other reading, AO with the O allele on a B background, is rare, and is why the sample
is flagged.

Accuracy: rs8176719 and rs8176746 are the functional variants. The tag SNPs are much less accurate, above all in African
ancestry (median r2 with the deletion 0.44 in AFR, 0.87 in EUR, 0.95 in EAS). Only ABO is inferred, not the Rh type.
The program prints which SNP was used for how many samples.

## APOE genotype

Two SNPs, both needed (no tag SNPs):

| SNP | GRCh38 / GRCh37 (chr19) | Counted allele | Other |
|---|---|---|---|
| rs429358 | 44908684 / 45411941 | C | T |
| rs7412 | 44908822 / 45412079 | T | C |

The haplotypes, written as (rs429358, rs7412): e2 = T,T; e3 = T,C; e4 = C,C; e1 = C,T (very rare, not a usual allele).
With c the copies of C at rs429358 and t the copies of T at rs7412, the data are unphased, so the two haplotypes
are the ones that fit the counts:

| c \ t | 0 | 1 | 2 |
|---|---|---|---|
| 0 | e3/e3 | e2/e3 | e2/e2 |
| 1 | e3/e4 | e2/e4 (flag `ambiguous`) | e1/e2 (flag `rare`) |
| 2 | e4/e4 | e1/e4 (flag `rare`) | e1/e1 (flag `rare`) |

The only real ambiguity is a sample heterozygous at both SNPs: e2/e4 or e1/e3. It is called e2/e4 because e1 is extremely rare;
the sample is marked `ambiguous` (phased data would settle it). A genotype with e1 is marked `rare`; check the data
for such samples. `E4_copies` is the number of e4 alleles. A sample with a missing genotype at either SNP is NA.

The two SNPs are typed on most arrays, are in sequencing data, and are usually imputed well; rs429358 is the one that is
harder to type, so check the missing rate of the two SNPs.

### `--apoe-partial`

By default a missing SNP is an error. With `--apoe-partial`, a data set with only one of the two SNPs gives what that SNP alone tells,
with `Flag` = `partial` for every sample:

| SNP in the data | `APOE` | Meaning |
|---|---|---|
| rs7412 only | `no_e2`, `e2_het`, `e2_hom` | T copies 0, 1, 2: no e2, one e2, two e2 (e1, extremely rare, also carries T); e3 and e4 are not told apart (`E4_copies` is NA) |
| rs429358 only | `no_e4`, `e4_het`, `e4_hom` | C copies 0, 1, 2: no e4, one e4, two e4 (e1, extremely rare, also carries C); e2 and e3 are not told apart (`E4_copies` = C copies) |

If neither SNP is in the data it is an error. No proxy SNP is used for a missing SNP.

## Sex (`--infersex`)

Unlike the two options above, this one uses many SNPs: the X and Y chromosomes of the data (chromosome `X`/`23` and `Y`/`24`,
any build). It writes `<prefix>_sex.txt` with `IID`, `RecordedSex` (the .fam/.psam sex; NA when not recorded, as in a VCF),
`InferredSex` (M, F or NA), `Status` (OK when the inferred sex is the recorded one, PROBLEM when not, NA when one of the two
is not known), `F_X`, `X_het`, `X_snps`, `Y_rate`, `Y_snps` and `Flag`.

SNPs: the pseudoautosomal regions (PAR), where the two sexes both have two copies, are left out. They are the regions of
GRCh38 (X 10,001-2,781,479 and 155,701,383-156,030,895; Y 10,001-2,781,479 and 56,887,903-57,217,415) and of GRCh37 (X 60,001-2,699,520
and 154,931,044-155,260,560; Y 10,001-2,649,520 and 59,034,050-59,363,566); the regions of the two builds are all left out, so the
build need not be known (the cost is a few SNPs of the non-PAR region of the other build). Chromosome XY / 25 is not used. At least 20
SNPs of X or of Y are needed. When there are more than 50000, that number of SNPs, evenly spaced, are used.

**X: the F coefficient** (as in PLINK `--check-sex`). Only the SNPs genotyped for at least 80% of the samples and of allele
frequency 0.05 to 0.95 are used (the frequency is that of all the samples). For a sample, with `O` its heterozygous genotypes
and `E` the expected number of them, the sum over its genotyped SNPs of `2p(1-p) * 2n/(2n-1)` (n the samples genotyped at the SNP):

    F = 1 - O / E

A male has one X chromosome (a haploid call is read as homozygous, a heterozygous call is an error), so `O` is near 0 and F near 1; a
female has F near 0. F > 0.8 is male, F < 0.2 is female, and between them ambiguous (`X_ambiguous`). A sample needs 20 usable X SNPs.

**Y: the genotype rate.** The proportion of the Y SNPs (outside PAR) with a genotype. On genotyping arrays a female has few
genotypes at the Y probes, and a male most. The samples are sorted by this proportion and split at the largest gap between two
neighbours, if the gap is at least 0.2; the group with the higher proportions is male. When there is no such gap (all the samples of
one sex, or no difference), Y is not used. Sequencing data where every sample has a call at every Y site (a female with
reference calls) cannot be told by this rule.

**The call.**

| X | Y | InferredSex | Flag |
|---|---|---|---|
| male or female | same, or not used | that sex | NA |
| male or female | the other sex | NA | `X_Y_conflict` (an aneuploidy such as XXY or XO, a sample swap or mixture, or a poor sample) |
| ambiguous or not available | male or female | the Y call | `Y_only` or `X_ambiguous_Y_used` |
| ambiguous | not used | NA | `X_ambiguous` |
| not available | not used | NA | |

`Status` is PROBLEM when the sex of the .fam/.psam file is not the inferred one (samples of unknown sex, 0, are NA).

Cautions: F depends on the allele frequencies of the samples and is less reliable for a few samples (below about 20),
for samples of mixed ancestries or many close relatives, and for data with few X SNPs; the program prints the number of SNPs
used. The rule recognizes the sex chromosome content, not the gender or the sex of a person with a chromosome disorder.
