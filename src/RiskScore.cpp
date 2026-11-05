////////////////////////////////////////////////////////////////////////
// RiskScore.cpp
// (c) 2026-2026 Wei-Min Chen
//
// This file is distributed as part of the secretscores source code
// package and may not be redistributed in any form, without prior
// written permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile secretscores.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Polygenic score calculation.
//
// Text model file: a PGS Catalog scoring file (tab-delimited, header lines start
// with #).  Only these columns are used:
//     effect_allele  effect_weight  hm_chr  hm_pos
//     other_allele  and/or  hm_inferOtherAllele
// (other_allele is used when it has a value, else hm_inferOtherAllele).
// The header lines (starting with #) of a PGS Catalog file are read for the trait annotation:
//     #pgs_id=  #pgs_name=  #trait_reported=  #trait_mapped=  #trait_efo=
// which go to the "trait" and "trait_efo" columns of the coverage file.
// Binary model file (name ends in .bin, made with --makebin): see BinaryInput.cpp.
//
// A model line needs an other allele, and it has to differ from the effect
// allele; lines without one are dropped (hm_inferOtherAllele often supplies it).
//
// A model SNP matches a SNP of the genotype data when chromosome and base-pair
// position are the same and the alleles are the same: effect_allele = A1 with
// other_allele = A2, or effect_allele = A2 with other_allele = A1 (the .bim
// allele columns).  Alleles are compared as they are; there is no strand
// flipping and no allele-frequency based decision.
//
// effect_weight always stays the weight of the effect allele (the PGS Catalog
// definition): score = sum of effect_weight x (number of effect alleles), counting
// allele A1 or allele A2 of the data, whichever is the effect allele.  No allele
// frequency decides which SNP matches.  Missing genotypes: by default (mean imputation) a
// missing genotype counts as 2 x the frequency of the effect allele, from --readfreq (a SNP
// that is not in the file gets no imputation, with a warning) or, without --readfreq, from the
// allelefrequency_effect column of the model; with no frequency, or with --no-mean-imputation,
// it contributes 0.
//
// The X chromosome is scored (see XChrModel in ScoreEngine.h): the sex of each person comes from
// the .fam file.  The genotype codes are used as they are stored in the .bed file, in which a
// male's hemizygous call outside the PAR is a homozygous genotype, so it counts as two copies of
// the effect allele (--allowX 2, as in PLINK 2) or as one (--allowX 1); a
// heterozygous genotype of a male outside the PAR is impossible and counts as missing (0).  The PAR
// has to be coded XY (25) in the .bim file, as PLINK requires (--split-x), and is scored like an
// autosome.  Females and people of unknown sex are scored as diploid.  The X SNPs of a model are
// added to its score after its other SNPs (each group in the order of the model file).  The Y
// chromosome and the mitochondrial DNA are not scored.
//
// Several models can be scored in one run.  They are processed in batches whose
// genotypes fit in GENOTYPE_MEMORY_LIMIT (half of the RAM available when the
// genotypes are about to be loaded), so the .bed file is read once for all
// the SNPs of a batch.  The genotypes are stored individual-major (see ReadPLINK.cpp),
// and the people are scored in parallel.
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <string>
#include <vector>
#include <utility>
#include <unordered_set>
#include <atomic>
#include <time.h>
#include <stdarg.h>
#include "ScoreEngine.h"
#include "Error.h"
#include "MemoryInfo.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

// genotypes held in memory at once (bytes); a single model larger than this is
// still loaded whole.  Set before the genotypes are loaded: half of the available RAM.
static long long GENOTYPE_MEMORY_LIMIT = 2LL * 1024 * 1024 * 1024;   // used if the RAM cannot be detected

// the models with no more SNPs than this get a table of their SNPs on screen
#define SNP_TABLE_LIMIT 100

// --allowX (see ScoreEngine.h)
int XChrModel = 0;
bool MeanImputation = true;

// the chromosomes that are not scored, for messages
static const char *NotScoredChromosomes()
{
   return XChrModel == 0? "X, Y, MT": "Y, MT";
}

// A score name has no white space, comma or double quote (the score file is tab-delimited, the
// coverage file comma-delimited, and SECRET uses the names as column names); empty is SCORE.
static std::string SafeScoreName(const std::string & name)
{
   std::string nm = name;
   for(size_t i = 0; i < nm.size(); i++)
      if(nm[i] == ' ' || nm[i] == '\t' || nm[i] == '\r' || nm[i] == '\n' || nm[i] == ',' || nm[i] == '"') nm[i] = '_';
   if(nm.empty()) nm = "SCORE";
   return nm;
}

// a finite number, from the whole token
static bool ParseNumber(const std::string & token, double & value)
{
   const char *s = token.c_str();
   char *end = NULL;
   value = strtod(s, &end);
   return (end != s) && (*end == '\0') && (value == value) && (fabs(value) < 1e300);
}

// a positive integer, from the whole token
static bool ParsePosition(const std::string & token, long & value)
{
   const char *s = token.c_str();
   char *end = NULL;
   value = strtol(s, &end, 10);
   return (end != s) && (*end == '\0') && value > 0 && value < 2000000000L;
}

static bool IsMissingValue(const std::string & s)
{
   if(s.empty() || s == "." ) return true;
   if(s.size() == 2 && (s[0]=='N' || s[0]=='n') && (s[1]=='A' || s[1]=='a')) return true;
   return false;
}

static std::string Upper(const std::string & s)
{
   std::string u = s;
   for(size_t i = 0; i < u.size(); i++)
      if(u[i] >= 'a' && u[i] <= 'z') u[i] = (char)(u[i] - 'a' + 'A');
   return u;
}

// Split a line into fields.  A tab-delimited line keeps its empty fields (PGS Catalog files
// leave cells empty), so that a column is always the same field; otherwise white space
// separates the fields.
static void SplitFields(const char *s, bool tabs, std::vector<std::string> & fields)
{
   fields.clear();
   const char *p = s;
   if(tabs){
      for(;;){
         const char *q = p;
         while(*q && *q != '\t' && *q != '\r' && *q != '\n') q++;
         fields.push_back(std::string(p, q - p));
         if(*q != '\t') break;
         p = q + 1;
      }
   }else{
      while(*p){
         while(*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
         const char *q = p;
         while(*q && *q != ' ' && *q != '\t' && *q != '\r' && *q != '\n') q++;
         if(q > p) fields.push_back(std::string(p, q - p));
         p = q;
      }
   }
}

static bool SameName(const std::string & a, const char *b)
{
   return Upper(a) == Upper(std::string(b));
}

static std::string Trim(const std::string & s)
{
   size_t b = 0, e = s.size();
   while(b < e && (s[b] == ' ' || s[b] == '\t' || s[b] == '\r' || s[b] == '\n')) b++;
   while(e > b && (s[e-1] == ' ' || s[e-1] == '\t' || s[e-1] == '\r' || s[e-1] == '\n')) e--;
   return s.substr(b, e - b);
}

// A header line of a PGS Catalog scoring file, "#key=value" (e.g. "#trait_reported=Breast
// cancer"): keep the annotation that is used.  Other comment lines ("###PGS CATALOG SCORING
// FILE ...", "# Row Headers", free text) have no "key=value" with a one-word key and are ignored.
static void ParseHeaderLine(const char *s, TextModelStats & stats)
{
   while(*s == '#') s++;
   const char *eq = strchr(s, '=');
   if(eq == NULL) return;
   const std::string key = Trim(std::string(s, eq - s));
   if(key.empty() || key.find_first_of(" \t") != std::string::npos) return;
   const std::string value = Trim(std::string(eq + 1));
   if(value.empty()) return;
   if(SameName(key, "pgs_id")) stats.pgsId = value;
   else if(SameName(key, "pgs_name")) stats.pgsName = value;
   else if(SameName(key, "trait_reported")) stats.traitReported = value;
   else if(SameName(key, "trait_mapped")) stats.traitMapped = value;
   else if(SameName(key, "trait_efo")) stats.traitEfo = value;
}

std::string ModelTrait(const TextModelStats & stats)
{
   if(!stats.traitReported.empty()) return stats.traitReported;
   if(!stats.traitMapped.empty()) return stats.traitMapped;
   return stats.traitEfo;
}

// one field of a comma-delimited file: in double quotes when it has a comma, quote or line break
static std::string CsvField(const std::string & value)
{
   if(value.find_first_of(",\"\r\n") == std::string::npos) return value;
   std::string q = "\"";
   for(size_t i = 0; i < value.size(); i++){
      if(value[i] == '"') q += '"';
      if(value[i] == '\r' || value[i] == '\n') q += ' ';
      else q += value[i];
   }
   return q + "\"";
}

static std::string CsvCount(long long n)
{
   if(n < 0) return "NA";
   char buf[32];
   snprintf(buf, sizeof(buf), "%lld", n);
   return buf;
}

void WriteCoverageFile(const std::string & filename, const std::vector<CoverageRow> & rows, bool withGenotypes)
{
   FILE *fp = fopen(filename.c_str(), "wt");
   if(fp == NULL) error("Cannot open %s to write", filename.c_str());
   // the trait and trait_efo columns are last
   const bool withImputed = withGenotypes && MeanImputation;
   fprintf(fp, "SCORE,SNPs,%sNoCoord,NoOtherAllele,SameAllele,NonAutosomal%s%s,trait,trait_efo\n",
      withGenotypes? "Found,Used,coverage,": "", withGenotypes? ",NotInData,AlleleMismatch": "",
      withImputed? ",PercentImputed": "");
   for(size_t k = 0; k < rows.size(); k++){
      const CoverageRow & r = rows[k];
      fprintf(fp, "%s,%s,", CsvField(r.name).c_str(), CsvCount(r.fileSNPs).c_str());
      if(withGenotypes){
         char coverage[32];
         if(r.used >= 0 && r.fileSNPs > 0) snprintf(coverage, sizeof(coverage), "%.4f", (double) r.used / r.fileSNPs);
         else snprintf(coverage, sizeof(coverage), "NA");
         fprintf(fp, "%s,%s,%s,", CsvCount(r.found).c_str(), CsvCount(r.used).c_str(), coverage);
      }
      fprintf(fp, "%s,%s,%s,%s", CsvCount(r.noCoord).c_str(), CsvCount(r.noOtherAllele).c_str(),
         CsvCount(r.sameAllele).c_str(), CsvCount(r.nonAutosomal).c_str());
      if(withGenotypes)
         fprintf(fp, ",%s,%s", CsvCount(r.notInData).c_str(), CsvCount(r.alleleMismatch).c_str());
      if(withImputed){
         if(r.percentImputed >= 0) fprintf(fp, ",%.4f", r.percentImputed);
         else fprintf(fp, ",NA");
      }
      fprintf(fp, ",%s,%s\n", CsvField(r.trait.empty()? "NA": r.trait).c_str(),
         CsvField(r.traitEfo.empty()? "NA": r.traitEfo).c_str());
   }
   fclose(fp);
}

// Parse a PGS Catalog text model without genotype data (see ScoreEngine.h).
void ParseTextModel(const std::string & file, TextModelStats & stats, const SNPCallback & onSNP)
{
   LineReader reader;
   FILE *test = fopen(file.c_str(), "rb");
   if(test == NULL) ThrowModelError("Cannot open %s to read", file.c_str());
   fclose(test);
   if(!reader.Open(file.c_str())) ThrowModelError("Cannot open %s to read", file.c_str());
   stats.fileSNPs = stats.malformed = stats.noCoord = stats.noOtherAllele = stats.sameAllele = 0;
   stats.hasFreq = false;
   stats.pgsId = stats.pgsName = stats.traitReported = stats.traitMapped = stats.traitEfo = "";
   std::vector<std::string> f;
   int colEA = -1, colWT = -1, colCHR = -1, colPOS = -1, colOA = -1, colHMOA = -1, colAF = -1;
   bool haveHeader = false, tabs = true;
   size_t needCols = 0;
   const char *s;
   size_t length;
   while(reader.Next(s, length)){
      if(s[0] == '\0') continue;
      if(s[0] == '#'){
         if(!haveHeader) ParseHeaderLine(s, stats);      // the trait annotation of a PGS Catalog file
         continue;
      }
      if(!haveHeader){
         haveHeader = true;
         tabs = (strchr(s, '\t') != NULL);
         SplitFields(s, tabs, f);
         for(size_t i = 0; i < f.size(); i++){
            if(SameName(f[i], "effect_allele")) colEA = (int) i;
            else if(SameName(f[i], "effect_weight")) colWT = (int) i;
            else if(SameName(f[i], "hm_chr")) colCHR = (int) i;
            else if(SameName(f[i], "hm_pos")) colPOS = (int) i;
            else if(SameName(f[i], "other_allele")) colOA = (int) i;
            else if(SameName(f[i], "hm_inferOtherAllele")) colHMOA = (int) i;
            else if(SameName(f[i], "allelefrequency_effect")) colAF = (int) i;
         }
         if(colEA == -1 || colWT == -1 || colCHR == -1 || colPOS == -1 || (colOA == -1 && colHMOA == -1)){
            std::string missing, found;
            if(colEA == -1) missing += " effect_allele";
            if(colWT == -1) missing += " effect_weight";
            if(colCHR == -1) missing += " hm_chr";
            if(colPOS == -1) missing += " hm_pos";
            if(colOA == -1 && colHMOA == -1) missing += " (other_allele or hm_inferOtherAllele)";
            for(size_t i = 0; i < f.size() && i < 40; i++) found += (i? ", ": "") + f[i];
            if(f.size() > 40) found += ", ...";
            ThrowModelError("The header of the model file %s needs the columns\n"
               "  effect_allele, effect_weight, hm_chr, hm_pos\n"
               "  and other_allele or hm_inferOtherAllele (PGS Catalog harmonized scoring file)\n"
               "  Missing:%s\n  Columns found (%d): %s", file.c_str(), missing.c_str(), (int) f.size(), found.c_str());
         }
         needCols = (size_t) std::max(std::max(colEA, colWT), std::max(colCHR, colPOS)) + 1;
         if(colOA + 1 > (int) needCols) needCols = (size_t) colOA + 1;
         if(colHMOA + 1 > (int) needCols) needCols = (size_t) colHMOA + 1;
         if(colAF + 1 > (int) needCols) needCols = (size_t) colAF + 1;
         stats.hasFreq = (colAF != -1);
         continue;
      }
      SplitFields(s, tabs, f);
      // a trailing empty cell can be cut off the end of the line
      if(f.size() < needCols) f.resize(needCols);
      double wt;
      std::string ea = Upper(f[colEA]);
      if(IsMissingValue(ea) || !ParseNumber(f[colWT], wt)){stats.malformed++; continue;}
      stats.fileSNPs++;
      long pos;
      int chr = IsMissingValue(f[colCHR])? 0: ChromosomeCode(f[colCHR].c_str(), f[colCHR].size());
      if(chr == 0 || !ParsePosition(f[colPOS], pos)){stats.noCoord++; continue;}
      // the other allele: other_allele, else hm_inferOtherAllele; a line without one is dropped
      std::string oa;
      if(colOA != -1 && !IsMissingValue(f[colOA])) oa = Upper(f[colOA]);
      else if(colHMOA != -1 && !IsMissingValue(f[colHMOA])) oa = Upper(f[colHMOA]);
      if(oa.empty()){stats.noOtherAllele++; continue;}
      if(oa == ea){stats.sameAllele++; continue;}
      // allelefrequency_effect: the frequency of the effect allele; NaN when missing or not in [0, 1]
      double af = std::numeric_limits<double>::quiet_NaN();
      if(colAF != -1 && !IsMissingValue(f[colAF])){
         double v;
         if(ParseNumber(f[colAF], v) && v >= 0.0 && v <= 1.0) af = v;
      }
      onSNP(chr, pos, ea, oa, wt, af);
   }
   {
      std::string why;
      if(reader.Failed(why))
         ThrowModelError("Cannot read %s completely (%s): the file may be truncated or damaged", file.c_str(), why.c_str());
   }
   if(!haveHeader) ThrowModelError("The model file %s has no header line", file.c_str());
}

// Read one text model file and match its SNPs to the genotype data (chromosome, position, alleles).
void ScoreEngine::ReadModelFile(ScoreModel & model)
{
   model.rowIdx.clear();
   model.weights.clear();
   model.switched.clear();
   model.freq.clear();
   model.nonAutosomal = model.notInData = model.alleleMismatch = 0;
   ScoreModel *md = &model;
   TextModelStats stats;
   ParseTextModel(model.file, stats,
      [this, md](int chr, long pos, const std::string & ea, const std::string & oa, double wt, double af){
         if(!ScoredChromosome(chr)){md->nonAutosomal++; return;}
         const unsigned int eaCode = AlleleCode(ea.c_str(), ea.size());
         const unsigned int oaCode = AlleleCode(oa.c_str(), oa.size());
         bool sw = false, atPosition = false;
         const int hit = MatchSNP(chr, pos, eaCode, oaCode, sw, atPosition);
         if(hit == -1){
            if(atPosition) md->alleleMismatch++;
            else md->notInData++;
            return;
         }
         md->rowIdx.push_back(hit);
         md->weights.push_back(wt);
         md->switched.push_back(sw? 1: 0);
         md->freq.push_back(MeanImputation? af: std::numeric_limits<double>::quiet_NaN());
      });
   model.fileSNPs = stats.fileSNPs;
   model.malformed = stats.malformed;
   model.noCoord = stats.noCoord;
   model.noOtherAllele = stats.noOtherAllele;
   model.sameAllele = stats.sameAllele;
   model.trait = ModelTrait(stats);
   model.traitEfo = stats.traitEfo;
}

// Tabulate for each matched SNP what every genotype code contributes to the score: the effect
// weight w times the copies of the effect allele (2, 1 or 0), and 0 for a missing genotype.
// A missing genotype contributes 2 x the effect allele frequency x w with mean imputation, else 0.
// The effect allele is A1 (copies of A1 are 2, 1, 0 for the codes 0, 2, 3) or A2 (copies of A2
// are 0, 1, 2 for the codes 0, 2, 3; SWITCHED).  The weight is never changed.
void ScoreEngine::PrepareModel(ScoreModel & model, bool verbose)
{
   const int snpcount = (int) model.rowIdx.size();
   model.vslot.clear();
   model.vval.clear();
   model.xslot.clear();
   model.xval.clear();
   model.xvalM.clear();
   model.vslot.reserve(snpcount);
   model.vval.reserve(4 * (size_t) snpcount);
   model.freqFromFile = model.freqFromModel = model.freqNone = 0;
   model.impMode = 0;
   model.imputedGenotypes = 0;
   model.vhas.clear();
   model.xhas.clear();
   int freqSkipped = 0;          // SNPs that --readfreq has no frequency for: no imputation
   const double NaN = std::numeric_limits<double>::quiet_NaN();
   const bool table = verbose && snpcount <= SNP_TABLE_LIMIT;
   if(table){
      printf("\nThe %d SNPs of the model that are in the genotype data\n", snpcount);
      printf("%5s %11s %5s %5s %5s %12s %9s %8s\n", "CHR", "POS", "A1", "A2", "EA", "WEIGHT", "ACTION", "EAF");
   }
   for(int m = 0; m < snpcount; m++){
      const int row = model.rowIdx[m];
      const bool sw = model.switched[m] != 0;
      const double w = model.weights[m];
      // the frequency of the effect allele for mean imputation of a missing genotype: from --readfreq
      // (given for allele A1 of the data) when it has the SNP, else from the model; NaN when none
      double eaf = NaN;
      if(MeanImputation){
         if(haveReadFreq){
            // with --readfreq every SNP is expected to be in the file: for a SNP that is not, the
            // imputation is skipped (a missing genotype counts as 0; the SNP is still scored)
            const double b = bimFreqA1[row];
            if(b >= 0.0 && b <= 1.0){eaf = sw? 1.0 - b: b; model.freqFromFile++;}
            else{freqSkipped++; model.freqNone++;}
         }
         else{
            const double mf = model.freq.empty()? NaN: model.freq[m];
            if(mf >= 0.0 && mf <= 1.0){eaf = mf; model.freqFromModel++;}
            else{eaf = NaN; model.freqNone++;}
         }
      }
      const bool haveEaf = (eaf >= 0.0 && eaf <= 1.0);
      if(table){
         char a1[8], a2[8];
         const unsigned int c1 = bimA1[row], c2 = bimA2[row];
         snprintf(a1, sizeof(a1), "%s", c1 < 256? std::string(1, (char) c1).c_str(): "<m>");
         snprintf(a2, sizeof(a2), "%s", c2 < 256? std::string(1, (char) c2).c_str(): "<m>");
         char fq[16];
         if(haveEaf) snprintf(fq, sizeof(fq), "%.4f", eaf); else snprintf(fq, sizeof(fq), "NA");
         printf("%5d %11d %5s %5s %5s %12.6g %9s %8s\n", (int) bimChr[row], bimPos[row], a1, a2,
            sw? a2: a1, w, sw? "SWITCHED": "OK", fq);
      }
      // a missing genotype counts as 2 x frequency x weight (0 without a frequency)
      const double miss = haveEaf? 2 * eaf * w: 0.0;
      if(bimChr[row] == SEXCHR){
         model.xhas.push_back(haveEaf? 1: 0);
         // X chromosome: females and people of unknown sex are diploid, as the autosomes
         model.xslot.push_back(slotOfRow[row]);
         model.xval.push_back(sw? 0.0: 2 * w);
         model.xval.push_back(miss);
         model.xval.push_back(w);
         model.xval.push_back(sw? 2 * w: 0.0);
         // males: code 0 and 3 are the two homozygous calls, which stand for the one allele of
         // the hemizygous male; it counts as 2 copies (--allowX 2) or 1 (--allowX 1);
         // code 2 (heterozygous) is impossible outside the PAR and is taken as missing
         const double one = (XChrModel == 1)? w: 2 * w;
         const double missM = haveEaf? ((XChrModel == 1)? eaf: 2 * eaf) * w: 0.0;
         model.xvalM.push_back(sw? 0.0: one);
         model.xvalM.push_back(missM);
         model.xvalM.push_back(missM);          // heterozygous: impossible, taken as missing
         model.xvalM.push_back(sw? one: 0.0);
         continue;
      }
      model.vslot.push_back(slotOfRow[row]);
      model.vhas.push_back(haveEaf? 1: 0);
      model.vval.push_back(sw? 0.0: 2 * w);      // code 0: two copies of A1
      model.vval.push_back(miss);                // code 1: missing
      model.vval.push_back(w);                   // code 2: one copy of A1
      model.vval.push_back(sw? 2 * w: 0.0);      // code 3: no copy of A1
   }
   if(freqSkipped > 0){
      printf("WARNING - %s: the --readfreq file has no frequency for %d SNP%s; no mean imputation for %s (still scored)\n",
         model.name.c_str(), freqSkipped, freqSkipped == 1? "": "s", freqSkipped == 1? "it": "them");
      fflush(stdout);
   }
   model.freqSkipped = freqSkipped;
   {
      const int withFreq = model.freqFromFile + model.freqFromModel;
      const int all = (int) (model.vslot.size() + model.xslot.size());
      model.impMode = (withFreq == 0)? 0: (withFreq == all? 1: 2);
      if(model.impMode != 2){
         std::vector<char>().swap(model.vhas);
         std::vector<char>().swap(model.xhas);
      }
   }
   model.used = (int) (model.vslot.size() + model.xslot.size());
   model.scored = (model.used > 0);
}

// Score the prepared models first..last-1.  People are independent, so they are split among
// the threads; for each person all the models are scored while the person's genotypes are in
// cache, and each person adds up the SNPs of a model in the order of the model file.
void ScoreEngine::ScoreBatch(ScoreModel *models, int first, int last,
   std::vector<std::vector<double> > & scores, std::vector<double> * missSNP)
{
   // One parallel loop over all the people (dynamic schedule in small chunks, so that every thread
   // stays busy whatever the number of people).  Progress is reported by thread 0 only, every 30
   // seconds, from a shared counter of the people done, so tracking costs next to nothing.
   const time_t start = time(NULL);
   time_t lastReport = start;
   long long hetMale = 0;        // heterozygous genotypes of males at X SNPs (taken as missing)
   std::atomic<int> peopleDone(0);
   int chunk = idCount / (8 * (threadCount > 0? threadCount: 1));
   if(chunk < 1) chunk = 1;
   if(chunk > 64) chunk = 64;
#ifdef _OPENMP
   #pragma omp parallel num_threads(threadCount) if(threadCount > 1) reduction(+:hetMale)
#endif
   {
#ifdef _OPENMP
      const bool reporter = (omp_get_thread_num() == 0);
#else
      const bool reporter = true;
#endif
      std::vector<long long> impLocal(last - first, 0);   // imputed genotypes of this thread, per model
#ifdef _OPENMP
      #pragma omp for schedule(dynamic, chunk) nowait
#endif
      for(int i = 0; i < idCount; i++){
         const unsigned char *gi = G_IND[i];
         const bool male = (sex[i] == 1);
         for(int k = first; k < last; k++){
            ScoreModel & md = models[k];
            if(!md.scored) continue;
            const int V = (int) md.vslot.size();
            const int *slot = md.vslot.data();
            const double *val = md.vval.data();
            double sum = 0.0;
            int ms = 0;
            for(int v = 0; v < V; v++){
               const int s = slot[v];
               const int g = (gi[s >> 2] >> ((s & 3) * 2)) & 3;
               sum += val[4*v + g];
               ms += (g == 1);
            }
            // the X chromosome: a male's table takes a heterozygous genotype as missing
            const int XV = (int) md.xslot.size();
            if(XV > 0){
               const int *xs = md.xslot.data();
               const double *xv = male? md.xvalM.data(): md.xval.data();
               for(int v = 0; v < XV; v++){
                  const int s = xs[v];
                  const int g = (gi[s >> 2] >> ((s & 3) * 2)) & 3;
                  sum += xv[4*v + g];
                  const int het = male && (g == 2);
                  ms += (g == 1) || het;
                  hetMale += het;
               }
            }
            // imputed genotypes: all the missing ones when every SNP has a frequency; with some SNPs
            // having none, the missing ones at the SNPs with one are counted in a second pass
            if(md.impMode == 1) impLocal[k - first] += ms;
            else if(md.impMode == 2){
               int imp = 0;
               for(int v = 0; v < V; v++){
                  const int s = slot[v];
                  imp += (((gi[s >> 2] >> ((s & 3) * 2)) & 3) == 1) & md.vhas[v];
               }
               for(int v = 0; v < XV; v++){
                  const int s = md.xslot[v];
                  const int g = (gi[s >> 2] >> ((s & 3) * 2)) & 3;
                  imp += ((g == 1) || (male && g == 2)) & md.xhas[v];
               }
               impLocal[k - first] += imp;
            }
            scores[k][i] = sum;
            if(missSNP && k == first)
               (*missSNP)[i] = (double) ms / md.used;
         }
         const int done = peopleDone.fetch_add(1, std::memory_order_relaxed) + 1;
         if(reporter){
            const time_t now = time(NULL);
            if(now - lastReport >= 30 && done < idCount){
               const double elapsed = difftime(now, start);
               const double frac = (double) done / idCount;
               printf("  Scoring: %d of %d individuals (%.0f%%) done, %.0f s elapsed, about %.0f s left\n",
                  done, idCount, 100.0 * frac, elapsed, elapsed * (1 - frac) / frac);
               fflush(stdout);
               lastReport = now;
            }
         }
      }
#ifdef _OPENMP
      #pragma omp critical
#endif
      for(int k = first; k < last; k++) models[k].imputedGenotypes += impLocal[k - first];
   }
   maleHetX += hetMale;
}

void ThrowModelError(const char *format, ...)
{
   char message[4096];
   va_list ap;
   va_start(ap, format);
   vsnprintf(message, sizeof(message), format, ap);
   va_end(ap);
   throw ModelError(message);
}

std::string OneLine(const std::string & message)
{
   std::string out;
   bool blank = false;
   for(size_t i = 0; i < message.size(); i++){
      const char c = message[i];
      if(c == ' ' || c == '\t' || c == '\n' || c == '\r'){blank = true; continue;}
      if(blank && !out.empty()) out += ' ';
      blank = false;
      out += c;
   }
   return out;
}

void WriteSkippedFile(const std::string & filename, const std::vector<std::string> & files,
   const std::vector<std::string> & reasons)
{
   FILE *fp = fopen(filename.c_str(), "wt");
   if(fp == NULL){
      printf("WARNING - cannot open %s to write the list of skipped models\n", filename.c_str());
      return;
   }
   fprintf(fp, "FILE\tREASON\n");
   for(size_t i = 0; i < files.size(); i++)
      fprintf(fp, "%s\t%s\n", files[i].c_str(), OneLine(reasons[i]).c_str());
   fclose(fp);
}

int ResolveThreads(int cpus)
{
   int threads = 1;
#ifdef _OPENMP
   if(cpus > 0)
      threads = cpus;
   else{
      int availableCoreCount = omp_get_max_threads();
      threads = availableCoreCount == 1? 1: availableCoreCount / 2;
   }
#else
   (void) cpus;
#endif
   return threads;
}

// score name from a model file name: no directory, no .gz, no final extension
std::string ModelNameFromFile(const std::string & file)
{
   std::string name = file;
   size_t slash = name.find_last_of("/\\");
   if(slash != std::string::npos) name = name.substr(slash + 1);
   if(name.size() > 3 && name.substr(name.size() - 3) == ".gz") name = name.substr(0, name.size() - 3);
   size_t dot = name.find_last_of('.');
   if(dot != std::string::npos && dot > 0) name = name.substr(0, dot);
   return name;
}

void ScoreEngine::PrintGeneticRiskScore(const std::string & binfile, std::vector<std::string> & modelFiles,
   std::vector<std::string> & modelNames)
{
   int modelCount = (int) modelFiles.size();      // models that cannot be read are removed below
   if(modelCount == 0) error("Please use --readmodel <file> or --readlist <file> to specify score models.");

   // ---- the files ----
   if(binfile.find(',') != std::string::npos)
      error("Multiple genotype files cannot be read currently.");
   plinkinput = binfile;
   if(binfile.size() > 4 && binfile.substr(binfile.size() - 4) == ".bed")
      plinkinput = binfile.substr(0, binfile.size() - 4);
   printf("Loading genotype data in PLINK binary format...\n");
   ReadPlinkFam(plinkinput);
   ReadPlinkBim(plinkinput);
   if(bimCount == 0)
      error("No SNPs were found in %s.bim", plinkinput.c_str());
   if(!readFreqFile.empty()){
      if(MeanImputation) ReadFreqFile(readFreqFile);
      else printf("--readfreq %s is not read: the frequencies are only used for mean imputation (off with --no-mean-imputation)\n",
         readFreqFile.c_str());
   }

   // ---- names of the scores ----
   // names have no white space (the score file is tab-delimited, and SECRET uses them as column names)
   ScoreModel *models = new ScoreModel[modelCount];
   std::vector<std::string> usedNames;
   for(int k = 0; k < modelCount; k++){
      models[k].file = modelFiles[k];
      models[k].binary = IsBinaryModelName(modelFiles[k]);
      const std::string nm = SafeScoreName(modelNames[k].empty()? ModelNameFromFile(modelFiles[k]): modelNames[k]);
      for(int j = 0; j < k; j++)
         if(models[j].name == nm)
            error("Scores %s and %s have the same name %s.\n  Give them different names in the second column of the --readlist file.",
               models[j].file.c_str(), modelFiles[k].c_str(), nm.c_str());
      models[k].name = nm;
      models[k].scored = false;
      models[k].used = 0;
      models[k].freqSkipped = 0;
   }

   // ---- the models ----
   const bool single = (modelCount == 1);
   if(!single){printf("Read in %d models...\n", modelCount); fflush(stdout);}
   const time_t readStart = time(NULL);
   time_t lastReadReport = readStart;
   int modelsRead = 0;
   // the models are independent of each other: with several threads each thread reads whole
   // models (dynamic schedule, since the sizes differ a lot); a single model is read serially
   const int threads = ResolveThreads(CoreCount);
   #pragma omp parallel for num_threads(threads) schedule(dynamic, 1) if(threads > 1 && !single)
   for(int k = 0; k < modelCount; k++){
      ScoreModel & md = models[k];
      if(single){printf("Read in %s%s...\n", md.file.c_str(), md.binary? " (binary model)": ""); fflush(stdout);}
      // a model that cannot be read is skipped with a note (an exception cannot leave a parallel
      // loop, so it is caught here); a single model that fails ends the program as before
      try{
         if(md.binary) ReadBinaryModelFile(md);
         else ReadModelFile(md);
      }catch(const ModelError & e){
         md.failed = e.what();
         if(md.failed.empty()) md.failed = "unknown error";
         std::vector<int>().swap(md.rowIdx);
         std::vector<double>().swap(md.weights);
         std::vector<char>().swap(md.switched);
         std::vector<double>().swap(md.freq);
      }
      if(single && !md.failed.empty()) error("%s", md.failed.c_str());
      if(!single){
         #pragma omp critical(secretscores_progress)
         {
            modelsRead++;
            if(!md.failed.empty()){
               printf("  WARNING - model %s is skipped: %s\n", md.file.c_str(), OneLine(md.failed).c_str());
               fflush(stdout);
            }
            if(time(NULL) - lastReadReport >= 30){
               lastReadReport = time(NULL);
               printf("  %d of %d models read (%.0f s)\n", modelsRead, modelCount, difftime(lastReadReport, readStart));
               fflush(stdout);
            }
         }
      }
      if(single){
         printf("  %d SNPs in the model, %d of them match the genotype data by chromosome, position and alleles\n",
            md.fileSNPs, (int) md.rowIdx.size());
         if(md.malformed)
            printf("  %d lines without a usable effect_allele and effect_weight were skipped\n", md.malformed);
         if(md.noCoord)
            printf("  %d SNPs have no usable hm_chr / hm_pos\n", md.noCoord);
         if(md.noOtherAllele)
            printf("  %d SNPs have no other allele (other_allele / hm_inferOtherAllele) and are dropped\n", md.noOtherAllele);
         if(md.sameAllele)
            printf("  %d SNPs have the same other allele as effect allele and are dropped\n", md.sameAllele);
         if(md.nonAutosomal)
            printf("  %d SNPs are on chromosomes that are not scored (%s) and are skipped\n", md.nonAutosomal, NotScoredChromosomes());
         if(md.notInData)
            printf("  %d SNPs are not in the genotype data (chromosome and position)\n", md.notInData);
         if(md.alleleMismatch)
            printf("  %d SNPs are in the genotype data but the alleles do not match\n", md.alleleMismatch);
         if(md.rowIdx.size() == 0) error("No SNPs were extracted");
      }
   }

   // models that could not be read are left out of everything that follows; they are listed in
   // <prefix>_scores.skipped
   {
      std::vector<std::string> skippedFiles, skippedReasons;
      int kept = 0;
      for(int k = 0; k < modelCount; k++){
         if(models[k].failed.empty()){
            if(kept != k) std::swap(models[kept], models[k]);
            kept++;
         }else{
            skippedFiles.push_back(models[k].file);
            skippedReasons.push_back(models[k].failed);
         }
      }
      if(!skippedFiles.empty()){
         const std::string skippedfile = prefix + "_scores.skipped";
         WriteSkippedFile(skippedfile, skippedFiles, skippedReasons);
         printf("  %d of %d models could not be read and are skipped; see %s\n", (int) skippedFiles.size(),
            modelCount, skippedfile.c_str());
         if(kept == 0) error("None of the %d models could be read.", modelCount);
         modelCount = kept;
      }
   }

   if(!single){
      long long matched = 0;
      for(int k = 0; k < modelCount; k++) matched += (long long) models[k].rowIdx.size();
      printf("  %d models read: %lld SNPs match the genotype data in total (%.0f s)\n", modelCount, matched,
         difftime(time(NULL), readStart));
      fflush(stdout);
   }

   // ---- options in effect ----
   printf("\nOptions in effect:\n");
   if(!optionReadModel.empty())
      printf("\t--readmodel %s\n", optionReadModel.c_str());
   if(!optionReadList.empty())
      printf("\t--readlist %s\n", optionReadList.c_str());
   if(!readFreqFile.empty())
      printf("\t--readfreq %s%s\n", readFreqFile.c_str(), MeanImputation? "": " (not used: --no-mean-imputation)");
   if(!MeanImputation)
      printf("\t--no-mean-imputation\n");
   if(CoreCount > 0)
      printf("\t--cpus %d\n", CoreCount);
   if(XChrModel != 0)
      printf("\t--allowX %d\n", XChrModel);
   if(readBlockMB > 0)
      printf("\t--blocksize %d\n", readBlockMB);
   if(prefix != "secretscores")
      printf("\t--prefix %s\n", prefix.c_str());
   printf("\n%d model%s, %d thread%s\n\n", modelCount, modelCount == 1? "": "s", threads, threads == 1? "": "s");
   fflush(stdout);

   // ---- the X chromosome: the sex of the people matters for the SNPs of the models on it ----
   if(XChrModel != 0){
      long long xEntries = 0;
      for(int k = 0; k < modelCount; k++)
         for(size_t m = 0; m < models[k].rowIdx.size(); m++)
            if(bimChr[models[k].rowIdx[m]] == SEXCHR) xEntries++;
      if(xEntries > 0){
         int males = 0, females = 0, unknown = 0;
         for(int i = 0; i < idCount; i++){
            if(sex[i] == 1) males++;
            else if(sex[i] == 2) females++;
            else unknown++;
         }
         printf("X chromosome: %lld SNPs of the models are on it and are scored (--allowX %d: %s)\n",
            xEntries, XChrModel, XChrModel == 1? "males 0/1, females 0/1/2": "males and females 0/1/2");
         printf("  Sex in %s.fam: %d male%s, %d female%s, %d unknown (scored as females)\n", plinkinput.c_str(),
            males, males == 1? "": "s", females, females == 1? "": "s", unknown);
         if(males == 0)
            printf("  Warning: no sample is male in the .fam file, so no heterozygous genotype can be recognized as\n"
                   "  impossible at X SNPs; all the samples are scored as diploid, as stored in the .bed file\n");
         printf("  The pseudo-autosomal regions need to be coded XY (25) in the .bim file (plink --split-x);\n"
                "  heterozygous genotypes of males at X SNPs (chromosome 23) are taken as missing\n\n");
         fflush(stdout);
      }
   }

   // ---- memory: half of the RAM that is available now holds genotypes ----
   {
      double available = AvailableMemory();
      if(available > 0){
         GENOTYPE_MEMORY_LIMIT = (long long)(available / 2);
         String a = MemoryInfo(available);
         String b = MemoryInfo((double) GENOTYPE_MEMORY_LIMIT);
         printf("Available RAM is %s; GENOTYPE_MEMORY_LIMIT is set to half of it, %s\n",
            (const char*) a, (const char*) b);
      }else{
         String b = MemoryInfo((double) GENOTYPE_MEMORY_LIMIT);
         printf("The available RAM cannot be detected; GENOTYPE_MEMORY_LIMIT is %s\n", (const char*) b);
      }
   }

   // ---- the .bed file: SNP-major or individual-major, kept in memory individual-major ----
   threadCount = threads;
   FILE *fp = OpenPLINKBinary();

   // ---- score the models, batch by batch ----
   slotOfRow.assign(bimCount, -1);
   std::vector<char> needed(bimCount, 0);
   std::vector<std::vector<double> > scores(modelCount, std::vector<double>(idCount, 0.0));
   std::vector<double> missSNP;                // only reported for a single model
   if(single) missSNP.assign(idCount, 0.0);
   int batchStart = 0, batchNumber = 0;
   while(batchStart < modelCount){
      // Models are added to the batch until it has maxBatchModels models, or has at least
      // minModels models and the union of their SNPs ("touched") has reached
      // TARGET_TOUCHED_SNPS_PER_BATCH SNPs (the model that reaches it is kept).  The memory limit
      // comes first: a model that would make the genotypes (idCount people x ceil(SNPs/4) bytes)
      // exceed GENOTYPE_MEMORY_LIMIT is left for the next batch, even below minModels.
      const int minModels = std::min(MIN_MODELS_PER_BATCH, maxBatchModels);
      std::vector<int> touched;
      int batchEnd = batchStart;
      while(batchEnd < modelCount && batchEnd - batchStart < maxBatchModels){
         size_t before = touched.size();
         std::vector<int> & rows = models[batchEnd].rowIdx;
         for(size_t m = 0; m < rows.size(); m++)
            if(!needed[rows[m]]){
               needed[rows[m]] = 1;
               touched.push_back(rows[m]);
            }
         if(batchEnd > batchStart &&
            (long long)((touched.size() + 3) / 4) * idCount > GENOTYPE_MEMORY_LIMIT){
            for(size_t t = before; t < touched.size(); t++)
               needed[touched[t]] = 0;
            touched.resize(before);
            break;
         }
         batchEnd++;
         if(batchEnd - batchStart >= minModels && touched.size() >= (size_t) TARGET_TOUCHED_SNPS_PER_BATCH)
            break;
      }
      if((long long)((touched.size() + 3) / 4) * idCount > GENOTYPE_MEMORY_LIMIT){
         String need = MemoryInfo((double)((touched.size() + 3) / 4) * idCount);
         String lim = MemoryInfo((double) GENOTYPE_MEMORY_LIMIT);
         printf("Warning: the genotypes of %s need %s, more than GENOTYPE_MEMORY_LIMIT (%s); loading them anyway.\n",
            models[batchStart].name.c_str(), (const char*) need, (const char*) lim);
      }
      std::vector<int> sortedRows(touched);
      std::sort(sortedRows.begin(), sortedRows.end());
      for(size_t t = 0; t < sortedRows.size(); t++)
         slotOfRow[sortedRows[t]] = (int) t;
      batchNumber++;
      if(!single){
         printf("Batch %d: models %d-%d, loading the genotypes of %d SNPs for %d individuals...\n",
            batchNumber, batchStart + 1, batchEnd, (int) sortedRows.size(), idCount);
         fflush(stdout);
      }
      const time_t loadStart = time(NULL);
      ReadPLINKBinary(fp, sortedRows);
      if(single)
         printf("  PLINK binary genotypes at %d SNPs are loaded.\n", loadedCount);
      else
         printf("  Genotypes loaded (%.0f s)\n", difftime(time(NULL), loadStart));
      long long entries = 0;
      for(int k = batchStart; k < batchEnd; k++){
         if(single) printf("Compute polygenic scores...\n");
         PrepareModel(models[k], single);
         entries += models[k].used;
      }
      // right before the scores are computed: say so, now (the computing can take very long)
      {
         String tables = MemoryInfo((double) entries * (sizeof(int) + 4 * sizeof(double)));
         printf("Computing the scores of %d model%s (%lld SNP entries, model tables %s) for %d individuals, %d thread%s...\n",
            batchEnd - batchStart, batchEnd - batchStart == 1? "": "s", entries, (const char*) tables,
            idCount, threads, threads == 1? "": "s");
         fflush(stdout);
      }
      const time_t scoreStart = time(NULL);
      ScoreBatch(models, batchStart, batchEnd, scores, single? &missSNP: NULL);
      if(!single){
         printf("  Batch %d scored (%.0f s)\n", batchNumber, difftime(time(NULL), scoreStart));
         fflush(stdout);
      }
      FreeGenotypes(sortedRows);
      // the tables of the batch are not needed any more
      for(int k = batchStart; k < batchEnd; k++){
         std::vector<int>().swap(models[k].vslot);
         std::vector<double>().swap(models[k].vval);
         std::vector<int>().swap(models[k].xslot);
         std::vector<double>().swap(models[k].xval);
         std::vector<double>().swap(models[k].xvalM);
      }
      for(size_t t = 0; t < touched.size(); t++) needed[touched[t]] = 0;
      batchStart = batchEnd;
   }
   fclose(fp);
   if(maleHetX > 0)
      printf("%lld heterozygous genotypes of males at X-chromosome SNPs (outside the PAR) were taken as missing\n"
             "  (counted for each SNP of each model; if there are many, code the PAR as XY with plink --split-x)\n", maleHetX);

   int scoredCount = 0;
   for(int k = 0; k < modelCount; k++)
      if(models[k].scored) scoredCount++;
   if(!MeanImputation)
      printf("Mean imputation is off (--no-mean-imputation): missing genotypes count as 0\n");
   else{
      long long fromFile = 0, fromModel = 0, none = 0;
      for(int k = 0; k < modelCount; k++){
         fromFile += models[k].freqFromFile;
         fromModel += models[k].freqFromModel;
         none += models[k].freqNone;
      }
      printf("Mean imputation of missing genotypes (2 x effect allele frequency): %lld SNP entries use the frequency from --readfreq, %lld from\n"
             "  the allelefrequency_effect of the model, %lld have none (missing genotypes count as 0 there)\n", fromFile, fromModel, none);
      if(none > 0 && fromFile + fromModel == 0)
         printf("  Note: no allele frequencies are available (use --readfreq, or a model with allelefrequency_effect)\n");
   }
   if(single){
      if(!models[0].scored)
         error("None of the %d SNPs of the model could be scored", (int) models[0].rowIdx.size());
      printf("  %d SNPs are scored\n", models[0].used);
      printf("Note: Action OK treats allele A1 (in the data) as the effect allele (in the model)\n");
      printf("      Action SWITCHED treats allele A2 (in the data) as the effect allele (in the model)\n\n");
   }else{
      printf("\n%15s %9s %9s %9s %9s\n", "SCORE", "SNPs", "Found", "Used", "Coverage");
      int listed = 0;
      for(int k = 0; k < modelCount; k++){
         if(listed < 20 || !models[k].scored){
            const int found = (int) models[k].rowIdx.size() + models[k].alleleMismatch;
            printf("%15s %9d %9d %9d %9.4f%s\n", models[k].name.c_str(), models[k].fileSNPs,
               found, models[k].used, models[k].fileSNPs? (double) models[k].used / models[k].fileSNPs: 0.0,
               models[k].scored? "": "   (no SNP scored: not written)");
         }
         listed++;
      }
      if(modelCount > 20) printf("  ... %d more models; see the coverage file\n", modelCount - 20);
      if(scoredCount == 0)
         error("None of the %d models could be scored (no SNPs found)", modelCount);
      printf("%d of %d models were scored\n\n", scoredCount, modelCount);
   }

   // ---- one line per model: SNPs = SNPs of the model, Found = of them at a chromosome and
   //      position of the data, Used = of them also with matching alleles (scored),
   //      coverage = Used / SNPs ----
   std::string coveragefile = prefix + "_scores.coverage";
   {
      std::vector<CoverageRow> rows(modelCount);
      for(int k = 0; k < modelCount; k++){
         CoverageRow & r = rows[k];
         r.name = models[k].name;
         r.file = models[k].file;
         r.trait = models[k].trait;
         r.traitEfo = models[k].traitEfo;
         r.fileSNPs = models[k].fileSNPs;
         r.found = (long long) models[k].rowIdx.size() + models[k].alleleMismatch;
         r.used = models[k].used;
         if(models[k].scored)
            r.percentImputed = 100.0 * (double) models[k].imputedGenotypes / ((double) models[k].used * idCount);
         r.noCoord = models[k].noCoord;
         r.noOtherAllele = models[k].noOtherAllele;
         r.sameAllele = models[k].sameAllele;
         r.nonAutosomal = models[k].nonAutosomal;
         r.notInData = models[k].notInData;
         r.alleleMismatch = models[k].alleleMismatch;
      }
      WriteCoverageFile(coveragefile, rows, true);
   }

   // ---- the score file for SECRET: IID, then one column per model, people in .fam order ----
   // SECRET matches the first column to the IID of the phenotype data.  When IIDs repeat
   // across families, FID_IID is written (SECRET uses that form when the IDs have a _).
   bool useFID = false;
   {
      std::unordered_set<std::string> seen;
      for(int i = 0; i < idCount; i++)
         if(!seen.insert(IID[i]).second){useFID = true; break;}
   }
   if(useFID)
      printf("Some IIDs appear in more than one family: IDs are written as FID_IID.\n");
   std::string scorefile = prefix + "_scores.txt";
   fp = fopen(scorefile.c_str(), "wt");
   if(fp == NULL) error("Cannot open %s to write", scorefile.c_str());
   fprintf(fp, "IID");
   for(int k = 0; k < modelCount; k++)
      if(models[k].scored) fprintf(fp, "\t%s", models[k].name.c_str());
   fprintf(fp, "\n");
   for(int i = 0; i < idCount; i++){
      if(useFID) fprintf(fp, "%s_%s", FID[i].c_str(), IID[i].c_str());
      else fprintf(fp, "%s", IID[i].c_str());
      for(int k = 0; k < modelCount; k++)
         if(models[k].scored) fprintf(fp, "\t%.8g", scores[k][i]);
      fprintf(fp, "\n");
   }
   fclose(fp);

   if(single){
      // ---- the score of one model, with the fraction of its SNPs that were genotyped ----
      std::string riskscorefile = prefix + "grs.txt";
      fp = fopen(riskscorefile.c_str(), "wt");
      if(fp == NULL) error("Cannot open %s to write", riskscorefile.c_str());
      fprintf(fp, "%15s %15s %7s %14s\n", "FID", "IID", "InfoSNP", "Score");
      for(int i = 0; i < idCount; i++)
         fprintf(fp, "%15s %15s %7.3lf %14.8g\n", FID[i].c_str(), IID[i].c_str(), 1 - missSNP[i], scores[0][i]);
      fclose(fp);
      printf("Polygenic scores are saved in file %s\n", riskscorefile.c_str());
   }
   printf("Scores of %d model%s for %d samples are saved in file %s (the score file for SECRET)\n",
      scoredCount, scoredCount==1? "": "s", idCount, scorefile.c_str());
   printf("Coverage of the models is saved in file %s\n", coveragefile.c_str());
   delete []models;
   printf("Polygenic score calculation ends at %s", currentTime());
}
