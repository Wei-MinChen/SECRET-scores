////////////////////////////////////////////////////////////////////////
// BinaryInput.cpp
// (c) 2024-2026 Wei-Min Chen
//
// This file is distributed as part of the SECRETS source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile SECRETS.
//
// All computer programs have bugs. Use this file at your own risk.
//
// SECRETS binary models (.bin) for secretscores: made from PGS Catalog text
// models by --makebin (no genotypes needed), read back for scoring.
//
// ---- The format, version 4 (little-endian; version 3 is the same without frequencies) ----
//
//   HEADER (headerSize bytes; the data start at byte headerSize):
//   byte 0         version = 4 (a polygenic score model with allele frequencies), or 3 (without)
//   bytes 1..4     uint32 headerSize: the size of the whole header in bytes, counted from
//                  the first byte of the file (version byte and this field included).
//                  It is 109 without a trait, and 109 + 8 + L with one.
//   bytes 5..108   26 x uint32: number of SNPs on chromosomes 1..26
//                  (23 = X, 24 = Y, 25 = XY, 26 = MT)
//   optional trait block (bytes 109..headerSize-1; written when the text model has a trait
//   annotation in its header lines; a header without it is still valid):
//      4 bytes "TRAI", uint32 length L, L bytes of text: the trait (trait_reported, else
//      trait_mapped of the PGS Catalog header; empty when it has neither), then, when there is
//      a trait_efo, a tab and the trait_efo.  A header with only trait_efo is stored as
//      "<TAB>trait_efo" (empty trait before the tab).  It is not used for scoring; it goes to the "trait" and "trait_efo"
//      columns of the coverage file.
//   Because the trait block is inside the header, whose size is given in bytes 1..4, nothing
//   after the header can be mistaken for it: the data are exactly 8 bytes per SNP (version 3) or
//   12 bytes per SNP (version 4), and the file size must be headerSize + 8 x (total number of
//   SNPs) for version 3 and headerSize + 12 x (total number of SNPs) for version 4.
//
//   DATA, for each chromosome that has SNPs, in order of chromosome:
//      count x uint32    position << 4 | allele code
//      count x float32   effect_weight
//      count x float32   frequency of the effect allele (version 4 only; right after the weights;
//                        the allelefrequency_effect column of the text model; NaN when a SNP has none)
//
//   allele code (4 bits) = effect allele * 4 + other allele, with A = 0, T = 1,
//   C = 2, G = 3.  The 12 codes with two different alleles give the ordered pair
//   (effect allele, other allele); every SNP has both alleles, and they differ.
//   The 4 codes with the same letter twice (0, 5, 10, 15) are not valid: a reader
//   drops such an entry (an early version 2 file used them for "no other allele").
//   Positions need 28 bits (up to 268,435,455).  Within a chromosome the SNPs are
//   in order of position (ties keep their order in the text model).
//   The weight is exactly the effect_weight of the text model, the weight of the
//   effect allele; it is never negated, and no SNP is dropped because its weight
//   is zero.  Only SNPs whose alleles are single letters A, C, G, T can be stored.
//
//   --makebin writes version 4 when the text model has an allelefrequency_effect column, else
//   version 3.  The frequencies are used for mean imputation of missing genotypes (see below).
//
//   Version 2 (secretscores 1.0.0 and 1.0.1) is the same without the headerSize field:
//   byte 0 version = 2, bytes 1..104 the 26 SNP counts, the data from byte 105, and the
//   trait block, if any, AFTER the data (found only from the file size).  It is still read,
//   but --makebin writes version 3.
//
//   Versions 0 and 1 of earlier secretscores/SECRETS stored the pair in a fixed
//   order and negated the weight when the effect allele came second.  That loses
//   which allele the weight belongs to, so those files cannot give PGS Catalog
//   scores; they are refused, and have to be made again from the text models.
//
// ---- --makebin: text model -> binary model ----
//   For every line of the PGS Catalog scoring file with a usable effect_allele,
//   effect_weight, hm_chr and hm_pos (the same lines that are scored):
//     1. drop it when it has no other allele (other_allele, else hm_inferOtherAllele),
//        or when the other allele is the effect allele;
//     2. skip it when the effect allele or the other allele is not one letter from
//        A, C, G, T, or when the position is above 268,435,455 (all counted);
//     3. code = effect allele * 4 + other allele, and keep (position << 4 | code,
//        weight as float32);
//     4. put the SNPs in chromosome order and, within a chromosome, in order of
//        position (stable), count them, and write the file.
//
// ---- Scoring with a model (text or binary: the same algorithm) ----
//   For every SNP of the model, find the SNP of the genotype data at the same
//   chromosome and position (.bim) whose alleles are the effect and other allele:
//     effect allele = A1 and other allele = A2   ->  OK
//     effect allele = A2 and other allele = A1   ->  SWITCHED
//   Alleles are compared as they are: no strand flip, no frequencies.  A model
//   SNP without an other allele (different from the effect allele) is dropped.
//   With c the number of effect alleles of a person at the SNP (OK: copies of
//   A1; SWITCHED: copies of A2 = 2 - copies of A1), and a missing genotype
//   counted as 0:
//        score = sum over the SNPs of  effect_weight x c
//   which is the PGS Catalog definition.  Weights are float32 in a binary model,
//   about 7 significant digits, so scores can differ from the text model's in
//   the 7th digit.
//   Missing genotypes: by default (mean imputation) a missing genotype counts as 2 x the
//   frequency of the effect allele (--readfreq gives the frequency of the data, which is used
//   first; else the frequency of the model), with no frequency as 0; --no-mean-imputation counts
//   every missing genotype as 0.
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <limits>
#include <vector>
#include "ScoreEngine.h"
#include "Error.h"

static const int BIN_VERSION = 4;                     // written by --makebin when the model has frequencies
static const int BIN_VERSION_NOFREQ = 3;              // written by --makebin without frequencies
static const int BIN_VERSION_LEGACY = 2;              // still read: no header size, trailer after the data
static const long long FIXED_HEADER_SIZE = 1 + 4 + 26 * 4;     // version, headerSize, 26 SNP counts = 109
static const char LETTERS[4] = { 'A', 'T', 'C', 'G' };
static const long MAX_BIN_POSITION = 0xFFFFFFFL;      // 28 bits

bool IsBinaryModelName(const std::string & file)
{
   if(file.size() <= 4) return false;
   const char *e = file.c_str() + file.size() - 4;
   return e[0] == '.' && (e[1]=='b' || e[1]=='B') && (e[2]=='i' || e[2]=='I') && (e[3]=='n' || e[3]=='N');
}

static int LetterCode(const std::string & allele)
{
   if(allele.size() != 1) return -1;
   switch(allele[0]){
      case 'A': case 'a': return 0;
      case 'T': case 't': return 1;
      case 'C': case 'c': return 2;
      case 'G': case 'g': return 3;
   }
   return -1;
}

// Errors of the binary model reader go to the screen (stdout) as everywhere in secretscores,
// but to stderr for --viewbin, whose stdout is the text model (e.g. redirected to a file).
static bool errorToStderr = false;

// Closes a file when a ModelError leaves the function that opened it.
struct FileCloser
{
   FILE *fp;
   explicit FileCloser(FILE *f) : fp(f) {}
   ~FileCloser() { if(fp) fclose(fp); }
   FILE *Release() { FILE *f = fp; fp = NULL; return f; }
};

static void Fail(const char *format, ...)
{
   char message[2048];
   va_list ap;
   va_start(ap, format);
   vsnprintf(message, sizeof(message), format, ap);
   va_end(ap);
   if(errorToStderr){
      fprintf(stderr, "\nFATAL ERROR - \n%s\n\n", message);
      exit(EXIT_FAILURE);
   }
   ThrowModelError("%s", message);      // a damaged model: skipped, the others go on
}

// Read the trait block "TRAI", uint32 L, L bytes of text (trait, tab, trait_efo) at the
// current position of fp, which must hold exactly `available` bytes.
static void ReadTraitBlock(FILE *fp, long long available, const char *filename, std::string & trait,
   std::string & traitEfo)
{
   char magic[4];
   unsigned int len = 0;
   if(available < 8 || fread(magic, 1, 4, fp) < 4 || memcmp(magic, "TRAI", 4) != 0
      || fread(&len, sizeof(unsigned int), 1, fp) < 1 || (long long) len != available - 8)
      Fail("The trait block in the header of %s is damaged", filename);
   std::vector<char> text(len + 1, 0);
   if(len > 0 && fread(&text[0], 1, len, fp) < len)
      Fail("Not enough bytes for the trait in %s", filename);
   const std::string all(&text[0], len);
   const size_t tab = all.find('\t');
   trait = all.substr(0, tab);
   if(tab != std::string::npos) traitEfo = all.substr(tab + 1);
}

// Open a binary model file and check it: the version, the header size, the SNP counts, and a
// size that agrees with them (a damaged file must not drive allocations).  Returns the file
// positioned at the first chromosome's data; snpCounts are the numbers of SNPs on chromosomes
// 1..26; trait and traitEfo come from the optional trait block (empty without one).
static FILE *OpenBinaryModel(const char *filename, unsigned int snpCounts[26], std::string & trait,
   std::string & traitEfo, bool & hasFreq)
{
   FILE *fp = fopen(filename, "rb");
   if(fp == NULL) Fail("Cannot open %s to read", filename);
   FileCloser guard(fp);
   int version = fgetc(fp);
   if(version == 0 || version == 1)
      Fail("%s is a binary model of an earlier format (version %d), which does not keep the weight\n"
         "  of the effect allele. Make it again from its text model with --makebin.", filename, version);
   if(version != BIN_VERSION && version != BIN_VERSION_NOFREQ && version != BIN_VERSION_LEGACY)
      Fail("The SECRETS format in %s cannot be recognized. The first byte is %x", filename, version);
   trait.clear();
   traitEfo.clear();
   hasFreq = (version == BIN_VERSION);
   const long long bytesPerSNP = hasFreq? 12: 8;
   FSEEK64(fp, 0, SEEK_END);
   const long long size = FTELL64(fp);

   if(version == BIN_VERSION_LEGACY){
      // version 2: SNP counts at byte 1, the data at byte 105, the optional trailer after the data
      const long long dataStart = 1 + 26 * 4;
      FSEEK64(fp, 1, SEEK_SET);
      if(fread(snpCounts, sizeof(unsigned int), 26, fp) < 26)
         Fail("Not enough words for SNP counts at each chromosome in %s", filename);
      long long expected = dataStart;
      for(int c = 0; c < 26; c++) expected += 8LL * snpCounts[c];
      if(size < expected)
         Fail("The size of %s (%lld bytes) does not agree with its SNP counts (%lld bytes expected)",
            filename, size, expected);
      if(size > expected){
         FSEEK64(fp, expected, SEEK_SET);
         ReadTraitBlock(fp, size - expected, filename, trait, traitEfo);
      }
      FSEEK64(fp, dataStart, SEEK_SET);
      return guard.Release();
   }

   // version 3: bytes 1..4 give the size of the header, which ends with the optional trait block
   unsigned int headerSize = 0;
   FSEEK64(fp, 1, SEEK_SET);
   if(fread(&headerSize, sizeof(unsigned int), 1, fp) < 1)
      Fail("Not enough bytes for the header size in %s", filename);
   if((long long) headerSize < FIXED_HEADER_SIZE)
      Fail("The header size in %s (%u bytes) is smaller than the minimum of %lld bytes", filename,
         headerSize, FIXED_HEADER_SIZE);
   if((long long) headerSize > size)
      Fail("The header size in %s (%u bytes) is larger than the file (%lld bytes)", filename,
         headerSize, size);
   if(fread(snpCounts, sizeof(unsigned int), 26, fp) < 26)
      Fail("Not enough words for SNP counts at each chromosome in %s", filename);
   long long expected = headerSize;
   for(int c = 0; c < 26; c++) expected += bytesPerSNP * snpCounts[c];
   if(size != expected)
      Fail("The size of %s (%lld bytes) does not agree with its header size and SNP counts (%lld bytes expected)",
         filename, size, expected);
   if(headerSize > FIXED_HEADER_SIZE)
      ReadTraitBlock(fp, (long long) headerSize - FIXED_HEADER_SIZE, filename, trait, traitEfo);
   FSEEK64(fp, headerSize, SEEK_SET);
   return guard.Release();
}

// Read one binary model file and match its SNPs to the genotype data (chromosome, position,
// alleles), exactly as ReadModelFile does for a text model.
void ScoreEngine::ReadBinaryModelFile(ScoreModel & model)
{
   const char *filename = model.file.c_str();
   model.rowIdx.clear();
   model.weights.clear();
   model.switched.clear();
   model.fileSNPs = model.malformed = model.noCoord = model.nonAutosomal = 0;
   model.noOtherAllele = model.sameAllele = 0;
   model.notInData = model.alleleMismatch = 0;
   model.freq.clear();
   unsigned int snpCounts[26];
   bool hasFreq = false;
   FILE *fp = OpenBinaryModel(filename, snpCounts, model.trait, model.traitEfo, hasFreq);
   FileCloser guard(fp);          // closed at the end, or when a damaged file throws

   std::vector<unsigned int> pos;
   std::vector<float> wt, fq;
   const double NaN = std::numeric_limits<double>::quiet_NaN();
   for(int c = 0; c < 26; c++){
      const size_t count = snpCounts[c];
      if(count == 0) continue;
      pos.resize(count);
      wt.resize(count);
      if(fread(&pos[0], sizeof(unsigned int), count, fp) < count)
         ThrowModelError("Not enough bytes for positions/alleles at the %dth chromosome in %s", c + 1, filename);
      if(fread(&wt[0], sizeof(float), count, fp) < count)
         ThrowModelError("Not enough bytes for weights at the %dth chromosome in %s", c + 1, filename);
      if(hasFreq){
         // the frequencies are only read for mean imputation; otherwise they are skipped
         if(MeanImputation){
            fq.resize(count);
            if(fread(&fq[0], sizeof(float), count, fp) < count)
               ThrowModelError("Not enough bytes for allele frequencies at the %dth chromosome in %s", c + 1, filename);
         }else FSEEK64(fp, (long long) count * (long long) sizeof(float), SEEK_CUR);
      }
      const int chr = c + 1;
      if(!ScoredChromosome(chr)){
         model.fileSNPs += (int) count;
         model.nonAutosomal += (int) count;
         continue;
      }
      for(size_t m = 0; m < count; m++){
         const double w = wt[m];
         const unsigned int code = pos[m] & 0xF;
         if(!(fabs(w) < 1e30)){model.malformed++; continue;}     // also catches NaN
         model.fileSNPs++;
         if((code >> 2) == (code & 3)){model.noOtherAllele++; continue;}   // not a valid code
         const long position = (long)(pos[m] >> 4);
         if(position <= 0){model.noCoord++; continue;}
         const unsigned int ea = (unsigned int)(unsigned char) LETTERS[code >> 2];
         const unsigned int oa = (unsigned int)(unsigned char) LETTERS[code & 3];
         bool sw = false, atPosition = false;
         const int hit = MatchSNP(chr, position, ea, oa, sw, atPosition);
         if(hit == -1){
            if(atPosition) model.alleleMismatch++;
            else model.notInData++;
            continue;
         }
         model.rowIdx.push_back(hit);
         model.weights.push_back(w);
         model.switched.push_back(sw? 1: 0);
         const double af = (hasFreq && MeanImputation)? (double) fq[m]: NaN;
         model.freq.push_back((af >= 0.0 && af <= 1.0)? af: NaN);
      }
   }
}

// ---------------------------------------------------------------------------------------
// --makebin
// ---------------------------------------------------------------------------------------

// Append printf-style text to a string (a message is built whole and printed at once, so that
// the messages of models converted by different threads do not mix).
static void Appendf(std::string & out, const char *format, ...)
{
   char buffer[1024];
   va_list ap;
   va_start(ap, format);
   int n = vsnprintf(buffer, sizeof(buffer), format, ap);
   va_end(ap);
   if(n < 0) return;
   if((size_t) n < sizeof(buffer)){out.append(buffer, n); return;}
   std::vector<char> big((size_t) n + 1);
   va_start(ap, format);
   vsnprintf(&big[0], big.size(), format, ap);
   va_end(ap);
   out.append(&big[0], n);
}

struct BinSNP
{
   unsigned int item;            // position << 4 | allele code
   float weight;
   float freq;                   // frequency of the effect allele, NaN: none
};

static bool ByPosition(const BinSNP & a, const BinSNP & b) { return (a.item >> 4) < (b.item >> 4); }

static void MakeOneBinaryModel(const std::string & textFile, const std::string & outFile, const std::string & name,
   CoverageRow & cover)
{
   std::vector<BinSNP> snps[26];
   TextModelStats stats;
   int unsupportedAlleles = 0, positionTooLarge = 0, nonAutosomal = 0;
   ParseTextModel(textFile, stats,
      [&](int chr, long pos, const std::string & ea, const std::string & oa, double wt, double af){
         if(!ScoredChromosome(chr)) nonAutosomal++;     // written, but not scored
         const int e = LetterCode(ea), o = LetterCode(oa);       // both alleles exist and differ
         if(e < 0 || o < 0 || e == o || !(fabs(wt) < 3.0e38)){unsupportedAlleles++; return;}
         if(pos > MAX_BIN_POSITION){positionTooLarge++; return;}
         BinSNP s;
         s.item = ((unsigned int) pos << 4) | (unsigned int)(e << 2 | o);
         s.weight = (float) wt;
         s.freq = (af >= 0.0 && af <= 1.0)? (float) af: std::numeric_limits<float>::quiet_NaN();
         snps[chr - 1].push_back(s);
      });
   unsigned int counts[26];
   long long written = 0, withFreq = 0;
   const bool hasFreq = stats.hasFreq;       // the model has an allelefrequency_effect column: version 4
   for(int c = 0; c < 26; c++){
      for(size_t m = 0; m < snps[c].size(); m++) if(snps[c][m].freq == snps[c][m].freq) withFreq++;
      std::stable_sort(snps[c].begin(), snps[c].end(), ByPosition);
      counts[c] = (unsigned int) snps[c].size();
      written += counts[c];
   }
   // the trait annotation: the trait block of the header
   // the trait is trait_reported, else trait_mapped, else empty: it does NOT fall back to
   // trait_efo, so a model with only trait_efo is stored as "\tEFO_xxxx" (empty trait)
   std::string trait = !stats.traitReported.empty()? stats.traitReported: stats.traitMapped;
   std::string efo = stats.traitEfo;
   if(trait.size() > 4096) trait.resize(4096);
   if(efo.size() > 4096) efo.resize(4096);
   // trait\tefo: a tab in either would break the split
   for(size_t i = 0; i < trait.size(); i++) if(trait[i] == '\t') trait[i] = ' ';
   for(size_t i = 0; i < efo.size(); i++) if(efo[i] == '\t') efo[i] = ' ';
   std::string traitText;
   bool haveTrait = !trait.empty() || !efo.empty();
   if(haveTrait) traitText = efo.empty()? trait: trait + "\t" + efo;
   const unsigned int traitLen = (unsigned int) traitText.size();
   const unsigned int headerSize = (unsigned int) FIXED_HEADER_SIZE + (haveTrait? 8 + traitLen: 0);
   const int version = hasFreq? BIN_VERSION: BIN_VERSION_NOFREQ;

   FILE *fp = fopen(outFile.c_str(), "wb");
   if(fp == NULL) ThrowModelError("Cannot open %s to write", outFile.c_str());
   fputc(version, fp);
   fwrite(&headerSize, sizeof(unsigned int), 1, fp);
   fwrite(counts, sizeof(unsigned int), 26, fp);
   if(haveTrait){
      fwrite("TRAI", 1, 4, fp);
      fwrite(&traitLen, sizeof(unsigned int), 1, fp);
      fwrite(traitText.c_str(), 1, traitLen, fp);
   }
   std::vector<unsigned int> item;
   std::vector<float> weight, freq;
   for(int c = 0; c < 26; c++){
      if(counts[c] == 0) continue;
      item.resize(counts[c]);
      weight.resize(counts[c]);
      freq.resize(counts[c]);
      for(size_t m = 0; m < counts[c]; m++){
         item[m] = snps[c][m].item;
         weight[m] = snps[c][m].weight;
         freq[m] = snps[c][m].freq;
      }
      fwrite(&item[0], sizeof(unsigned int), counts[c], fp);
      fwrite(&weight[0], sizeof(float), counts[c], fp);
      if(hasFreq) fwrite(&freq[0], sizeof(float), counts[c], fp);      // right after the weights
   }
   fclose(fp);
   cover.name = name;
   cover.file = textFile;
   cover.trait = trait;
   cover.traitEfo = efo;
   cover.fileSNPs = stats.fileSNPs;
   cover.noCoord = stats.noCoord;
   cover.noOtherAllele = stats.noOtherAllele;
   cover.sameAllele = stats.sameAllele;
   cover.nonAutosomal = nonAutosomal;
   std::string out;
   Appendf(out, "%s: %d SNPs read, %lld written to %s", textFile.c_str(), stats.fileSNPs, written, outFile.c_str());
   if(hasFreq) Appendf(out, " (with allele frequencies: %lld of %lld SNPs have one)", withFreq, written);
   if(!trait.empty()) Appendf(out, " (trait: %s%s%s)", trait.c_str(), efo.empty()? "": ", ", efo.c_str());
   else if(!efo.empty()) Appendf(out, " (trait_efo: %s)", efo.c_str());
   std::string notWritten;
   char part[160];
   if(stats.noCoord){snprintf(part, sizeof(part), "%s%d without hm_chr/hm_pos", notWritten.empty()? "": ", ", stats.noCoord); notWritten += part;}
   if(stats.noOtherAllele){snprintf(part, sizeof(part), "%s%d without other allele", notWritten.empty()? "": ", ", stats.noOtherAllele); notWritten += part;}
   if(stats.sameAllele){snprintf(part, sizeof(part), "%s%d with other allele = effect allele", notWritten.empty()? "": ", ", stats.sameAllele); notWritten += part;}
   if(unsupportedAlleles){snprintf(part, sizeof(part), "%s%d with alleles other than single A/C/G/T", notWritten.empty()? "": ", ", unsupportedAlleles); notWritten += part;}
   if(positionTooLarge){snprintf(part, sizeof(part), "%s%d with position above %ld", notWritten.empty()? "": ", ", positionTooLarge, MAX_BIN_POSITION); notWritten += part;}
   if(stats.malformed){snprintf(part, sizeof(part), "%s%d lines without effect_allele/effect_weight", notWritten.empty()? "": ", ", stats.malformed); notWritten += part;}
   if(!notWritten.empty()) Appendf(out, " (not written: %s)", notWritten.c_str());
   Appendf(out, "\n");
   if(written == 0) Appendf(out, "  Warning: %s has no SNPs.\n", outFile.c_str());
   #pragma omp critical(secretscores_output)
   {
      fputs(out.c_str(), stdout);
      fflush(stdout);
   }
}

void MakeBinaryModels(const std::vector<std::string> & files, const std::vector<std::string> & names,
   const std::string & prefix, int threads)
{
   const size_t n = files.size();
   std::vector<std::string> outFiles(n);
   for(size_t i = 0; i < n; i++){
      if(IsBinaryModelName(files[i])) continue;       // skipped below, nothing is written for it
      std::string name = names[i].empty()? ModelNameFromFile(files[i]): names[i];
      outFiles[i] = name + ".bin";
      for(size_t j = 0; j < n; j++)
         if(files[j] == outFiles[i])
            error("The binary model %s would overwrite the model file %s that is also used as input.\n  Give it a different name in the second column of the --readlist file.",
               outFiles[i].c_str(), files[j].c_str());
      for(size_t j = 0; j < i; j++)
         if(outFiles[j] == outFiles[i])
            error("Models %s and %s would both be written to %s.\n  Give them different names in the second column of the --readlist file.",
               files[j].c_str(), files[i].c_str(), outFiles[i].c_str());
   }
   std::vector<size_t> todo;
   for(size_t i = 0; i < n; i++){
      if(IsBinaryModelName(files[i])){
         printf("%s is a binary model already and is skipped.\n", files[i].c_str());
         continue;
      }
      todo.push_back(i);
   }
   const int total = (int) todo.size();
   std::vector<CoverageRow> allRows(todo.size());
   std::vector<std::string> failed(todo.size());
   // each model is converted by one thread; the dynamic schedule balances the very different sizes.
   // A model that cannot be converted is skipped with a note (an exception cannot leave a parallel
   // loop, so it is caught here); a single model that fails ends the program.
   #pragma omp parallel for num_threads(threads) schedule(dynamic, 1) if(threads > 1 && total > 1)
   for(int j = 0; j < total; j++){
      const size_t i = todo[j];
      const std::string name = outFiles[i].substr(0, outFiles[i].size() - 4);
      try{
         MakeOneBinaryModel(files[i], outFiles[i], name, allRows[j]);
      }catch(const ModelError & e){
         failed[j] = e.what();
         if(failed[j].empty()) failed[j] = "unknown error";
         if(total == 1) error("%s", failed[j].c_str());
         #pragma omp critical(secretscores_output)
         {
            printf("WARNING - model %s is skipped: %s\n", files[i].c_str(), OneLine(failed[j]).c_str());
            fflush(stdout);
         }
      }
   }
   // the models that were made, in the order of the input; the others go to the skipped file
   std::vector<CoverageRow> rows;
   std::vector<std::string> skippedFiles, skippedReasons;
   for(int j = 0; j < total; j++){
      if(failed[j].empty()) rows.push_back(allRows[j]);
      else{
         skippedFiles.push_back(files[todo[j]]);
         skippedReasons.push_back(failed[j]);
      }
   }
   const int made = (int) rows.size();
   if(!skippedFiles.empty()){
      const std::string skippedfile = prefix + "_scores.skipped";
      WriteSkippedFile(skippedfile, skippedFiles, skippedReasons);
      printf("%d of %d models could not be converted and are skipped; see %s\n", (int) skippedFiles.size(),
         total, skippedfile.c_str());
      if(made == 0) error("None of the %d models could be converted.", total);
   }
   printf("%d binary model%s made.\n", made, made == 1? "": "s");
   // one coverage file for all the models made, with the trait of each
   if(made > 0){
      const std::string coveragefile = prefix + "_scores.coverage";
      WriteCoverageFile(coveragefile, rows, false);
      printf("Coverage of the models is saved in file %s\n", coveragefile.c_str());
   }
}

// ---------------------------------------------------------------------------------------
// --viewbin: binary model -> text model on the screen (stdout), the reverse of --makebin
// ---------------------------------------------------------------------------------------

// The chromosome label of a code 1..26, as a text model writes it.
static const char *ChromosomeLabel(int chr, char *buffer)
{
   switch(chr){
      case 23: return "X";
      case 24: return "Y";
      case 25: return "XY";
      case 26: return "MT";
   }
   snprintf(buffer, 8, "%d", chr);
   return buffer;
}

void ViewBinaryModel(const std::string & file)
{
   errorToStderr = true;
   unsigned int snpCounts[26];
   std::string trait, traitEfo;
   bool hasFreq = false;
   FILE *fp = OpenBinaryModel(file.c_str(), snpCounts, trait, traitEfo, hasFreq);
   // header lines as in a PGS Catalog file, so that the text can be read back with --readmodel
   // and --makebin; the trait is written as trait_reported (the binary model keeps the one
   // trait that --makebin chose: trait_reported, else trait_mapped, else trait_efo)
   printf("###SECRETS binary model converted to a text model by secretscores --viewbin\n");
   if(!trait.empty()) printf("#trait_reported=%s\n", trait.c_str());
   if(!traitEfo.empty()) printf("#trait_efo=%s\n", traitEfo.c_str());
   printf("hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight%s\n", hasFreq? "\tallelefrequency_effect": "");
   std::vector<unsigned int> pos;
   std::vector<float> wt, fq;
   long long invalid = 0;
   for(int c = 0; c < 26; c++){
      const size_t count = snpCounts[c];
      if(count == 0) continue;
      pos.resize(count);
      wt.resize(count);
      if(fread(&pos[0], sizeof(unsigned int), count, fp) < count)
         Fail("Not enough bytes for positions/alleles at the %dth chromosome in %s", c + 1, file.c_str());
      if(fread(&wt[0], sizeof(float), count, fp) < count)
         Fail("Not enough bytes for weights at the %dth chromosome in %s", c + 1, file.c_str());
      if(hasFreq){
         fq.resize(count);
         if(fread(&fq[0], sizeof(float), count, fp) < count)
            Fail("Not enough bytes for allele frequencies at the %dth chromosome in %s", c + 1, file.c_str());
      }
      char label[8];
      const char *chr = ChromosomeLabel(c + 1, label);
      for(size_t m = 0; m < count; m++){
         const unsigned int code = pos[m] & 0xF;
         if((code >> 2) == (code & 3)){invalid++; continue;}     // not a valid code: dropped, as when scoring
         // %.9g: all the digits of a float32 weight, so that --makebin gives the same file back
         printf("%s\t%u\t%c\t%c\t%.9g", chr, pos[m] >> 4, LETTERS[code >> 2], LETTERS[code & 3], (double) wt[m]);
         if(hasFreq){
            if(fq[m] >= 0.0f && fq[m] <= 1.0f) printf("\t%.9g", (double) fq[m]);
            else printf("\tNA");
         }
         printf("\n");
      }
   }
   fclose(fp);
   if(invalid > 0)
      fprintf(stderr, "Note: %lld entries of %s with an invalid allele code were not written.\n", invalid, file.c_str());
}
