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
// ---- The format, version 2 (little-endian) ----
//
//   byte 0         version = 2 (a polygenic score model)
//   bytes 1..104   26 x uint32: number of SNPs on chromosomes 1..26
//                  (23 = X, 24 = Y, 25 = XY, 26 = MT)
//   then, for each chromosome that has SNPs, in order of chromosome:
//      count x uint32    position << 4 | allele code
//      count x float32   effect_weight
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
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "ScoreEngine.h"
#include "Error.h"

static const int BIN_VERSION = 2;
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

// Read one binary model file and match its SNPs to the genotype data (chromosome, position,
// alleles), exactly as ReadModelFile does for a text model.
void ScoreEngine::ReadBinaryModelFile(ScoreModel & model)
{
   const char *filename = model.file.c_str();
   FILE *fp = fopen(filename, "rb");
   if(fp == NULL) error("Cannot open %s to read", filename);
   model.rowIdx.clear();
   model.weights.clear();
   model.switched.clear();
   model.fileSNPs = model.malformed = model.noCoord = model.nonAutosomal = 0;
   model.noOtherAllele = model.sameAllele = 0;
   model.notInData = model.alleleMismatch = 0;

   int version = fgetc(fp);
   if(version == 0 || version == 1)
      error("%s is a binary model of an earlier format (version %d), which does not keep the weight\n"
         "  of the effect allele. Make it again from its text model with --makebin.", filename, version);
   if(version != BIN_VERSION)
      error("The SECRETS format in %s cannot be recognized. The first byte is %x", filename, version);
   unsigned int snpCounts[26];
   if(fread(snpCounts, sizeof(unsigned int), 26, fp) < 26)
      error("Not enough words for SNP counts at each chromosome in %s", filename);
   // the size of the file has to agree with the counts (a damaged file must not drive allocations)
   long long expected = 1 + 26 * 4;
   for(int c = 0; c < 26; c++) expected += 8LL * snpCounts[c];
   FSEEK64(fp, 0, SEEK_END);
   long long size = FTELL64(fp);
   if(size != expected)
      error("The size of %s (%lld bytes) does not agree with its SNP counts (%lld bytes expected)",
         filename, size, expected);
   FSEEK64(fp, 1 + 26 * 4, SEEK_SET);

   std::vector<unsigned int> pos;
   std::vector<float> wt;
   for(int c = 0; c < 26; c++){
      const size_t count = snpCounts[c];
      if(count == 0) continue;
      pos.resize(count);
      wt.resize(count);
      if(fread(&pos[0], sizeof(unsigned int), count, fp) < count)
         error("Not enough bytes for positions/alleles at the %dth chromosome in %s", c + 1, filename);
      if(fread(&wt[0], sizeof(float), count, fp) < count)
         error("Not enough bytes for weights at the %dth chromosome in %s", c + 1, filename);
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
      }
   }
   fclose(fp);
}

// ---------------------------------------------------------------------------------------
// --makebin
// ---------------------------------------------------------------------------------------

struct BinSNP
{
   unsigned int item;            // position << 4 | allele code
   float weight;
};

static bool ByPosition(const BinSNP & a, const BinSNP & b) { return (a.item >> 4) < (b.item >> 4); }

static void MakeOneBinaryModel(const std::string & textFile, const std::string & outFile)
{
   std::vector<BinSNP> snps[26];
   TextModelStats stats;
   int unsupportedAlleles = 0, positionTooLarge = 0;
   ParseTextModel(textFile, stats,
      [&](int chr, long pos, const std::string & ea, const std::string & oa, double wt){
         const int e = LetterCode(ea), o = LetterCode(oa);       // both alleles exist and differ
         if(e < 0 || o < 0 || e == o || !(fabs(wt) < 3.0e38)){unsupportedAlleles++; return;}
         if(pos > MAX_BIN_POSITION){positionTooLarge++; return;}
         BinSNP s;
         s.item = ((unsigned int) pos << 4) | (unsigned int)(e << 2 | o);
         s.weight = (float) wt;
         snps[chr - 1].push_back(s);
      });
   unsigned int counts[26];
   long long written = 0;
   for(int c = 0; c < 26; c++){
      std::stable_sort(snps[c].begin(), snps[c].end(), ByPosition);
      counts[c] = (unsigned int) snps[c].size();
      written += counts[c];
   }
   FILE *fp = fopen(outFile.c_str(), "wb");
   if(fp == NULL) error("Cannot open %s to write", outFile.c_str());
   fputc(BIN_VERSION, fp);
   fwrite(counts, sizeof(unsigned int), 26, fp);
   std::vector<unsigned int> item;
   std::vector<float> weight;
   for(int c = 0; c < 26; c++){
      if(counts[c] == 0) continue;
      item.resize(counts[c]);
      weight.resize(counts[c]);
      for(size_t m = 0; m < counts[c]; m++){
         item[m] = snps[c][m].item;
         weight[m] = snps[c][m].weight;
      }
      fwrite(&item[0], sizeof(unsigned int), counts[c], fp);
      fwrite(&weight[0], sizeof(float), counts[c], fp);
   }
   fclose(fp);
   printf("%s: %d SNPs read, %lld written to %s", textFile.c_str(), stats.fileSNPs, written, outFile.c_str());
   std::string notWritten;
   char part[160];
   if(stats.noCoord){snprintf(part, sizeof(part), "%s%d without hm_chr/hm_pos", notWritten.empty()? "": ", ", stats.noCoord); notWritten += part;}
   if(stats.noOtherAllele){snprintf(part, sizeof(part), "%s%d without other allele", notWritten.empty()? "": ", ", stats.noOtherAllele); notWritten += part;}
   if(stats.sameAllele){snprintf(part, sizeof(part), "%s%d with other allele = effect allele", notWritten.empty()? "": ", ", stats.sameAllele); notWritten += part;}
   if(unsupportedAlleles){snprintf(part, sizeof(part), "%s%d with alleles other than single A/C/G/T", notWritten.empty()? "": ", ", unsupportedAlleles); notWritten += part;}
   if(positionTooLarge){snprintf(part, sizeof(part), "%s%d with position above %ld", notWritten.empty()? "": ", ", positionTooLarge, MAX_BIN_POSITION); notWritten += part;}
   if(stats.malformed){snprintf(part, sizeof(part), "%s%d lines without effect_allele/effect_weight", notWritten.empty()? "": ", ", stats.malformed); notWritten += part;}
   if(!notWritten.empty()) printf(" (not written: %s)", notWritten.c_str());
   printf("\n");
   if(written == 0) printf("  Warning: %s has no SNPs.\n", outFile.c_str());
}

void MakeBinaryModels(const std::vector<std::string> & files, const std::vector<std::string> & names)
{
   const size_t n = files.size();
   std::vector<std::string> outFiles(n);
   for(size_t i = 0; i < n; i++){
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
   int made = 0;
   for(size_t i = 0; i < n; i++){
      if(IsBinaryModelName(files[i])){
         printf("%s is a binary model already and is skipped.\n", files[i].c_str());
         continue;
      }
      MakeOneBinaryModel(files[i], outFiles[i]);
      made++;
   }
   printf("%d binary model%s made.\n", made, made == 1? "": "s");
}
