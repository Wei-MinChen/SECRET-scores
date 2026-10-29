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
// frequency is used: a missing genotype contributes 0.
//
// Several models can be scored in one run.  They are processed in batches whose
// genotypes fit in GENOTYPE_MEMORY_LIMIT (half of the RAM available when the
// genotypes are about to be loaded), so the .bed file is read once for all
// the SNPs of a batch.  The genotypes are stored individual-major (see ReadPLINK.cpp),
// and the people are scored in parallel.
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <string>
#include <vector>
#include <utility>
#include <unordered_set>
#include <time.h>
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

// Parse a PGS Catalog text model without genotype data (see ScoreEngine.h).
void ParseTextModel(const std::string & file, TextModelStats & stats, const SNPCallback & onSNP)
{
   LineReader reader;
   FILE *test = fopen(file.c_str(), "rb");
   if(test == NULL) error("Cannot open %s to read", file.c_str());
   fclose(test);
   if(!reader.Open(file.c_str())) error("Cannot open %s to read", file.c_str());
   stats.fileSNPs = stats.malformed = stats.noCoord = stats.noOtherAllele = stats.sameAllele = 0;
   std::vector<std::string> f;
   int colEA = -1, colWT = -1, colCHR = -1, colPOS = -1, colOA = -1, colHMOA = -1;
   bool haveHeader = false, tabs = true;
   size_t needCols = 0;
   const char *s;
   size_t length;
   while(reader.Next(s, length)){
      if(s[0] == '\0' || s[0] == '#') continue;
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
         }
         if(colEA == -1 || colWT == -1 || colCHR == -1 || colPOS == -1 || (colOA == -1 && colHMOA == -1))
            error("The header of the model file %s needs the columns\n"
               "  effect_allele, effect_weight, hm_chr, hm_pos\n"
               "  and other_allele or hm_inferOtherAllele (PGS Catalog harmonized scoring file)", file.c_str());
         needCols = (size_t) std::max(std::max(colEA, colWT), std::max(colCHR, colPOS)) + 1;
         if(colOA + 1 > (int) needCols) needCols = (size_t) colOA + 1;
         if(colHMOA + 1 > (int) needCols) needCols = (size_t) colHMOA + 1;
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
      onSNP(chr, pos, ea, oa, wt);
   }
   if(!haveHeader) error("The model file %s has no header line", file.c_str());
}

// Read one text model file and match its SNPs to the genotype data (chromosome, position, alleles).
void ScoreEngine::ReadModelFile(ScoreModel & model)
{
   model.rowIdx.clear();
   model.weights.clear();
   model.switched.clear();
   model.nonAutosomal = model.notInData = model.alleleMismatch = 0;
   ScoreModel *md = &model;
   TextModelStats stats;
   ParseTextModel(model.file, stats,
      [this, md](int chr, long pos, const std::string & ea, const std::string & oa, double wt){
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
      });
   model.fileSNPs = stats.fileSNPs;
   model.malformed = stats.malformed;
   model.noCoord = stats.noCoord;
   model.noOtherAllele = stats.noOtherAllele;
   model.sameAllele = stats.sameAllele;
}

// Tabulate for each matched SNP what every genotype code contributes to the score: the effect
// weight w times the copies of the effect allele (2, 1 or 0), and 0 for a missing genotype.
// The effect allele is A1 (copies of A1 are 2, 1, 0 for the codes 0, 2, 3) or A2 (copies of A2
// are 0, 1, 2 for the codes 0, 2, 3; SWITCHED).  The weight is never changed.
void ScoreEngine::PrepareModel(ScoreModel & model, bool verbose)
{
   const int snpcount = (int) model.rowIdx.size();
   model.vslot.clear();
   model.vval.clear();
   model.vslot.reserve(snpcount);
   model.vval.reserve(4 * (size_t) snpcount);
   const bool table = verbose && snpcount <= SNP_TABLE_LIMIT;
   if(table){
      printf("\nThe %d SNPs of the model that are in the genotype data\n", snpcount);
      printf("%5s %11s %5s %5s %5s %12s %9s\n", "CHR", "POS", "A1", "A2", "EA", "WEIGHT", "ACTION");
   }
   for(int m = 0; m < snpcount; m++){
      const int row = model.rowIdx[m];
      const bool sw = model.switched[m] != 0;
      const double w = model.weights[m];
      if(table){
         char a1[8], a2[8];
         const unsigned int c1 = bimA1[row], c2 = bimA2[row];
         snprintf(a1, sizeof(a1), "%s", c1 < 256? std::string(1, (char) c1).c_str(): "<m>");
         snprintf(a2, sizeof(a2), "%s", c2 < 256? std::string(1, (char) c2).c_str(): "<m>");
         printf("%5d %11d %5s %5s %5s %12.6g %9s\n", (int) bimChr[row], bimPos[row], a1, a2,
            sw? a2: a1, w, sw? "SWITCHED": "OK");
      }
      model.vslot.push_back(slotOfRow[row]);
      model.vval.push_back(sw? 0.0: 2 * w);      // code 0: two copies of A1
      model.vval.push_back(0.0);                 // code 1: missing
      model.vval.push_back(w);                   // code 2: one copy of A1
      model.vval.push_back(sw? 2 * w: 0.0);      // code 3: no copy of A1
   }
   model.used = (int) model.vslot.size();
   model.scored = (model.used > 0);
}

// Score the prepared models first..last-1.  People are independent, so they are split among
// the threads; for each person all the models are scored while the person's genotypes are in
// cache, and each person adds up the SNPs of a model in the order of the model file.
void ScoreEngine::ScoreBatch(ScoreModel *models, int first, int last,
   std::vector<std::vector<double> > & scores, std::vector<double> * missSNP)
{
   // The people are done in 100 chunks, so that progress can be reported (and flushed to the
   // screen or the log file) every 30 seconds with an estimate of the time left.
   const int chunkSize = idCount > 100? (idCount + 99) / 100: 1;
   const time_t start = time(NULL);
   time_t lastReport = start;
   for(int c0 = 0; c0 < idCount; c0 += chunkSize){
      const int c1 = (c0 + chunkSize < idCount)? c0 + chunkSize: idCount;
#ifdef _OPENMP
      #pragma omp parallel for num_threads(threadCount) schedule(static) if(threadCount > 1)
#endif
      for(int i = c0; i < c1; i++){
         const unsigned char *gi = G_IND[i];
         for(int k = first; k < last; k++){
            ScoreModel & md = models[k];
            if(!md.scored) continue;
            const int V = md.used;
            const int *slot = &md.vslot[0];
            const double *val = &md.vval[0];
            double sum = 0.0;
            int ms = 0;
            for(int v = 0; v < V; v++){
               const int s = slot[v];
               const int g = (gi[s >> 2] >> ((s & 3) * 2)) & 3;
               sum += val[4*v + g];
               ms += (g == 1);
            }
            scores[k][i] = sum;
            if(missSNP && k == first)
               (*missSNP)[i] = (double) ms / md.used;
         }
      }
      const time_t now = time(NULL);
      if(c1 < idCount && now - lastReport >= 30){
         const double elapsed = difftime(now, start);
         const double done = (double) c1 / idCount;
         printf("  Scoring: %d of %d individuals (%.0f%%) done, %.0f s elapsed, about %.0f s left\n",
            c1, idCount, 100.0 * done, elapsed, elapsed * (1 - done) / done);
         fflush(stdout);
         lastReport = now;
      }
   }
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
   const int modelCount = (int) modelFiles.size();
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

   // ---- names of the scores ----
   // names have no white space (the score file is tab-delimited, and SECRET uses them as column names)
   ScoreModel *models = new ScoreModel[modelCount];
   std::vector<std::string> usedNames;
   for(int k = 0; k < modelCount; k++){
      models[k].file = modelFiles[k];
      models[k].binary = IsBinaryModelName(modelFiles[k]);
      std::string nm = modelNames[k].empty()? ModelNameFromFile(modelFiles[k]): modelNames[k];
      for(size_t i = 0; i < nm.size(); i++)
         if(nm[i] == ' ' || nm[i] == '\t' || nm[i] == ',' || nm[i] == '"') nm[i] = '_';
      if(nm.empty()) nm = "SCORE";
      for(int j = 0; j < k; j++)
         if(models[j].name == nm)
            error("Scores %s and %s have the same name %s.\n  Give them different names in the second column of the --readlist file.",
               models[j].file.c_str(), modelFiles[k].c_str(), nm.c_str());
      models[k].name = nm;
      models[k].scored = false;
      models[k].used = 0;
   }

   // ---- the models ----
   const bool single = (modelCount == 1);
   if(!single){printf("Read in %d models...\n", modelCount); fflush(stdout);}
   const time_t readStart = time(NULL);
   time_t lastReadReport = readStart;
   for(int k = 0; k < modelCount; k++){
      ScoreModel & md = models[k];
      if(single){printf("Read in %s%s...\n", md.file.c_str(), md.binary? " (binary model)": ""); fflush(stdout);}
      if(!single && k > 0 && time(NULL) - lastReadReport >= 30){
         lastReadReport = time(NULL);
         printf("  %d of %d models read (%.0f s)\n", k, modelCount, difftime(lastReadReport, readStart));
         fflush(stdout);
      }
      if(md.binary) ReadBinaryModelFile(md);
      else ReadModelFile(md);
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
            printf("  %d SNPs are not on an autosome and are not scored\n", md.nonAutosomal);
         if(md.notInData)
            printf("  %d SNPs are not in the genotype data (chromosome and position)\n", md.notInData);
         if(md.alleleMismatch)
            printf("  %d SNPs are in the genotype data but the alleles do not match\n", md.alleleMismatch);
         if(md.rowIdx.size() == 0) error("No SNPs were extracted");
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
   int threads = 1;
#ifdef _OPENMP
   if(CoreCount > 0)
      threads = CoreCount;
   else{
      int availableCoreCount = omp_get_max_threads();
      threads = availableCoreCount == 1? 1: availableCoreCount / 2;
   }
#endif
   printf("\nOptions in effect:\n");
   if(!optionReadModel.empty())
      printf("\t--readmodel %s\n", optionReadModel.c_str());
   if(!optionReadList.empty())
      printf("\t--readlist %s\n", optionReadList.c_str());
   if(CoreCount > 0)
      printf("\t--cpus %d\n", CoreCount);
   if(prefix != "secretscores")
      printf("\t--prefix %s\n", prefix.c_str());
   printf("\n%d model%s, %d thread%s\n\n", modelCount, modelCount == 1? "": "s", threads, threads == 1? "": "s");
   fflush(stdout);

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
      // models are added to the batch while the union of their SNPs fits in memory
      // (idCount people x ceil(SNPs/4) bytes)
      std::vector<int> touched;
      int batchEnd = batchStart;
      while(batchEnd < modelCount){
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
      }
      for(size_t t = 0; t < touched.size(); t++) needed[touched[t]] = 0;
      batchStart = batchEnd;
   }
   fclose(fp);

   int scoredCount = 0;
   for(int k = 0; k < modelCount; k++)
      if(models[k].scored) scoredCount++;
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
   fp = fopen(coveragefile.c_str(), "wt");
   if(fp == NULL) error("Cannot open %s to write", coveragefile.c_str());
   fprintf(fp, "SCORE\tFILE\tSNPs\tFound\tUsed\tcoverage\tNoCoord\tNoOtherAllele\tSameAllele\tNonAutosomal\tNotInData\tAlleleMismatch\n");
   for(int k = 0; k < modelCount; k++){
      const int found = (int) models[k].rowIdx.size() + models[k].alleleMismatch;
      char coverage[32];
      if(models[k].fileSNPs > 0) snprintf(coverage, sizeof(coverage), "%.4f", (double) models[k].used / models[k].fileSNPs);
      else snprintf(coverage, sizeof(coverage), "NA");
      fprintf(fp, "%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", models[k].name.c_str(), models[k].file.c_str(),
         models[k].fileSNPs, found, models[k].used, coverage, models[k].noCoord, models[k].noOtherAllele,
         models[k].sameAllele, models[k].nonAutosomal, models[k].notInData, models[k].alleleMismatch);
   }
   fclose(fp);

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
