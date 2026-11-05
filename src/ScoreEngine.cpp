////////////////////////////////////////////////////////////////////////
// ScoreEngine.cpp
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
// Reading of the PLINK .fam and .bim files, in large blocks (fread-style),
// so that files of biobank size (millions of people, 100 million SNPs) load
// quickly and the .bim file is held compactly in memory.
////////////////////////////////////////////////////////////////////////

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <limits>
#include <algorithm>
#include <math.h>
#include <unordered_map>
#include "ScoreEngine.h"
#include "Error.h"

ScoreEngine::ScoreEngine()
{
   idCount = 0;
   bimCount = 0;
   bimSorted = true;
   G_IND = NULL;
   genoBytes = 0;
   bedMode = -1;
   threadCount = 1;
   loadedCount = 0;
   maleHetX = 0;
   haveReadFreq = false;
   readBlockBytes = 256LL << 20;
   prefix = "secretscores";
   CoreCount = 0;
   readBlockMB = 0;
   maxBatchModels = MAX_MODELS_PER_BATCH;
}

ScoreEngine::~ScoreEngine()
{
   if(G_IND){
      for(int i = 0; i < idCount; i++)
         delete []G_IND[i];
      delete []G_IND;
   }
}

char *ScoreEngine::currentTime()
{
   time_t timer;
   struct tm *tblock;
   timer = time(NULL);
   tblock = localtime(&timer);
   return(asctime(tblock));
}

// ---------------------------------------------------------------------------------------
// LineReader: gzread() reads plain files as well as gzipped ones, in blocks of 4 MB; the
// lines are found with memchr in the block, with no per-character calls.
// ---------------------------------------------------------------------------------------
LineReader::LineReader()
{
   gz = NULL;
   begin = end = 0;
   finished = false;
   bad = false;
}

LineReader::~LineReader()
{
   Close();
}

bool LineReader::Open(const char *filename)
{
   Close();
   gz = gzopen(filename, "rb");
   if(gz == NULL) return false;
   gzbuffer(gz, 1 << 20);
   buffer.assign(4 << 20, 0);
   begin = end = 0;
   finished = false;
   bad = false;
   badMessage.clear();
   return true;
}

void LineReader::Close()
{
   if(gz) gzclose(gz);
   gz = NULL;
}

bool LineReader::Next(const char *& line, size_t & length)
{
   if(gz == NULL || bad) return false;
   for(;;){
      char *nl = (end > begin)? (char *) memchr(&buffer[begin], '\n', end - begin): NULL;
      if(nl == NULL && !finished){
         // keep the unfinished line, read more after it (the block grows for a very long line)
         if(begin > 0){
            memmove(&buffer[0], &buffer[begin], end - begin);
            end -= begin;
            begin = 0;
         }
         if(end + 1 >= buffer.size()) buffer.resize(buffer.size() * 2);
         const unsigned int wanted = (unsigned int)(buffer.size() - end - 1);
         int got = gzread(gz, &buffer[end], wanted);
         if(got <= 0) finished = true;
         else end += got;
         if(got < (int) wanted){
            // a short read is the end of the file only when zlib says so: a truncated or damaged
            // .gz file must not pass for a shorter file
            int errnum = 0;
            const char *message = gzerror(gz, &errnum);
            if(errnum != Z_OK && errnum != Z_STREAM_END){
               bad = true;
               badMessage = (message && message[0])? message: "read error";
               return false;
            }
         }
         continue;
      }
      if(nl == NULL){                       // the last line has no line break
         if(end == begin) return false;
         nl = &buffer[end];
      }
      *nl = '\0';
      char *start = &buffer[begin];
      size_t n = nl - start;
      if(n > 0 && start[n - 1] == '\r') start[--n] = '\0';
      begin = (nl - &buffer[0]) + 1;
      if(begin > end) begin = end;
      line = start;
      length = n;
      return true;
   }
}

// ---------------------------------------------------------------------------------------

int ChromosomeCode(const char *label, size_t length)
{
   if(length > 3 && (label[0]=='c' || label[0]=='C') && (label[1]=='h' || label[1]=='H')
      && (label[2]=='r' || label[2]=='R')){
      label += 3;
      length -= 3;
   }
   if(length == 0) return 0;
   if(length == 1 && (label[0]=='X' || label[0]=='x')) return SEXCHR;
   if(length == 1 && (label[0]=='Y' || label[0]=='y')) return SEXCHR+1;
   if(length == 2 && (label[0]=='X' || label[0]=='x') && (label[1]=='Y' || label[1]=='y')) return SEXCHR+2;
   if(length == 2 && (label[0]=='M' || label[0]=='m') && (label[1]=='T' || label[1]=='t')) return SEXCHR+3;
   int chr = 0;
   for(size_t i = 0; i < length; i++){
      if(label[i] < '0' || label[i] > '9') return 0;
      chr = chr * 10 + (label[i] - '0');
      if(chr > SEXCHR+3) return 0;
   }
   return (chr < 1 || chr > SEXCHR+3)? 0: chr;
}

unsigned int AlleleCode(const char *label, size_t length)
{
   if(length == 0) return NO_ALLELE;
   if(length == 1)
      return (label[0] >= 'a' && label[0] <= 'z')? (unsigned int)(label[0] - 'a' + 'A'): (unsigned int)(unsigned char) label[0];
   unsigned int h = 2166136261u;                 // FNV-1a, on the upper case letters
   for(size_t i = 0; i < length; i++){
      unsigned char c = (unsigned char) label[i];
      if(c >= 'a' && c <= 'z') c = (unsigned char)(c - 'a' + 'A');
      h = (h ^ c) * 16777619u;
   }
   return 0x80000000u | (h % 0x7fffffffu);
}

// ---------------------------------------------------------------------------------------

static bool IsSpace(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

// the first n white-space separated tokens of a line; returns how many there are (up to n)
static int Tokens(const char *line, const char **start, size_t *length, int n)
{
   int count = 0;
   const char *p = line;
   while(count < n){
      while(*p && IsSpace(*p)) p++;
      if(*p == '\0') break;
      const char *q = p;
      while(*q && !IsSpace(*q)) q++;
      start[count] = p;
      length[count] = q - p;
      count++;
      p = q;
   }
   return count;
}

static void SplitFieldsPublic(const char *s, std::vector<std::string> & fields)
{
   fields.clear();
   const char *p = s;
   while(*p){
      while(IsSpace(*p)) p++;
      const char *q = p;
      while(*q && !IsSpace(*q)) q++;
      if(q > p) fields.push_back(std::string(p, q - p));
      p = q;
   }
}

// open name+extension, or name+extension+".gz" when there is no such file
static bool OpenText(LineReader & reader, const std::string & name, const char *extension, std::string & opened)
{
   opened = name + extension;
   FILE *test = fopen(opened.c_str(), "rb");
   if(test == NULL){
      opened += ".gz";
      test = fopen(opened.c_str(), "rb");
      if(test == NULL) return false;
   }
   fclose(test);
   return reader.Open(opened.c_str());
}

void ScoreEngine::ReadPlinkFam(const std::string & filename)
{
   FID.clear();
   IID.clear();
   sex.clear();
   LineReader reader;
   std::string opened;
   if(!OpenText(reader, filename, ".fam", opened))
      error("Pedigree file %s.fam cannot be opened", filename.c_str());
   printf("Read in PLINK fam file %s...\n", opened.c_str());
   const char *line;
   size_t length;
   const char *tk[6];
   size_t tl[6];
   while(reader.Next(line, length)){
      // a .fam line has 6 fields; FID, IID and the sex (field 5) are used (the order of the
      // people is the file's).  Sex: 1 male, 2 female, anything else unknown.
      if(Tokens(line, tk, tl, 6) < 6) continue;
      FID.push_back(std::string(tk[0], tl[0]));
      IID.push_back(std::string(tk[1], tl[1]));
      sex.push_back((tl[4] == 1 && tk[4][0] == '1')? 1: (tl[4] == 1 && tk[4][0] == '2')? 2: 0);
   }
   {
      std::string why;
      if(reader.Failed(why)) error("Cannot read %s completely (%s): the file may be truncated or damaged", opened.c_str(), why.c_str());
   }
   idCount = (int) IID.size();
   printf("  PLINK fam loaded: %d samples\n", idCount);
}

// Read every row of the .bim file.  Row r of the .bim file is row r of the .bed
// file, whatever its chromosome, so all rows are kept and the SNPs that cannot
// be scored (not on chromosomes 1-22, X or XY) are only counted.
// 64-bit FNV-1a hash of a SNP ID (case-sensitive); the IDs of a .bim file are matched to those of
// a --readfreq file through it, so that 100 million IDs need 8 bytes each instead of strings.
static unsigned long long IdHash(const char *label, size_t length)
{
   unsigned long long h = 1469598103934665603ULL;
   for(size_t i = 0; i < length; i++)
      h = (h ^ (unsigned char) label[i]) * 1099511628211ULL;
   return h;
}

void ScoreEngine::ReadPlinkBim(const std::string & filename)
{
   bimChr.clear();
   bimPos.clear();
   bimA1.clear();
   bimA2.clear();
   bimIdHash.clear();
   const bool keepIds = !readFreqFile.empty() && MeanImputation;     // IDs are only needed to match --readfreq
   LineReader reader;
   std::string opened;
   if(!OpenText(reader, filename, ".bim", opened))
      error("Map file %s.bim cannot be opened", filename.c_str());
   printf("Read in PLINK bim file %s...\n", opened.c_str());
   // a plain file: count its lines first (fast), so that the arrays are allocated once
   if(opened.size() < 3 || opened.substr(opened.size() - 3) != ".gz"){
      FILE *fp = fopen(opened.c_str(), "rb");
      if(fp){
         std::vector<char> block(8 << 20);
         long long lines = 0;
         size_t got, last = '\n';
         while((got = fread(&block[0], 1, block.size(), fp)) > 0){
            const char *p = &block[0], *e = p + got;
            while((p = (const char *) memchr(p, '\n', e - p)) != NULL){lines++; p++;}
            last = (unsigned char) block[got - 1];
         }
         if(last != '\n') lines++;
         fclose(fp);
         if(lines > 0 && lines < INT_MAX){
            bimChr.reserve((size_t) lines);
            bimPos.reserve((size_t) lines);
            bimA1.reserve((size_t) lines);
            bimA2.reserve((size_t) lines);
            if(keepIds) bimIdHash.reserve((size_t) lines);
         }
      }
   }
   const char *line;
   size_t length;
   const char *tk[6];
   size_t tl[6];
   long long counts[SEXCHR + 4];
   for(int c = 0; c < SEXCHR + 4; c++) counts[c] = 0;
   long long unrecognized = 0;
   bimSorted = true;
   long long previousKey = -1;
   while(reader.Next(line, length)){
      if(Tokens(line, tk, tl, 6) < 6) continue;
      int chr = ChromosomeCode(tk[0], tl[0]);
      long pos = 0;
      {
         // the position: digits only, 0 when there are none
         const char *p = tk[3], *e = tk[3] + tl[3];
         long long v = 0;
         bool ok = (p < e);
         for(; p < e; p++){
            if(*p < '0' || *p > '9'){ok = false; break;}
            v = v * 10 + (*p - '0');
            if(v > INT_MAX){ok = false; break;}
         }
         pos = ok? (long) v: 0;
      }
      bimChr.push_back((unsigned char) chr);
      bimPos.push_back((int) pos);
      bimA1.push_back(AlleleCode(tk[4], tl[4]));
      bimA2.push_back(AlleleCode(tk[5], tl[5]));
      if(keepIds) bimIdHash.push_back(IdHash(tk[1], tl[1]));
      counts[chr]++;
      long long key = ((long long) chr << 32) | (long long) pos;
      if(key < previousKey) bimSorted = false;
      previousKey = key;
   }
   {
      std::string why;
      if(reader.Failed(why)) error("Cannot read %s completely (%s): the file may be truncated or damaged", opened.c_str(), why.c_str());
   }
   reader.Close();
   bimCount = (int) bimPos.size();
   // autosomes and XY; X is counted apart, and is scored unless --allowX 0
   long long autosomes = counts[SEXCHR + 2];
   for(int c = 1; c < SEXCHR; c++) autosomes += counts[c];
   const long long scorable = autosomes + (ScoredChromosome(SEXCHR)? counts[SEXCHR]: 0);
   unrecognized = counts[0];
   printf("  Genotype data consist of %lld autosome SNPs", autosomes);
   if(counts[SEXCHR+2]) printf(" (including %lld XY SNPs)", counts[SEXCHR+2]);
   if(counts[SEXCHR]) printf(", %lld X-chromosome SNPs%s", counts[SEXCHR], ScoredChromosome(SEXCHR)? "": " (not scored)");
   if(counts[SEXCHR+1]) printf(", %lld Y-chromosome SNPs", counts[SEXCHR+1]);
   if(counts[SEXCHR+3]) printf(", %lld mitochondrial SNPs", counts[SEXCHR+3]);
   if(unrecognized) printf(", %lld SNPs on unrecognized chromosomes", unrecognized);
   printf("\n");
   if(scorable < bimCount)
      printf("  %lld other SNPs cannot be scored.\n", (long long) bimCount - scorable);
   // SNPs that are not in (chromosome, position) order (rare) are looked up through a sorted index
   posIndex.clear();
   if(!bimSorted){
      printf("  The .bim file is not sorted by chromosome and position; building an index...\n");
      for(int r = 0; r < bimCount; r++)
         if(ScoredChromosome(bimChr[r]) && bimPos[r] > 0)
            posIndex.push_back(std::make_pair(((long long) bimChr[r] << 32) | (long long) bimPos[r], r));
      std::sort(posIndex.begin(), posIndex.end());
   }
   printf("  PLINK maps loaded: %d SNPs\n", bimCount);
}

// --readfreq.  The first line (a leading # is dropped) names the columns:
//   SNP or ID                    the SNP ID of the .bim file
//   A1 or ALT                    the allele whose frequency is given
//   MAF or ALT_FREQS             that frequency (PLINK .frq: A1 and MAF; PLINK 2 .afreq: ALT and ALT_FREQS)
// A line with several ALT alleles (a comma in the frequency), a missing or invalid frequency, or the
// ID "." is skipped.  The allele has to be A1 or A2 of the .bim row; if it is A2, A1 gets 1 - frequency.
void ScoreEngine::ReadFreqFile(const std::string & filename)
{
   haveReadFreq = false;
   bimFreqA1.assign((size_t) bimCount, std::numeric_limits<float>::quiet_NaN());
   if(bimIdHash.size() != (size_t) bimCount)
      error("The SNP IDs of the .bim file were not kept for --readfreq");
   LineReader reader;
   FILE *test = fopen(filename.c_str(), "rb");
   if(test == NULL) error("Cannot open %s to read (--readfreq)", filename.c_str());
   fclose(test);
   if(!reader.Open(filename.c_str())) error("Cannot open %s to read (--readfreq)", filename.c_str());
   printf("Read in allele frequencies %s...\n", filename.c_str());

   // the .bim IDs, sorted by hash, to find the rows of an ID
   std::vector<std::pair<unsigned long long, int> > idIndex;
   idIndex.reserve((size_t) bimCount);
   for(int r = 0; r < bimCount; r++) idIndex.push_back(std::make_pair(bimIdHash[r], r));
   std::sort(idIndex.begin(), idIndex.end());

   const char *line;
   size_t length;
   int colId = -1, colAllele = -1, colFreq = -1;
   bool haveHeader = false;
   long long lines = 0, used = 0, notInData = 0, alleleMismatch = 0, invalid = 0;
   std::vector<std::string> f;
   while(reader.Next(line, length)){
      const char *p = line;
      while(*p == ' ' || *p == '\t') p++;
      if(*p == '\0') continue;
      if(!haveHeader){
         while(*p == '#') p++;
         SplitFieldsPublic(p, f);
         for(size_t i = 0; i < f.size(); i++){
            std::string u = f[i];
            for(size_t j = 0; j < u.size(); j++) if(u[j] >= 'a' && u[j] <= 'z') u[j] = (char)(u[j] - 'a' + 'A');
            if(u == "SNP" || u == "ID") colId = (int) i;
            else if(u == "A1" || u == "ALT") colAllele = (int) i;
            else if(u == "MAF" || u == "ALT_FREQS") colFreq = (int) i;
         }
         if(colId < 0 || colAllele < 0 || colFreq < 0)
            error("The header of the --readfreq file %s needs the columns\n"
               "  SNP, A1, MAF (PLINK --freq, .frq) or ID, ALT, ALT_FREQS (PLINK 2 --freq, .afreq)", filename.c_str());
         haveHeader = true;
         continue;
      }
      SplitFieldsPublic(p, f);
      const size_t need = (size_t) std::max(colId, std::max(colAllele, colFreq)) + 1;
      if(f.size() < need) continue;
      lines++;
      const std::string & id = f[colId];
      if(id == ".") continue;
      char *end = NULL;
      const double fr = strtod(f[colFreq].c_str(), &end);
      if(end == f[colFreq].c_str() || *end != '\0' || !(fr >= 0.0 && fr <= 1.0)){invalid++; continue;}
      const unsigned int allele = AlleleCode(f[colAllele].c_str(), f[colAllele].size());
      const unsigned long long h = IdHash(id.c_str(), id.size());
      std::vector<std::pair<unsigned long long, int> >::const_iterator it =
         std::lower_bound(idIndex.begin(), idIndex.end(), std::make_pair(h, -1));
      bool found = false, matched = false;
      for(; it != idIndex.end() && it->first == h; ++it){
         found = true;
         const int r = it->second;
         if(allele == bimA1[r]){bimFreqA1[r] = (float) fr; matched = true;}
         else if(allele == bimA2[r]){bimFreqA1[r] = (float)(1.0 - fr); matched = true;}
      }
      if(!found) notInData++;
      else if(!matched) alleleMismatch++;
      else used++;
   }
   {
      std::string why;
      if(reader.Failed(why)) error("Cannot read %s completely (%s): the file may be truncated or damaged", filename.c_str(), why.c_str());
   }
   reader.Close();
   std::vector<unsigned long long>().swap(bimIdHash);
   if(!haveHeader) error("The --readfreq file %s is empty", filename.c_str());
   printf("  %lld frequencies read, %lld of them match SNPs of the genotype data", lines, used);
   if(notInData) printf("; %lld SNP IDs are not in the .bim file", notInData);
   if(alleleMismatch) printf("; %lld have alleles that do not match the .bim file", alleleMismatch);
   if(invalid) printf("; %lld lines skipped (missing, multi-allelic or invalid frequency)", invalid);
   printf("\n");
   if(used == 0) error("None of the frequencies in %s match the SNPs of the genotype data (by SNP ID and allele)", filename.c_str());
   haveReadFreq = true;
}

int ScoreEngine::MatchSNP(int chr, long pos, unsigned int ea, unsigned int oa, bool & sw, bool & atPosition) const
{
   atPosition = false;
   const long long key = ((long long) chr << 32) | (long long) pos;
   int hit = -1;
   sw = false;
   if(bimSorted){
      // binary search for the first row at the position, in the .bim arrays themselves
      size_t lo = 0, hi = (size_t) bimCount;
      while(lo < hi){
         size_t mid = lo + (hi - lo) / 2;
         long long k = ((long long) bimChr[mid] << 32) | (long long) bimPos[mid];
         if(k < key) lo = mid + 1;
         else hi = mid;
      }
      for(size_t r = lo; r < (size_t) bimCount; r++){
         if((((long long) bimChr[r] << 32) | (long long) bimPos[r]) != key) break;
         atPosition = true;
         const unsigned int a1 = bimA1[r], a2 = bimA2[r];
         if(ea == a1 && oa == a2){hit = (int) r; sw = false;}
         else if(ea == a2 && oa == a1){hit = (int) r; sw = true;}
         if(hit > -1) break;
      }
   }else{
      std::vector<std::pair<long long, int> >::const_iterator it =
         std::lower_bound(posIndex.begin(), posIndex.end(), std::make_pair(key, -1));
      for(; it != posIndex.end() && it->first == key; ++it){
         atPosition = true;
         const int r = it->second;
         const unsigned int a1 = bimA1[r], a2 = bimA2[r];
         if(ea == a1 && oa == a2){hit = r; sw = false;}
         else if(ea == a2 && oa == a1){hit = r; sw = true;}
         if(hit > -1) break;
      }
   }
   return hit;
}
