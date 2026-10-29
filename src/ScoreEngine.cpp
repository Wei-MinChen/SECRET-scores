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
#include <algorithm>
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
   prefix = "secretscores";
   CoreCount = 0;
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
   return true;
}

void LineReader::Close()
{
   if(gz) gzclose(gz);
   gz = NULL;
}

bool LineReader::Next(const char *& line, size_t & length)
{
   if(gz == NULL) return false;
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
         int got = gzread(gz, &buffer[end], (unsigned int)(buffer.size() - end - 1));
         if(got <= 0) finished = true;
         else end += got;
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
      // a .fam line has 6 fields; only FID and IID are used (the order of the people is the file's)
      if(Tokens(line, tk, tl, 6) < 6) continue;
      FID.push_back(std::string(tk[0], tl[0]));
      IID.push_back(std::string(tk[1], tl[1]));
   }
   idCount = (int) IID.size();
   printf("  PLINK fam loaded: %d samples\n", idCount);
}

// Read every row of the .bim file.  Row r of the .bim file is row r of the .bed
// file, whatever its chromosome, so all rows are kept and the SNPs that cannot
// be scored (not on chromosomes 1-22 or XY) are only counted.
void ScoreEngine::ReadPlinkBim(const std::string & filename)
{
   bimChr.clear();
   bimPos.clear();
   bimA1.clear();
   bimA2.clear();
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
      counts[chr]++;
      long long key = ((long long) chr << 32) | (long long) pos;
      if(key < previousKey) bimSorted = false;
      previousKey = key;
   }
   reader.Close();
   bimCount = (int) bimPos.size();
   long long autosomes = 0;
   for(int c = 1; c < SEXCHR + 4; c++)
      if(ScoredChromosome(c)) autosomes += counts[c];
   unrecognized = counts[0];
   printf("  Genotype data consist of %lld autosome SNPs", autosomes);
   if(counts[SEXCHR+2]) printf(" (including %lld XY SNPs)", counts[SEXCHR+2]);
   if(counts[SEXCHR]) printf(", %lld X-chromosome SNPs", counts[SEXCHR]);
   if(counts[SEXCHR+1]) printf(", %lld Y-chromosome SNPs", counts[SEXCHR+1]);
   if(counts[SEXCHR+3]) printf(", %lld mitochondrial SNPs", counts[SEXCHR+3]);
   if(unrecognized) printf(", %lld SNPs on unrecognized chromosomes", unrecognized);
   printf("\n");
   if(autosomes < bimCount)
      printf("  %lld other SNPs cannot be scored.\n", (long long) bimCount - autosomes);
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
