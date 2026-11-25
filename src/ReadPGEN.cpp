////////////////////////////////////////////////////////////////////////
// ReadPGEN.cpp
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
// --pfile <prefix>: the PLINK 2 fileset <prefix>.pgen, <prefix>.pvar and <prefix>.psam.
//
// .psam: the header line (starting with #) names the columns; FID (optional), IID and SEX are used.  Sex: 1 male,
//        2 female, anything else unknown.
// .pvar: lines starting with ## are skipped; the line starting with #CHROM names the columns; CHROM, POS, ID, REF
//        and ALT are used (without a #CHROM line the columns are in that order).  Plain or .gz (a .pvar.zst has to
//        be decompressed first).  A variant with several ALT alleles gives one SNP for each ALT allele (the allele
//        is A1, REF is A2); a variant whose ALT is "." gives none.
// .pgen: PgenReader.cpp.  Only hard calls are used; dosages are not read.
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ScoreEngine.h"
#include "PgenReader.h"
#include "Error.h"

// the tab-separated fields of a line (a .pvar / VCF line): pointers and lengths, up to n
static int TabFields(const char *line, const char **start, size_t *length, int n)
{
   int count = 0;
   const char *p = line;
   while(count < n){
      const char *q = p;
      while(*q && *q != '\t' && *q != '\r' && *q != '\n') q++;
      start[count] = p;
      length[count] = (size_t)(q - p);
      count++;
      if(*q != '\t') break;
      p = q + 1;
   }
   return count;
}

static bool NameIs(const char *p, size_t length, const char *name)
{
   size_t n = strlen(name);
   if(length > 0 && p[0] == '#'){p++; length--;}
   if(length != n) return false;
   for(size_t i = 0; i < n; i++){
      char c = p[i];
      if(c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
      if(c != name[i]) return false;
   }
   return true;
}

void ScoreEngine::ReadPsam(const std::string & prefix)
{
   LineReader reader;
   std::string opened;
   if(!OpenText(reader, prefix, ".psam", opened))
      error("Sample file %s.psam cannot be opened", prefix.c_str());
   printf("Read in PLINK 2 psam file %s...\n", opened.c_str());
   const char *line;
   size_t length;
   std::vector<std::string> fids, iids;
   std::vector<unsigned char> sexes;
   int colFid = -1, colIid = -1, colSex = -1;
   bool haveHeader = false;
   const char *tk[64];
   size_t tl[64];
   bool first = true;
   while(reader.Next(line, length)){
      if(first && strncmp(line, "\xEF\xBB\xBF", 3) == 0) line += 3;
      first = false;
      if(line[0] == '\0') continue;
      if(line[0] == '#'){
         if(line[1] == '#') continue;
         const int n = TabFields(line, tk, tl, 64);
         for(int c = 0; c < n; c++){
            if(NameIs(tk[c], tl[c], "FID")) colFid = c;
            else if(NameIs(tk[c], tl[c], "IID")) colIid = c;
            else if(NameIs(tk[c], tl[c], "SEX")) colSex = c;
         }
         if(colIid < 0) error("The header line of %s has no IID column", opened.c_str());
         haveHeader = true;
         continue;
      }
      if(!haveHeader)
         error("%s has no header line (a line that starts with #IID or #FID); a .psam file needs one", opened.c_str());
      const int n = TabFields(line, tk, tl, 64);
      if(n <= colIid) continue;
      fids.push_back(colFid >= 0 && colFid < n? std::string(tk[colFid], tl[colFid]): std::string("0"));
      iids.push_back(std::string(tk[colIid], tl[colIid]));
      unsigned char s = 0;
      if(colSex >= 0 && colSex < n && tl[colSex] == 1) s = (tk[colSex][0] == '1')? 1: (tk[colSex][0] == '2')? 2: 0;
      sexes.push_back(s);
   }
   {
      std::string why;
      if(reader.Failed(why)) error("Cannot read %s completely (%s): the file may be truncated or damaged", opened.c_str(), why.c_str());
   }
   reader.Close();
   if(!haveHeader) error("%s is empty", opened.c_str());
   sampleFileLabel = ".psam file";
   sampleHasFID = (colFid >= 0);
   sampleSexSource = (colSex >= 0)? ".psam": "";
   SelectSamples(fids, iids, sexes, "PLINK 2 psam");
}

void ScoreEngine::ReadPvar(const std::string & prefix)
{
   BeginVariantTable();
   rowRecord.clear();
   rowAlt.clear();
   recordAltCount.clear();
   std::string opened = prefix + ".pvar";
   {
      FILE *test = fopen(opened.c_str(), "rb");
      if(test == NULL){
         opened += ".gz";
         test = fopen(opened.c_str(), "rb");
      }
      if(test == NULL){
         if(FILE *z = fopen((prefix + ".pvar.zst").c_str(), "rb")){
            fclose(z);
            error("%s.pvar.zst is compressed with zstd, which is not read here.\n"
               "  Decompress it (plink2 --zst-decompress %s.pvar.zst > %s.pvar) or write a .pvar.gz file", prefix.c_str(), prefix.c_str(), prefix.c_str());
         }
         error("Variant file %s.pvar cannot be opened", prefix.c_str());
      }
      fclose(test);
   }
   LineReader reader;
   if(!reader.Open(opened.c_str())) error("Variant file %s cannot be opened", opened.c_str());
   printf("Read in PLINK 2 pvar file %s...\n", opened.c_str());
   const bool keepIds = wantVariantIds || (!readFreqFile.empty() && MeanImputation);
   const char *line;
   size_t length;
   int colChr = 0, colPos = 1, colId = 2, colRef = 3, colAlt = 4;
   const char *tk[256];
   size_t tl[256];
   bool first = true;
   long long records = 0;
   while(reader.Next(line, length)){
      if(first && strncmp(line, "\xEF\xBB\xBF", 3) == 0) line += 3;
      first = false;
      if(line[0] == '\0') continue;
      if(line[0] == '#'){
         if(line[1] == '#') continue;
         const int n = TabFields(line, tk, tl, 256);
         colChr = colPos = colId = colRef = colAlt = -1;
         for(int c = 0; c < n; c++){
            if(NameIs(tk[c], tl[c], "CHROM")) colChr = c;
            else if(NameIs(tk[c], tl[c], "POS")) colPos = c;
            else if(NameIs(tk[c], tl[c], "ID")) colId = c;
            else if(NameIs(tk[c], tl[c], "REF")) colRef = c;
            else if(NameIs(tk[c], tl[c], "ALT")) colAlt = c;
         }
         if(colChr < 0 || colPos < 0 || colRef < 0 || colAlt < 0)
            error("The header line of %s needs the columns CHROM, POS, REF and ALT", opened.c_str());
         continue;
      }
      const int n = TabFields(line, tk, tl, 256);
      const int need = std::max(std::max(colChr, colPos), std::max(std::max(colRef, colAlt), colId));
      if(n <= need) error("Line %lld of the variants of %s has fewer than %d columns", records + 1, opened.c_str(), need + 1);
      const int chr = ChromosomeCode(tk[colChr], tl[colChr]);
      const long pos = ParsePositionToken(tk[colPos], tl[colPos]);
      const unsigned int ref = AlleleCode(tk[colRef], tl[colRef]);
      const char *alts = tk[colAlt];
      const size_t altLength = tl[colAlt];
      int altCt = 0;
      if(!(altLength == 1 && alts[0] == '.')){
         size_t s = 0;
         while(s <= altLength){
            size_t e = s;
            while(e < altLength && alts[e] != ',') e++;
            altCt++;
            if(altCt > 254) error("Variant %lld of %s has more than 254 ALT alleles", records + 1, opened.c_str());
            rowRecord.push_back((int) records);
            rowAlt.push_back((unsigned char) altCt);
            AddVariant(chr, pos, AlleleCode(alts + s, e - s), ref);
            if(keepIds) bimIdHash.push_back(colId >= 0? IdHash(tk[colId], tl[colId]): 0);
            s = e + 1;
         }
      }
      recordAltCount.push_back((unsigned char)(altCt > 0? altCt: 1));
      records++;
      if(records >= 2147483647LL) error("Too many variants in %s", opened.c_str());
   }
   {
      std::string why;
      if(reader.Failed(why)) error("Cannot read %s completely (%s): the file may be truncated or damaged", opened.c_str(), why.c_str());
   }
   reader.Close();
   // rows are the first ALT allele of the variants of the same number: no need for the maps
   if(rowRecord.size() == (size_t) records){
      bool same = true;
      for(size_t r = 0; r < rowRecord.size() && same; r++)
         if(rowRecord[r] != (int) r || rowAlt[r] != 1) same = false;
      if(same){rowRecord.clear(); rowAlt.clear();}
   }
   pgenRecordCount = (int) records;
   FinishVariantTable(".pvar file");
   printf("  PLINK 2 variants loaded: %lld variants, %d SNPs\n", records, bimCount);
}

void ScoreEngine::OpenPgen(const std::string & prefix)
{
   delete pgenReader;
   pgenReader = new PgenReader();
   std::string err;
   const std::string name = prefix + ".pgen";
   printf("Read in PLINK 2 pgen file %s...\n", name.c_str());
   if(!pgenReader->Open(name, (uint32_t) famCount, (uint32_t) pgenRecordCount, err))
      error("%s", err.c_str());
   printf("  PLINK 2 hard calls (storage mode 0x%02x)\n", pgenReader->StorageMode());
}

void ScoreEngine::ReadPgenGenotypes(std::vector<int> & rows)
{
   const long long SNP_space = (famCount + 3) / 4;
   ReadSNPMajor(rows, [&](int k0, int nb, unsigned char *buf){
      std::string err;
      for(int j = 0; j < nb; j++){
         const int row = rows[k0 + j];
         const int rec = rowRecord.empty()? row: rowRecord[row];
         const unsigned alt = rowAlt.empty()? 1: rowAlt[row];
         if(!pgenReader->ReadRow((uint32_t) rec, recordAltCount[rec], alt, buf + (long long) j * SNP_space, err))
            error("%s (variant %d of the .pvar file, chromosome %d, position %d)", err.c_str(), rec + 1, (int) bimChr[row], bimPos[row]);
      }
   });
}
