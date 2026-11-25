////////////////////////////////////////////////////////////////////////
// ReadVCF.cpp
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
// --vcf <file>: a VCF file (plain or gzipped, .vcf or .vcf.gz; a BCF file is not read).
//
// The header line #CHROM ... gives the samples (the columns after FORMAT).  The IID of a sample is its name
// and its FID is 0; VCF has no sex, so every sample is of unknown sex (the X chromosome is scored as for females).
// Every record gives one SNP for each ALT allele (ALT_j is A1, REF is A2); a record whose ALT is "." gives none.
// The genotype of a sample is its GT field (the field named GT in FORMAT), "0/1", "1|1" ...:
//    both alleles given               the copies of the ALT allele (a phased genotype is not treated differently)
//    one allele (haploid)             as homozygous, as PLINK does
//    a missing allele ("./1", "./.")  missing
//    more than two alleles, no GT     missing
// Only hard calls are used; dosages (DS, GP) are not read.
// The file is read from start to end for every batch of models (only the lines of the SNPs of the batch are parsed).
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ScoreEngine.h"
#include "Error.h"

static int TabFieldsVCF(const char *line, const char **start, size_t *length, int n)
{
   int count = 0;
   const char *p = line;
   while(count < n){
      const char *q = p;
      while(*q && *q != '\t') q++;
      start[count] = p;
      length[count] = (size_t)(q - p);
      count++;
      if(*q != '\t') break;
      p = q + 1;
   }
   return count;
}

static bool IsBCF(const std::string & name)
{
   if(name.size() >= 4 && name.substr(name.size() - 4) == ".bcf") return true;
   FILE *fp = fopen(name.c_str(), "rb");
   if(fp == NULL) return false;
   unsigned char b[3] = {0, 0, 0};
   size_t got = fread(b, 1, 3, fp);
   fclose(fp);
   return got == 3 && b[0] == 'B' && b[1] == 'C' && b[2] == 'F';
}

void ScoreEngine::ReadVCFSamplesAndVariants(const std::string & filename)
{
   if(IsBCF(filename))
      error("%s is a BCF file, which is not read here. Convert it to VCF (bcftools view -Oz -o out.vcf.gz %s) or to PLINK 2 "
         "(plink2 --bcf %s --make-pgen --out <prefix>, then use --pfile)", filename.c_str(), filename.c_str(), filename.c_str());
   BeginVariantTable();
   rowRecord.clear();
   rowAlt.clear();
   recordAltCount.clear();
   FILE *test = fopen(filename.c_str(), "rb");
   if(test == NULL) error("Cannot open VCF file %s to read", filename.c_str());
   fclose(test);
   LineReader reader;
   if(!reader.Open(filename.c_str())) error("Cannot open VCF file %s to read", filename.c_str());
   printf("Read in VCF file %s...\n", filename.c_str());
   const bool keepIds = wantVariantIds || (!readFreqFile.empty() && MeanImputation);
   const char *line;
   size_t length;
   const char *tk[10];
   size_t tl[10];
   bool haveHeader = false;
   long long records = 0, lineNumber = 0;
   std::vector<std::string> fids, iids;
   while(reader.Next(line, length)){
      lineNumber++;
      if(line[0] == '\0') continue;
      if(line[0] == '#'){
         if(line[1] == '#') continue;
         // the #CHROM line: the samples are the columns after the ninth
         const char *p = line;
         int column = 0;
         for(;;){
            const char *q = p;
            while(*q && *q != '\t') q++;
            if(column >= 9) iids.push_back(std::string(p, (size_t)(q - p)));
            column++;
            if(*q != '\t') break;
            p = q + 1;
         }
         if(column < 8) error("The #CHROM line of %s has fewer than 8 columns", filename.c_str());
         if(iids.empty())
            error("%s has no samples (no genotype columns after FORMAT)", filename.c_str());
         fids.assign(iids.size(), "0");
         {
            std::vector<std::string> sorted(iids);
            std::sort(sorted.begin(), sorted.end());
            for(size_t i = 1; i < sorted.size(); i++)
               if(sorted[i] == sorted[i - 1]){
                  printf("  Warning: sample ID %s occurs more than once in %s\n", sorted[i].c_str(), filename.c_str());
                  break;
               }
         }
         sampleFileLabel = "VCF header";
         sampleHasFID = false;
         sampleSexSource = "";
         SelectSamples(fids, iids, std::vector<unsigned char>(iids.size(), 0), "VCF");
         haveHeader = true;
         continue;
      }
      if(!haveHeader) error("%s has a variant line before its #CHROM line (line %lld)", filename.c_str(), lineNumber);
      const int n = TabFieldsVCF(line, tk, tl, 5);
      if(n < 5) error("Line %lld of %s has fewer than 5 columns", lineNumber, filename.c_str());
      const int chr = ChromosomeCode(tk[0], tl[0]);
      const long pos = ParsePositionToken(tk[1], tl[1]);
      const unsigned int ref = AlleleCode(tk[3], tl[3]);
      const char *alts = tk[4];
      const size_t altLength = tl[4];
      int altCt = 0;
      if(altLength > 0 && !(altLength == 1 && alts[0] == '.')){
         size_t s = 0;
         while(s <= altLength){
            size_t e = s;
            while(e < altLength && alts[e] != ',') e++;
            altCt++;
            if(altCt > 254) error("Line %lld of %s has more than 254 ALT alleles", lineNumber, filename.c_str());
            rowRecord.push_back((int) records);
            rowAlt.push_back((unsigned char) altCt);
            AddVariant(chr, pos, AlleleCode(alts + s, e - s), ref);
            if(keepIds) bimIdHash.push_back(IdHash(tk[2], tl[2]));
            s = e + 1;
         }
      }
      recordAltCount.push_back((unsigned char)(altCt > 0? altCt: 1));
      records++;
      if(records >= 2147483647LL) error("Too many records in %s", filename.c_str());
   }
   {
      std::string why;
      if(reader.Failed(why)) error("Cannot read %s completely (%s): the file may be truncated or damaged", filename.c_str(), why.c_str());
   }
   reader.Close();
   if(!haveHeader) error("%s has no #CHROM line: it is not a VCF file", filename.c_str());
   if(rowRecord.size() == (size_t) records){
      bool same = true;
      for(size_t r = 0; r < rowRecord.size() && same; r++)
         if(rowRecord[r] != (int) r || rowAlt[r] != 1) same = false;
      if(same){rowRecord.clear(); rowAlt.clear();}
   }
   pgenRecordCount = (int) records;
   FinishVariantTable("VCF file");
   printf("  VCF records loaded: %lld records, %d SNPs, %d samples\n", records, bimCount, famCount);
}

// The genotypes of the SNPs `rows` (increasing): the file is read from the start, and the lines of the records that
// have a SNP of `rows` are parsed.
void ScoreEngine::ReadVCFGenotypes(std::vector<int> & rows)
{
   LineReader reader;
   if(!reader.Open(inputName.c_str())) error("Cannot open VCF file %s to read", inputName.c_str());
   const long long SNP_space = (famCount + 3) / 4;
   const char *line = NULL;
   size_t length = 0;
   long long record = -1;               // the record of `line`
   long long parsed = -1;               // the record whose genotypes are in a1 / a2
   long long lineNumber = 0;
   bool pastHeader = false;
   std::vector<unsigned char> a1(famCount), a2(famCount);        // alleles of the samples: 255 = missing
   ReadSNPMajor(rows, [&](int k0, int nb, unsigned char *buf){
      for(int j = 0; j < nb; j++){
         const int row = rows[k0 + j];
         const long long rec = rowRecord.empty()? row: rowRecord[row];
         const unsigned alt = rowAlt.empty()? 1: rowAlt[row];
         while(record < rec){
            if(!reader.Next(line, length)) error("%s has fewer records than when it was first read", inputName.c_str());
            lineNumber++;
            if(line[0] == '\0') continue;
            if(line[0] == '#'){pastHeader = true; continue;}
            record++;
         }
         (void) pastHeader;
         if(parsed != rec){
            parsed = -1;
            // the columns: skip the first 8, the ninth is FORMAT
            const char *p = line;
            for(int c = 0; c < 8; c++){
               while(*p && *p != '\t') p++;
               if(*p != '\t') error("Line %lld of %s has fewer than 10 columns", lineNumber, inputName.c_str());
               p++;
            }
            const char *fmt = p;
            while(*p && *p != '\t') p++;
            const char *fmtEnd = p;
            // the index of GT in FORMAT
            int gtIndex = -1, idx = 0;
            for(const char *q = fmt; q <= fmtEnd; ){
               const char *e = q;
               while(e < fmtEnd && *e != ':') e++;
               if(e - q == 2 && q[0] == 'G' && q[1] == 'T'){gtIndex = idx; break;}
               idx++;
               q = e + 1;
            }
            const unsigned maxAllele = recordAltCount[rec];
            if(gtIndex < 0){
               memset(&a1[0], 255, a1.size());
               memset(&a2[0], 255, a2.size());
            }else{
               int sampleIndex = 0;
               while(*p == '\t'){
                  p++;
                  if(sampleIndex >= famCount) error("Line %lld of %s has more samples than its header", lineNumber, inputName.c_str());
                  // the GT subfield of the sample
                  const char *g = p;
                  bool absent = false;
                  for(int f = 0; f < gtIndex && !absent; f++){
                     while(*g && *g != ':' && *g != '\t') g++;
                     if(*g != ':') absent = true;
                     else g++;
                  }
                  unsigned al[3] = {255, 255, 255};
                  int alleles = 0;
                  bool missing = false;
                  for(; !absent;){
                     if(*g == '.'){missing = true; g++;}
                     else if(*g >= '0' && *g <= '9'){
                        unsigned v = 0;
                        while(*g >= '0' && *g <= '9'){v = v * 10 + (unsigned)(*g - '0'); if(v > 100000) v = 100000; g++;}
                        if(v > maxAllele)
                           error("Line %lld of %s: sample %d has the allele %u in GT but the record has only %u ALT allele%s",
                              lineNumber, inputName.c_str(), sampleIndex + 1, v, maxAllele, maxAllele == 1? "": "s");
                        if(alleles < 3) al[alleles] = v;
                        alleles++;
                     }else break;
                     if(*g == '/' || *g == '|') g++;
                     else break;
                  }
                  if(missing || alleles == 0 || alleles > 2){a1[sampleIndex] = a2[sampleIndex] = 255;}
                  else{
                     a1[sampleIndex] = (unsigned char) al[0];
                     a2[sampleIndex] = (unsigned char)(alleles == 1? al[0]: al[1]);
                  }
                  sampleIndex++;
                  while(*p && *p != '\t') p++;
               }
               if(sampleIndex != famCount)
                  error("Line %lld of %s has %d samples but its header has %d", lineNumber, inputName.c_str(), sampleIndex, famCount);
            }
            parsed = rec;
         }
         unsigned char *out = buf + (long long) j * SNP_space;
         memset(out, 0, (size_t) SNP_space);
         for(int i = 0; i < famCount; i++){
            unsigned code;
            if(a1[i] == 255) code = 1;
            else{
               const int copies = (a1[i] == alt) + (a2[i] == alt);
               code = (copies == 2)? 0: (copies == 1)? 2: 3;
            }
            out[i >> 2] |= (unsigned char)(code << ((i & 3) * 2));
         }
      }
   });
   reader.Close();
}
