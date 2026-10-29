////////////////////////////////////////////////////////////////////////
// ReadPLINK.cpp
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
// Reading of the PLINK binary genotypes (.bed) into individual-major storage.
//
// A .bed file starts with 0x6c 0x1b and a mode byte:
//   0x00  individual-major: for each person in .fam order, the genotypes at all the
//         SNPs of the .bim file, 2 bits each, 4 SNPs per byte (ceil(M/4) bytes
//         per person)
//   0x01  SNP-major: for each SNP in .bim order, the genotypes of all the
//         people, 2 bits each, 4 people per byte (ceil(N/4) bytes per SNP)
// The 2-bit codes are the same: 00 = 2 copies of allele 1 (A1), 01 = missing,
// 10 = one copy of A1, 11 = no copy of A1.
//
// Whatever the file mode, the genotypes of the SNPs being scored are kept in
// memory as G_IND[i], person i's genotypes at those SNPs, so that every
// analysis runs person by person.
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include "ScoreEngine.h"
#include "Error.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

// Open the .bed file, check its header and size, and set bedMode.
FILE *ScoreEngine::OpenPLINKBinary()
{
   std::string bedfile = plinkinput + ".bed";
   FILE *fp = fopen(bedfile.c_str(), "rb");
   if(fp == NULL) error("Cannot open PLINK bed file %s", bedfile.c_str());
   printf("Read in PLINK bed file %s...\n", bedfile.c_str());
   int b0 = fgetc(fp);
   if(b0 != 0x6c)
      error("The PLINK format in %s cannot be recognized. The first byte is %x",
         bedfile.c_str(), b0);
   int b1 = fgetc(fp);
   if(b1 != 0x1b)
      error("The PLINK format in %s cannot be recognized. The second byte is %x",
         bedfile.c_str(), b1);
   int mode = fgetc(fp);
   if(mode != 0 && mode != 1)
      error("The storage mode of %s cannot be recognized. The third byte is %x\n"
         "  (0 = individual-major, 1 = SNP-major)", bedfile.c_str(), mode);
   bedMode = mode;
   // the size of the .bed file must agree with the .fam and .bim files
   FSEEK64(fp, 0, SEEK_END);
   long long bedsize = FTELL64(fp);
   long long expected = bedMode == 1?
      3 + (long long)((idCount+3)/4) * (long long) bimCount:
      3 + (long long)idCount * (long long)((bimCount+3)/4);
   if(bedsize != expected)
      error("The size of %s (%lld bytes) does not match %d samples in %s.fam and %d SNPs in %s.bim (%lld bytes expected in %s-major mode)",
         bedfile.c_str(), bedsize, idCount, plinkinput.c_str(), bimCount,
         plinkinput.c_str(), expected, bedMode == 1? "SNP": "individual");
   printf("  PLINK binary genotypes are in %s-major mode\n", bedMode == 1? "SNP": "individual");
   if(bedMode == 1){
      printf("  Note: individual-major format is faster for this program, because each\n"
             "  person's genotypes are read in one piece. It is recommended to save the genotypes in\n"
             "  individual-major format once, with PLINK 2:\n"
             "     plink2 --bfile %s --export ind-major-bed --out %s_indmajor\n"
             "  and to use %s_indmajor with -b.\n", plinkinput.c_str(), plinkinput.c_str(), plinkinput.c_str());
   }
   return fp;
}

void ScoreEngine::ReadPLINKBinary(FILE *fp, std::vector<int> & rows)
{
   if(bedMode == 0)
      ReadPLINKBinaryIndividualMajor(fp, rows);
   else
      ReadPLINKBinarySNPMajor(fp, rows);
}

// Individual-major file -> individual-major storage.  The file is read once, person by person
// (large sequential reads); the SNPs of `rows` are picked out of each person's genotypes
// and repacked, 4 per byte.  Bits that are not SNPs of `rows` (the unused bits at the
// end of a person's bytes) are never looked at.
void ScoreEngine::ReadPLINKBinaryIndividualMajor(FILE *fp, std::vector<int> & rows)
{
   const int S = (int) rows.size();
   loadedCount = S;
   genoBytes = (S + 3) / 4;
   G_IND = new unsigned char *[idCount > 0? idCount: 1];
   for(int i = 0; i < idCount; i++){
      G_IND[i] = new unsigned char[genoBytes > 0? genoBytes: 1];
      memset(G_IND[i], 0, genoBytes > 0? genoBytes: 1);
   }
   const long long fileBytes = (bimCount + 3) / 4;      // bytes of one person in the file
   long long blockPersons = (64LL * 1024 * 1024) / fileBytes;   // ~64 MB at a time
   if(blockPersons < 1) blockPersons = 1;
   if(blockPersons > idCount) blockPersons = idCount;
   unsigned char *buf = new unsigned char[(size_t)(blockPersons * fileBytes)];
   const int *r = &rows[0];
   FSEEK64(fp, 3, SEEK_SET);
   for(int first = 0; first < idCount; first += (int) blockPersons){
      int n = (int) blockPersons;
      if(first + n > idCount) n = idCount - first;
      if(fread(buf, sizeof(unsigned char), (size_t)(n * fileBytes), fp) < (size_t)(n * fileBytes))
         error("Not enough genotypes for the sample %s", IID[first].c_str());
#ifdef _OPENMP
      #pragma omp parallel for num_threads(threadCount) schedule(static) if(threadCount > 1 && n > 1)
#endif
      for(int p = 0; p < n; p++){
         const unsigned char *person = buf + (long long) p * fileBytes;
         unsigned char *out = G_IND[first + p];
         for(int k = 0; k < S; k++){
            int row = r[k];
            int code = (person[row >> 2] >> ((row & 3) * 2)) & 3;
            out[k >> 2] |= (unsigned char)(code << ((k & 3) * 2));
         }
      }
   }
   delete []buf;
}

// SNP-major file -> individual-major storage (a transposition).  The SNPs of `rows` are read in
// blocks; each block is turned into whole bytes (4 SNPs) of every person's genotypes.
void ScoreEngine::ReadPLINKBinarySNPMajor(FILE *fp, std::vector<int> & rows)
{
   const int S = (int) rows.size();
   loadedCount = S;
   genoBytes = (S + 3) / 4;
   G_IND = new unsigned char *[idCount > 0? idCount: 1];
   for(int i = 0; i < idCount; i++){
      G_IND[i] = new unsigned char[genoBytes > 0? genoBytes: 1];
      memset(G_IND[i], 0, genoBytes > 0? genoBytes: 1);
   }
   const long long SNP_space = (idCount + 3) / 4;       // bytes of one SNP in the file
   const int BLOCK = 128;                              // SNPs per block, a multiple of 4
   unsigned char *buf = new unsigned char[(size_t)(BLOCK * SNP_space)];
   FSEEK64(fp, 3, SEEK_SET);
   int previousRow = -1;
   for(int k0 = 0; k0 < S; k0 += BLOCK){
      int nb = S - k0 < BLOCK? S - k0: BLOCK;
      for(int j = 0; j < nb; j++){
         int row = rows[k0 + j];
         if(row - previousRow > 1)
            FSEEK64(fp, SNP_space * (long long)(row - previousRow - 1), SEEK_CUR);
         if(fread(buf + (long long) j * SNP_space, sizeof(unsigned char), (size_t) SNP_space, fp) < (size_t) SNP_space)
            error("Not enough genotypes at the %dth SNP (chromosome %d, position %d)\n", row + 1, (int) bimChr[row], bimPos[row]);
         previousRow = row;
      }
#ifdef _OPENMP
      #pragma omp parallel for num_threads(threadCount) schedule(static) if(threadCount > 1)
#endif
      for(int i = 0; i < idCount; i++){
         const int byteIdx = i >> 2, shift = (i & 3) * 2;
         unsigned char *out = G_IND[i] + (k0 >> 2);
         for(int g = 0; g < nb; g += 4){
            unsigned char v = 0;
            for(int j = 0; j < 4 && g + j < nb; j++){
               int code = (buf[(long long)(g + j) * SNP_space + byteIdx] >> shift) & 3;
               v |= (unsigned char)(code << (2 * j));
            }
            out[g >> 2] = v;
         }
      }
   }
   delete []buf;
}

void ScoreEngine::FreeGenotypes(std::vector<int> & rows)
{
   if(G_IND){
      for(int i = 0; i < idCount; i++)
         delete []G_IND[i];
      delete []G_IND;
   }
   G_IND = NULL;
   for(size_t m = 0; m < rows.size(); m++)
      slotOfRow[rows[m]] = -1;
   loadedCount = 0;
   genoBytes = 0;
}
