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
// Reading of the PLINK binary genotypes (.bed) into individual-major storage.  The transposition of
// SNP-major genotypes (ReadSNPMajor) is shared with the VCF and .pgen readers (ReadVCF.cpp, ReadPGEN.cpp).
//
// A .bed file starts with 0x6c 0x1b and a mode byte:
//   0x00  individual-major: for each person in .fam order, the genotypes at all the
//         SNPs of the .bim file, 2 bits each, 4 SNPs per byte (ceil(M/4) bytes
//         per person)
//   0x01  SNP-major: for each SNP in .bim order, the genotypes of all the
//         people, 2 bits each, 4 people per byte (ceil(N/4) bytes per SNP)
// N is the number of samples of the .fam file (famCount).  With --keep / --remove the file is
// still read as a whole, but only the genotypes of the samples that are kept (idCount of them,
// sampleRow / famToSample map between .fam rows and kept samples) are stored.
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
#include "PgenReader.h"
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
   // the size of the .bed file must agree with the .fam and .bim files (all the samples of the
   // .fam file, whether they are kept by --keep / --remove or not)
   FSEEK64(fp, 0, SEEK_END);
   long long bedsize = FTELL64(fp);
   long long expected = bedMode == 1?
      3 + (long long)((famCount+3)/4) * (long long) bimCount:
      3 + (long long)famCount * (long long)((bimCount+3)/4);
   if(bedsize != expected)
      error("The size of %s (%lld bytes) does not match %d samples in %s.fam and %d SNPs in %s.bim (%lld bytes expected in %s-major mode)",
         bedfile.c_str(), bedsize, famCount, plinkinput.c_str(), bimCount,
         plinkinput.c_str(), expected, bedMode == 1? "SNP": "individual");
   printf("  PLINK binary genotypes are in %s-major mode\n", bedMode == 1? "SNP": "individual");
   if(bedMode == 1){
      printf("  Note: individual-major format is faster for this program, because each\n"
             "  person's genotypes are read in one piece. It is recommended to save the genotypes in\n"
             "  individual-major format once, with PLINK 2:\n"
             "     plink2 --bfile %s --export ind-major-bed --out %s_indmajor\n"
             "  and to use %s_indmajor with --bfile.\n", plinkinput.c_str(), plinkinput.c_str(), plinkinput.c_str());
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

// Individual-major file -> individual-major storage.  The file is read once, in large sequential
// blocks of whole persons (readBlockBytes, 256 MB unless --blocksize says otherwise); the SNPs of
// `rows` are picked out of each person's genotypes and repacked, 4 per byte.  Bits that are not
// SNPs of `rows` (the unused bits at the end of a person's bytes) are never looked at.
//
// With a whole-genome .bim a person is tens of MB (100 M SNPs: 25 MB), so a block holds only a
// few persons, or one.  Two things keep all the threads busy and the disk busy:
//  - the work is split into tasks of (person, SNP_CHUNK = 65536 SNPs), not into persons, so a block of two
//    persons still gives every thread something to do;
//  - two blocks are held: while the other threads unpack block b, thread 0 reads block b+1, and
//    joins in when its read is done.  The reading is serial (one stream of the file), the
//    unpacking is in parallel.
// Every task writes whole bytes of one person's storage that no other task writes (SNP_CHUNK is a
// multiple of 4), so the threads need no synchronization within a block.
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
   if(S == 0 || idCount == 0) return;                   // nothing to pick out of the file

   const long long fileBytes = ((long long) bimCount + 3) / 4;      // bytes of one person in the file
   const long long blockBytes = readBlockMB > 0? (long long) readBlockMB << 20: readBlockBytes;
   long long blockPersons = blockBytes / fileBytes;
   if(blockPersons < 1) blockPersons = 1;
   if(blockPersons > famCount) blockPersons = famCount;
   const int blocks = (int) ((famCount + blockPersons - 1) / blockPersons);      // blocks of the whole file
#ifndef SNP_CHUNK_SIZE
   #define SNP_CHUNK_SIZE (1 << 16)                     // SNPs per task, a multiple of 4
#endif
   const int SNP_CHUNK = SNP_CHUNK_SIZE;
   const int chunks = (S + SNP_CHUNK - 1) / SNP_CHUNK;

   if(blocks > 1){
      char size[64];
      snprintf(size, sizeof(size), "%.0f MB", (double) blockPersons * fileBytes / (1024.0 * 1024.0));
      printf("  Reading the .bed file in %d blocks of %lld persons (%s each)\n", blocks, blockPersons, size);
      fflush(stdout);
   }
   unsigned char *buf[2] = {NULL, NULL};
   buf[0] = new unsigned char[(size_t)(blockPersons * fileBytes)];
   if(blocks > 1) buf[1] = new unsigned char[(size_t)(blockPersons * fileBytes)];
   const int *r = &rows[0];

   FSEEK64(fp, 3, SEEK_SET);
   // read block b into buf[b & 1]; false when the file is too short (returns the first person it lacks)
   int shortAt = -1;
   #define READ_BLOCK(b) \
      do{ \
         const int f0 = (int) ((b) * blockPersons); \
         const int np = (int) ((famCount - f0 < blockPersons)? famCount - f0: blockPersons); \
         const size_t want = (size_t)(np * fileBytes); \
         const size_t got = fread(buf[(b) & 1], sizeof(unsigned char), want, fp); \
         if(got < want) shortAt = f0 + (int)(got / fileBytes); \
      }while(0)

   READ_BLOCK(0);
   if(shortAt < 0){
      int nextTask = 0;                                 // the next task of the block being unpacked
#ifdef _OPENMP
      #pragma omp parallel num_threads(threadCount) if(threadCount > 1)
#endif
      {
         int tid = 0;
#ifdef _OPENMP
         tid = omp_get_thread_num();
#endif
         for(int b = 0; b < blocks; b++){
            const unsigned char *cur = buf[b & 1];
            const int first = (int) (b * blockPersons);
            const int n = (int) ((famCount - first < blockPersons)? famCount - first: blockPersons);
            const int tasks = n * chunks;
            if(tid == 0 && b + 1 < blocks) READ_BLOCK(b + 1);   // the next block, while the others unpack this one
            for(;;){
               int t;
#ifdef _OPENMP
               #pragma omp atomic capture
#endif
               t = nextTask++;
               if(t >= tasks) break;
               const int p = t / chunks;
               const int who = famToSample[first + p];          // the kept sample, -1 when --keep / --remove drop this person
               if(who < 0) continue;
               const int k0 = (t % chunks) * SNP_CHUNK;
               const int k1 = (k0 + SNP_CHUNK < S)? k0 + SNP_CHUNK: S;
               const unsigned char *person = cur + (long long) p * fileBytes;
               unsigned char *out = G_IND[who];
               int k = k0;
               for(; k + 4 <= k1; k += 4){
                  const int r0 = r[k], r1 = r[k+1], r2 = r[k+2], r3 = r[k+3];
                  out[k >> 2] = (unsigned char)(
                     ((person[r0 >> 2] >> ((r0 & 3) * 2)) & 3) |
                     (((person[r1 >> 2] >> ((r1 & 3) * 2)) & 3) << 2) |
                     (((person[r2 >> 2] >> ((r2 & 3) * 2)) & 3) << 4) |
                     (((person[r3 >> 2] >> ((r3 & 3) * 2)) & 3) << 6));
               }
               if(k < k1){                               // the last 1-3 SNPs, in one more byte
                  unsigned char v = 0;
                  for(int j = 0; k + j < k1; j++){
                     const int row = r[k + j];
                     v |= (unsigned char)(((person[row >> 2] >> ((row & 3) * 2)) & 3) << (2 * j));
                  }
                  out[k >> 2] = v;
               }
            }
            // everybody is done with block b, and block b+1 is in memory
#ifdef _OPENMP
            #pragma omp barrier
#endif
            if(tid == 0) nextTask = 0;
#ifdef _OPENMP
            #pragma omp barrier
#endif
            if(shortAt >= 0) break;                      // the file is too short: the same for every thread
         }
      }
   }
   #undef READ_BLOCK
   delete []buf[0];
   delete []buf[1];
   if(shortAt >= 0)
      error("Not enough genotypes for sample number %d of %s.fam", (shortAt < famCount? shortAt: famCount - 1) + 1, plinkinput.c_str());
}

// SNP-major file -> individual-major storage (a transposition).  The SNPs of `rows` are read in
// blocks; each block is turned into whole bytes (4 SNPs) of every person's genotypes.
void ScoreEngine::ReadPLINKBinarySNPMajor(FILE *fp, std::vector<int> & rows)
{
   const long long SNP_space = (famCount + 3) / 4;      // bytes of one SNP in the file (all the .fam samples)
   FSEEK64(fp, 3, SEEK_SET);
   int previousRow = -1;
   ReadSNPMajor(rows, [&](int k0, int nb, unsigned char *buf){
      for(int j = 0; j < nb; j++){
         int row = rows[k0 + j];
         if(row - previousRow > 1)
            FSEEK64(fp, SNP_space * (long long)(row - previousRow - 1), SEEK_CUR);
         if(fread(buf + (long long) j * SNP_space, sizeof(unsigned char), (size_t) SNP_space, fp) < (size_t) SNP_space)
            error("Not enough genotypes at the %dth SNP (chromosome %d, position %d)\n", row + 1, (int) bimChr[row], bimPos[row]);
         previousRow = row;
      }
   });
}

// The transposition of a SNP-major source: fetch(k0, nb, buf) puts the genotypes of the SNPs rows[k0..k0+nb-1]
// in buf, one SNP after the other, (famCount+3)/4 bytes each (PLINK 1 coding of the genotype of every sample
// of the sample file, 4 per byte: 00 = 2 copies of A1, 01 = missing, 10 = one copy, 11 = no copy); the bytes
// are turned into the individual-major storage of the samples that are used.  A .bed file, a VCF and a .pgen
// file differ in how fetch() makes its block, not in what is done with it.
void ScoreEngine::ReadSNPMajor(std::vector<int> & rows, const FetchSNPs & fetch)
{
   const int S = (int) rows.size();
   loadedCount = S;
   genoBytes = (S + 3) / 4;
   G_IND = new unsigned char *[idCount > 0? idCount: 1];
   for(int i = 0; i < idCount; i++){
      G_IND[i] = new unsigned char[genoBytes > 0? genoBytes: 1];
      memset(G_IND[i], 0, genoBytes > 0? genoBytes: 1);
   }
   const long long SNP_space = (famCount + 3) / 4;      // bytes of one SNP in the block (all the samples of the file)
   // the block is 128 SNPs, fewer when a SNP is very long (a biobank: 4 million people = 1 MB per SNP)
   int BLOCK = 128;
   while(BLOCK > 4 && (long long) BLOCK * SNP_space > (256LL << 20)) BLOCK /= 2;      // a multiple of 4
   unsigned char *buf = new unsigned char[(size_t)(BLOCK * SNP_space)];
   for(int k0 = 0; k0 < S; k0 += BLOCK){
      int nb = S - k0 < BLOCK? S - k0: BLOCK;
      fetch(k0, nb, buf);
#ifdef _OPENMP
      #pragma omp parallel for num_threads(threadCount) schedule(static) if(threadCount > 1)
#endif
      for(int i = 0; i < idCount; i++){
         const int famRow = sampleRow[i];                 // the row of sample i in the sample file
         const int byteIdx = famRow >> 2, shift = (famRow & 3) * 2;
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

// ---- the genotype data of any format: .bed (--bfile), VCF (--vcf) or .pgen (--pfile) ----

// the samples and the variants; sets plinkinput (the name of the data in messages) and inputName
void ScoreEngine::LoadSamplesAndVariants(const std::string & name)
{
   std::string base = name;
   if(genotypeFormat == GENO_BED){
      if(base.size() > 4 && base.substr(base.size() - 4) == ".bed") base = base.substr(0, base.size() - 4);
      plinkinput = inputName = base;
      ReadPlinkFam(base);
      ReadPlinkBim(base);
   }else if(genotypeFormat == GENO_PGEN){
      const char *extensions[3] = {".pgen", ".pvar", ".psam"};
      for(int e = 0; e < 3; e++){
         const size_t n = strlen(extensions[e]);
         if(base.size() > n && base.substr(base.size() - n) == extensions[e]){base = base.substr(0, base.size() - n); break;}
      }
      plinkinput = inputName = base;
      ReadPsam(base);
      ReadPvar(base);
   }else{
      plinkinput = inputName = base;
      ReadVCFSamplesAndVariants(base);
   }
}

void ScoreEngine::OpenGenotypes()
{
   if(genotypeFormat == GENO_BED) bedFile = OpenPLINKBinary();
   else if(genotypeFormat == GENO_PGEN) OpenPgen(inputName);
}

void ScoreEngine::CloseGenotypes()
{
   if(bedFile) fclose(bedFile);
   bedFile = NULL;
   delete pgenReader;
   pgenReader = NULL;
}

void ScoreEngine::ReadGenotypes(std::vector<int> & rows)
{
   if(genotypeFormat == GENO_BED) ReadPLINKBinary(bedFile, rows);
   else if(genotypeFormat == GENO_PGEN) ReadPgenGenotypes(rows);
   else ReadVCFGenotypes(rows);
}
