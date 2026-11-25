////////////////////////////////////////////////////////////////////////
// PgenReader.h
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
// A reader of the hard calls of a PLINK 2 .pgen file, written from the
// public PGEN format specification
// (https://github.com/chrchang/plink-ng/tree/master/pgen_spec); it does not
// use the pgenlib library.  See PgenReader.cpp for what is supported.
////////////////////////////////////////////////////////////////////////

#ifndef __PGENREADER_H__
#define __PGENREADER_H__

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>

class PgenReader
{
   public:
      PgenReader();
      ~PgenReader();

      // Open a .pgen file.  sampleCt and variantCt are the numbers of samples (.psam) and variants (.pvar) that
      // the file has to agree with.  False, with the reason in err, when the file cannot be read.
      bool Open(const std::string & filename, uint32_t sampleCt, uint32_t variantCt, std::string & err);

      uint32_t SampleCount() const { return N; }
      uint32_t VariantCount() const { return M; }
      int StorageMode() const { return mode; }

      // The hard calls of ALT allele `alt` (1 = the first ALT allele) of variant v, which has altCt ALT alleles
      // (from the .pvar), for all the samples, in PLINK 1 coding with allele A1 = that ALT allele: 4 samples per
      // byte, (N+3)/4 bytes, 00 = 2 copies of the ALT allele, 01 = missing, 10 = one copy, 11 = no copy.
      // Variants are best asked for in increasing order.  False, with the reason in err, when the record is damaged.
      bool ReadRow(uint32_t v, unsigned altCt, unsigned alt, unsigned char *out, std::string & err);

   private:
      PgenReader(const PgenReader &);
      PgenReader & operator = (const PgenReader &);

      FILE *ff;                          // the .pgen file
      FILE *hf;                          // the file with the header: ff, or the .pgen.pgi index
      int mode;                          // the storage mode (third byte)
      uint32_t N, M;                     // samples, variants
      uint32_t sampleBytes;              // (N + 3) / 4
      uint64_t fileSize;
      // fixed-width modes (0x01 - 0x04)
      bool fixedWidth;
      uint64_t fixedStart, fixedRecord;
      // variable-width modes
      unsigned ctrl;                     // the twelfth byte
      unsigned vrtypeStorage;            // ctrl & 15
      unsigned altCountBytes;            // bytes of an allele count per variant in the header (0 or 1)
      bool nonrefFlags;                  // the header stores the provisional-REF flags
      uint64_t bodyBase;                 // where the body of the first variant block starts in hf
      uint64_t bodyFull;                 // bytes of the body of a block of 65536 variants
      std::vector<uint64_t> blockFpos;   // where the records of each variant block start in the .pgen file
      // the variant block that is loaded: the types, the lengths and the file positions of its records
      long long loadedBlock;
      std::vector<unsigned char> types, altCts;
      std::vector<uint32_t> lengths;
      std::vector<uint64_t> positions;
      // the latest variant that is not LD-compressed (the base of the LD-compressed ones): its main track
      long long baseVariant;
      std::vector<unsigned char> baseCats, baseCatsCopy, work, record;

      bool LoadBlock(uint32_t block, std::string & err);
      bool ReadRecord(uint32_t v, std::string & err);
      // the main track of variant v (not LD-compressed) in cats (packed 2-bit categories 0 = REF-REF, 1 = REF-ALT1,
      // 2 = ALT1-ALT1, 3 = missing); *end is the end of the main track in `record`
      bool DecodeMainTrack(unsigned type, const unsigned char *p, const unsigned char *e, const std::vector<unsigned char> *base,
         std::vector<unsigned char> & cats, const unsigned char **end, std::string & err);
      bool DecodeMultiallelic(const unsigned char *cats, const unsigned char *p, const unsigned char *e, unsigned altCt,
         unsigned alt, unsigned char *out, std::string & err);
};

#endif
