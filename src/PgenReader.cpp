////////////////////////////////////////////////////////////////////////
// PgenReader.cpp
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
// A reader of the hard calls of a PLINK 2 .pgen file, written from the public PGEN format
// specification (pgen_spec in https://github.com/chrchang/plink-ng); it does not use pgenlib.
//
// Supported: the storage modes 0x01 (a PLINK 1 variant-major .bed), 0x02, 0x03 and 0x04 (fixed width), 0x10 and
// 0x11 (variable width) and 0x20 and 0x21 (variable width, with the index in the .pgen.pgi file); the variant
// record types of the main data track, 0 (raw 2-bit), 1 (1-bit with a difflist), 2 and 3 (LD-compressed, the
// second inverted), 4, 5, 6 and 7 (difflists from category 0, 2, 3, or empty); and the multiallelic hard-calls
// (bit 3 of the record type) with up to 254 ALT alleles.  The other tracks of a record (hard-call phase, dosages)
// are never read: only hard calls are used.  Not supported: a sample-major .bed (mode 0; PLINK 2 converts it).
//
// The main data track of a variant has 4 categories for every sample: 0 = homozygous REF, 1 = heterozygous
// REF-ALT1, 2 = homozygous ALT1, 3 = missing; ReadRow() turns it into the PLINK 1 coding of the copies of one
// ALT allele (00 = 2 copies, 01 = missing, 10 = 1 copy, 11 = none), 4 samples per byte.
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "PgenReader.h"

#ifdef _WIN32
  #define PG_FSEEK _fseeki64
#else
  #define PG_FSEEK fseeko
#endif

static const uint32_t VBLOCK = 65536;       // variants per variant block

static void Fail(std::string & err, const char *format, ...)
{
   char text[2048];
   va_list ap;
   va_start(ap, format);
   vsnprintf(text, sizeof(text), format, ap);
   va_end(ap);
   err = text;
}

// ---- the tables that turn 4 categories of a byte at once ----
struct PgenTables
{
   unsigned char toPlink[256];         // category -> copies of ALT1: 0 -> 11 (none), 1 -> 10, 2 -> 00, 3 -> 01 (missing)
   unsigned char toPlinkNone[256];     // an ALT allele after the first, without a patch: missing stays missing, the rest is "none"
   unsigned char invert[256];          // categories 0 and 2 swapped
   unsigned short oneBit[6][256];      // for the 6 pairs of categories of the 1-bit representation: 8 samples -> 16 bits
   PgenTables()
   {
      static const unsigned char plinkOf[4] = {3, 2, 0, 1};
      for(int b = 0; b < 256; b++){
         unsigned v = 0, w = 0, x = 0;
         for(int j = 0; j < 4; j++){
            const unsigned c = (b >> (2 * j)) & 3;
            v |= (unsigned) plinkOf[c] << (2 * j);
            w |= (unsigned)(c == 3? 1: 3) << (2 * j);
            x |= (unsigned)(c == 0? 2: c == 2? 0: c) << (2 * j);
         }
         toPlink[b] = (unsigned char) v;
         toPlinkNone[b] = (unsigned char) w;
         invert[b] = (unsigned char) x;
      }
      static const unsigned char lo[6] = {0, 0, 0, 1, 1, 2};
      static const unsigned char hi[6] = {1, 2, 3, 2, 3, 3};
      for(int k = 0; k < 6; k++)
         for(int b = 0; b < 256; b++){
            unsigned v = 0;
            for(int j = 0; j < 8; j++)
               v |= (unsigned)(((b >> j) & 1)? hi[k]: lo[k]) << (2 * j);
            oneBit[k][b] = (unsigned short) v;
         }
   }
};
static const PgenTables & Tables()
{
   static const PgenTables t;
   return t;
}

// the two most common categories of a 1-bit record: the code of the byte that starts it
static int OneBitPair(unsigned code)
{
   switch(code){
      case 1: return 0;     // 0, 1
      case 2: return 1;     // 0, 2
      case 3: return 2;     // 0, 3
      case 5: return 3;     // 1, 2
      case 6: return 4;     // 1, 3
      case 9: return 5;     // 2, 3
   }
   return -1;
}

static inline void SetCategory(unsigned char *cats, uint32_t i, unsigned c)
{
   const unsigned shift = (i & 3) * 2;
   cats[i >> 2] = (unsigned char)((cats[i >> 2] & ~(3u << shift)) | (c << shift));
}

static inline unsigned GetCategory(const unsigned char *cats, uint32_t i)
{
   return (cats[i >> 2] >> ((i & 3) * 2)) & 3;
}

// a base-128 varint
static bool ReadVarint(const unsigned char *& p, const unsigned char *e, uint32_t & v)
{
   v = 0;
   unsigned shift = 0;
   while(p < e){
      const unsigned char b = *p++;
      v |= (uint32_t)(b & 127) << shift;
      if(!(b & 128)) return true;
      shift += 7;
      if(shift > 28) return false;
   }
   return false;
}

// A difflist (see the specification): the increasing sample IDs, and the 2-bit categories with them when withCats.
static bool ParseDifflist(const unsigned char *& p, const unsigned char *e, uint32_t N, bool withCats,
   std::vector<uint32_t> & ids, std::vector<unsigned char> & cats, std::string & err)
{
   ids.clear();
   cats.clear();
   uint32_t L = 0;
   if(!ReadVarint(p, e, L)){Fail(err, "A difflist is damaged"); return false;}
   if(L == 0) return true;
   if(L > N){Fail(err, "A difflist has more samples (%u) than the file", L); return false;}
   const uint32_t G = (L + 63) / 64;
   // pgenlib uses 1 + floor(log2(N) / 8) bytes: 2 bytes already at N = 256 (the specification says N <= 256 for 1 byte)
   const unsigned w = (N < 256)? 1: (N < 65536)? 2: (N < 16777216)? 3: 4;
   const uint64_t need = (uint64_t) G * w + (G - 1) + (withCats? (L + 3) / 4: 0);
   if((uint64_t)(e - p) < need){Fail(err, "A difflist is cut short"); return false;}
   const unsigned char *firsts = p;
   p += (size_t) G * w;
   p += G - 1;                                     // the sizes of the groups: not needed to read the list in order
   const unsigned char *packed = p;
   if(withCats) p += (L + 3) / 4;
   ids.resize(L);
   for(uint32_t g = 0; g < G; g++){
      uint32_t first = 0;
      for(unsigned k = 0; k < w; k++) first |= (uint32_t) firsts[(size_t) g * w + k] << (8 * k);
      const uint32_t n = (L - g * 64 < 64)? L - g * 64: 64;
      ids[g * 64] = first;
      for(uint32_t k = 1; k < n; k++){
         uint32_t d;
         if(!ReadVarint(p, e, d)){Fail(err, "A difflist is cut short"); return false;}
         ids[g * 64 + k] = ids[g * 64 + k - 1] + d;
      }
   }
   for(uint32_t k = 0; k < L; k++)
      if(ids[k] >= N || (k > 0 && ids[k] <= ids[k - 1])){Fail(err, "A difflist has a sample ID out of order or out of range"); return false;}
   if(withCats){
      cats.resize(L);
      for(uint32_t k = 0; k < L; k++) cats[k] = (packed[k >> 2] >> (2 * (k & 3))) & 3;
   }
   return true;
}

static inline unsigned ReadBits(const unsigned char *p, uint64_t bitStart, unsigned width)
{
   unsigned v = 0;
   for(unsigned k = 0; k < width; k++){
      const uint64_t b = bitStart + k;
      v |= (unsigned)((p[b >> 3] >> (b & 7)) & 1) << k;
   }
   return v;
}

PgenReader::PgenReader()
{
   ff = hf = NULL;
   mode = 0;
   N = M = 0;
   sampleBytes = 0;
   fileSize = 0;
   fixedWidth = false;
   fixedStart = fixedRecord = 0;
   ctrl = vrtypeStorage = altCountBytes = 0;
   nonrefFlags = false;
   bodyBase = bodyFull = 0;
   loadedBlock = -1;
   baseVariant = -1;
}

PgenReader::~PgenReader()
{
   if(hf && hf != ff) fclose(hf);
   if(ff) fclose(ff);
}

// bytes of the part of the header that describes a variant block of n variants
static uint64_t BlockBodyBytes(unsigned vs, unsigned altCountBytes, bool nonrefFlags, uint64_t n)
{
   uint64_t bytes;
   if(vs >= 8){
      const unsigned bits = (vs == 8)? 1: (vs == 11 || vs == 14)? 4: 2;
      bytes = (n * bits + 7) / 8;
   }else{
      bytes = ((vs < 4)? (n + 1) / 2: n) + n * (1 + (vs & 3));
   }
   bytes += n * altCountBytes;
   if(nonrefFlags) bytes += (n + 7) / 8;
   return bytes;
}

bool PgenReader::Open(const std::string & filename, uint32_t sampleCt, uint32_t variantCt, std::string & err)
{
   ff = fopen(filename.c_str(), "rb");
   if(ff == NULL){Fail(err, "Cannot open %s to read", filename.c_str()); return false;}
   hf = ff;
   unsigned char head[3];
   if(fread(head, 1, 3, ff) < 3 || head[0] != 0x6c || head[1] != 0x1b){
      Fail(err, "%s is not a .pgen file: it does not start with the bytes 0x6c 0x1b", filename.c_str());
      return false;
   }
   mode = head[2];
   PG_FSEEK(ff, 0, SEEK_END);
   fileSize = (uint64_t) ftello(ff);
   PG_FSEEK(ff, 3, SEEK_SET);
   N = sampleCt;
   M = variantCt;
   sampleBytes = (N + 3) / 4;
   if(mode == 0){
      Fail(err, "%s is a sample-major PLINK 1 .bed file, not a .pgen file that can be read here.\n"
         "  Convert it with PLINK 2:  plink2 --bfile <prefix> --make-pgen --out <prefix2>   (or read it with --bfile)", filename.c_str());
      return false;
   }
   if(mode == 1){                                  // a PLINK 1 variant-major .bed: nothing else in the header
      fixedWidth = true;
      fixedStart = 3;
      fixedRecord = sampleBytes;
      if(fileSize != fixedStart + (uint64_t) M * fixedRecord){
         Fail(err, "The size of %s (%llu bytes) does not match %u samples and %u variants (%llu bytes expected)", filename.c_str(),
            (unsigned long long) fileSize, N, M, (unsigned long long)(fixedStart + (uint64_t) M * fixedRecord));
         return false;
      }
      return true;
   }
   if(mode != 2 && mode != 3 && mode != 4 && mode != 0x10 && mode != 0x11 && mode != 0x20 && mode != 0x21){
      Fail(err, "The storage mode of %s (the third byte, 0x%02x) is not known to this reader", filename.c_str(), mode);
      return false;
   }
   if(mode == 0x20 || mode == 0x21){               // the header is in the index file
      const std::string pgi = filename + ".pgi";
      hf = fopen(pgi.c_str(), "rb");
      if(hf == NULL){hf = ff; Fail(err, "%s needs its index file %s, which cannot be opened", filename.c_str(), pgi.c_str()); return false;}
      unsigned char h2[3];
      if(fread(h2, 1, 3, hf) < 3 || h2[0] != 0x6c || h2[1] != 0x1b || h2[2] != mode + 0x10){
         Fail(err, "%s is not the index file of %s", pgi.c_str(), filename.c_str());
         return false;
      }
   }
   uint32_t headM = 0, headN = 0;
   unsigned char c = 0;
   if(fread(&headM, 4, 1, hf) < 1 || fread(&headN, 4, 1, hf) < 1 || fread(&c, 1, 1, hf) < 1){
      Fail(err, "Not enough bytes for the header of %s", filename.c_str());
      return false;
   }
   if(headM != variantCt || headN != sampleCt){
      Fail(err, "%s has %u variants and %u samples, but the .pvar file has %u variants and the .psam file has %u samples",
         filename.c_str(), headM, headN, variantCt, sampleCt);
      return false;
   }
   ctrl = c;
   const unsigned nonrefStorage = ctrl >> 6;
   if(mode <= 4){
      // fixed width: every record has the same type
      if(ctrl & 63){Fail(err, "The header of %s is damaged (twelfth byte)", filename.c_str()); return false;}
      fixedWidth = true;
      fixedStart = 12 + (nonrefStorage == 3? ((uint64_t) M + 7) / 8: 0);
      fixedRecord = sampleBytes + (mode == 3? (uint64_t) 2 * N: mode == 4? (uint64_t) 4 * N: 0);
      if(fileSize != fixedStart + (uint64_t) M * fixedRecord){
         Fail(err, "The size of %s (%llu bytes) does not match %u samples and %u variants (%llu bytes expected)", filename.c_str(),
            (unsigned long long) fileSize, N, M, (unsigned long long)(fixedStart + (uint64_t) M * fixedRecord));
         return false;
      }
      return true;
   }
   vrtypeStorage = ctrl & 15;
   altCountBytes = (ctrl >> 4) & 3;
   nonrefFlags = (nonrefStorage == 3);
   if(altCountBytes > 1){Fail(err, "%s has more than 254 ALT alleles at a variant, which this reader does not support", filename.c_str()); return false;}
   if(vrtypeStorage >= 8){
      if(vrtypeStorage == 13 || vrtypeStorage == 15 || (vrtypeStorage & 7) == 5 || (vrtypeStorage & 7) == 7){
         Fail(err, "The header of %s (twelfth byte 0x%02x) is a storage layout that this reader does not support", filename.c_str(), ctrl);
         return false;
      }
      if(vrtypeStorage >= 12 && N != 1){Fail(err, "The header of %s describes one sample, not %u", filename.c_str(), N); return false;}
   }
   const uint32_t B = (M + VBLOCK - 1) / VBLOCK;
   blockFpos.resize(B);
   if(M == 0 || fread(&blockFpos[0], 8, B, hf) < B){
      Fail(err, "Not enough bytes for the variant block offsets of %s", filename.c_str());
      return false;
   }
   bodyBase = 12 + (uint64_t) 8 * B;
   bodyFull = BlockBodyBytes(vrtypeStorage, altCountBytes, nonrefFlags, VBLOCK);
   return true;
}

// Load the types, lengths and positions of the records of a variant block.
bool PgenReader::LoadBlock(uint32_t block, std::string & err)
{
   loadedBlock = -1;
   const uint32_t first = block * VBLOCK;
   const uint32_t n = (M - first < VBLOCK)? M - first: VBLOCK;
   const uint64_t bytes = BlockBodyBytes(vrtypeStorage, altCountBytes, nonrefFlags, n);
   std::vector<unsigned char> body((size_t) bytes);
   if(PG_FSEEK(hf, (long long)(bodyBase + (uint64_t) block * bodyFull), SEEK_SET) != 0
      || (bytes > 0 && fread(&body[0], 1, (size_t) bytes, hf) < (size_t) bytes)){
      Fail(err, "Not enough bytes for the header of variant block %u of the .pgen file", block);
      return false;
   }
   types.assign(n, 0);
   lengths.assign(n, 0);
   positions.assign(n, 0);
   const unsigned char *p = bodyBase == 0 || body.empty()? NULL: &body[0];
   if(vrtypeStorage >= 8){
      // fused record types and lengths in 1, 2 or 4 bits
      unsigned bw = (vrtypeStorage == 8)? 1: (vrtypeStorage == 11 || vrtypeStorage == 14)? 4: 2;
      unsigned char tt[16];
      uint32_t ll[16];
      for(int i = 0; i < 16; i++){tt[i] = 0; ll[i] = 0xFFFFFFFFu;}
      tt[0] = 5; ll[0] = 0;
      tt[1] = 0; ll[1] = sampleBytes;
      if(vrtypeStorage == 9){tt[2] = 8; tt[3] = 8; ll[2] = sampleBytes + 2; ll[3] = sampleBytes + 3;}
      else if(vrtypeStorage == 10){tt[2] = 16; tt[3] = 16; ll[2] = sampleBytes + 1; ll[3] = sampleBytes + 2;}
      else if(vrtypeStorage == 11){
         for(int i = 2; i < 8; i++){tt[i] = 8; ll[i] = sampleBytes + i; tt[i + 8] = 24; ll[i + 8] = sampleBytes + 1 + i;}
         tt[9] = 16; ll[9] = sampleBytes + 1;
         tt[8] = 0; ll[8] = 0xFFFFFFFFu;
      }else if(vrtypeStorage == 12){tt[2] = 0x45; ll[2] = 2; tt[3] = 0x40; ll[3] = 3;}
      else if(vrtypeStorage == 14){
         tt[3] = 0x10; ll[3] = 2; tt[4] = 0x45; ll[4] = 2; tt[5] = 0x40; ll[5] = 3;
         tt[12] = 0xc5; ll[12] = 4; tt[13] = 0xc0; ll[13] = 5; tt[15] = 0xe0; ll[15] = 6;
      }
      const unsigned mask = (1u << bw) - 1;
      uint64_t fpos = blockFpos[block];
      for(uint32_t i = 0; i < n; i++){
         const uint64_t bit = (uint64_t) i * bw;
         const unsigned x = (p[bit >> 3] >> (bit & 7)) & mask;
         if(ll[x] == 0xFFFFFFFFu){Fail(err, "Invalid record type and length code in the header of the .pgen file"); return false;}
         types[i] = tt[x];
         lengths[i] = ll[x];
         positions[i] = fpos;
         fpos += ll[x];
      }
      p += (n * bw + 7) / 8;
   }else{
      if(vrtypeStorage < 4){
         for(uint32_t i = 0; i < n; i++) types[i] = (p[i >> 1] >> (4 * (i & 1))) & 15;
         p += (n + 1) / 2;
      }else{
         memcpy(&types[0], p, n);
         p += n;
      }
      const unsigned lb = 1 + (vrtypeStorage & 3);
      uint64_t fpos = blockFpos[block];
      for(uint32_t i = 0; i < n; i++){
         uint32_t len = 0;
         for(unsigned k = 0; k < lb; k++) len |= (uint32_t) p[(size_t) i * lb + k] << (8 * k);
         lengths[i] = len;
         positions[i] = fpos;
         fpos += len;
      }
      p += (size_t) n * lb;
   }
   altCts.assign(n, 0);
   if(altCountBytes){
      for(uint32_t i = 0; i < n; i++) altCts[i] = p[i];
   }
   if(types[0] & 2 && (types[0] & 7) <= 3){        // the first variant of a block cannot be LD-compressed
      Fail(err, "Variant %u of the .pgen file is LD-compressed, which is not allowed at the start of a variant block", first);
      return false;
   }
   loadedBlock = block;
   return true;
}

bool PgenReader::ReadRecord(uint32_t v, std::string & err)
{
   (void) v; (void) err;
   return true;
}

// The main data track of a record that is not LD-compressed.
bool PgenReader::DecodeMainTrack(unsigned type, const unsigned char *p, const unsigned char *e,
   const std::vector<unsigned char> *, std::vector<unsigned char> & cats, const unsigned char **end, std::string & err)
{
   cats.assign(sampleBytes, 0);
   const unsigned t = type & 7;
   if(t == 0){
      if((uint64_t)(e - p) < sampleBytes){Fail(err, "A variant record is cut short"); return false;}
      memcpy(&cats[0], p, sampleBytes);
      *end = p + sampleBytes;
      return true;
   }
   if(t == 5){*end = p; return true;}
   std::vector<uint32_t> ids;
   std::vector<unsigned char> gens;
   if(t == 1){
      if(p >= e){Fail(err, "A variant record is cut short"); return false;}
      const int pair = OneBitPair(*p++);
      if(pair < 0){Fail(err, "A 1-bit variant record has an invalid pair of categories"); return false;}
      const uint32_t bitBytes = (N + 7) / 8;
      if((uint64_t)(e - p) < bitBytes){Fail(err, "A variant record is cut short"); return false;}
      const unsigned short *table = Tables().oneBit[pair];
      for(uint32_t b = 0; b < bitBytes; b++){
         const unsigned short v16 = table[p[b]];
         const uint32_t at = b * 2;
         cats[at] = (unsigned char)(v16 & 0xFF);
         if(at + 1 < sampleBytes) cats[at + 1] = (unsigned char)(v16 >> 8);
      }
      p += bitBytes;
   }else{
      const unsigned char fill = (t == 4)? 0x00: (t == 6)? 0xAA: 0xFF;       // t = 7: all category 3
      if(fill != 0) memset(&cats[0], fill, sampleBytes);
   }
   if(!ParseDifflist(p, e, N, true, ids, gens, err)) return false;
   for(size_t k = 0; k < ids.size(); k++) SetCategory(&cats[0], ids[k], gens[k]);
   *end = p;
   return true;
}

// The hard-calls of a record that has multiallelic patches (bit 3 of its type).  cats are the categories of the
// main track (REF-ALT1 and ALT1-ALT1 until patched); p..e is the multiallelic track and what comes after it.
bool PgenReader::DecodeMultiallelic(const unsigned char *cats, const unsigned char *p, const unsigned char *e,
   unsigned altCt, unsigned alt, unsigned char *out, std::string & err)
{
   if(altCt < 2 || altCt > 254){Fail(err, "A multiallelic record is read as a variant with %u ALT alleles", altCt); return false;}
   if(p >= e){Fail(err, "A variant record is cut short"); return false;}
   const unsigned fmt = *p++;
   const unsigned f1 = fmt & 15, f2 = fmt >> 4;
   // the two alleles of every sample: 255 = missing
   std::vector<unsigned char> a1(N), a2(N);
   std::vector<uint32_t> ids1, ids2;                 // the samples in category 1 and 2
   for(uint32_t i = 0; i < N; i++){
      switch(GetCategory(cats, i)){
         case 0: a1[i] = 0; a2[i] = 0; break;
         case 1: a1[i] = 0; a2[i] = 1; ids1.push_back(i); break;
         case 2: a1[i] = 1; a2[i] = 1; ids2.push_back(i); break;
         default: a1[i] = 255; a2[i] = 255; break;
      }
   }
   std::vector<uint32_t> ids;
   std::vector<unsigned char> none;
   // ---- the patches of the category 1 samples: REF-ALTk (k > 1) ----
   if(f1 != 15){
      std::vector<uint32_t> patched;
      if(f1 == 0){
         const uint32_t J = (uint32_t) ids1.size();
         const uint32_t bytes = (J + 7) / 8;
         if((uint64_t)(e - p) < bytes){Fail(err, "A multiallelic record is cut short"); return false;}
         for(uint32_t j = 0; j < J; j++) if((p[j >> 3] >> (j & 7)) & 1) patched.push_back(ids1[j]);
         p += bytes;
      }else if(f1 == 1){
         if(!ParseDifflist(p, e, N, false, patched, none, err)) return false;
      }else{Fail(err, "A multiallelic record has an unknown patch format (%u)", f1); return false;}
      const unsigned w = (altCt == 2)? 0: (altCt == 3)? 1: (altCt <= 5)? 2: (altCt <= 17)? 4: 8;
      const uint64_t bytes = ((uint64_t) patched.size() * w + 7) / 8;
      if((uint64_t)(e - p) < bytes){Fail(err, "A multiallelic record is cut short"); return false;}
      for(size_t k = 0; k < patched.size(); k++){
         const unsigned v = w? ReadBits(p, (uint64_t) k * w, w): 0;
         a2[patched[k]] = (unsigned char)(v + 2);
      }
      p += bytes;
   }
   // ---- the patches of the category 2 samples: ALT1-ALTk or ALTj-ALTk ----
   if(f2 != 15){
      std::vector<uint32_t> patched;
      if(f2 == 0){
         const uint32_t J = (uint32_t) ids2.size();
         const uint32_t bytes = (J + 7) / 8;
         if((uint64_t)(e - p) < bytes){Fail(err, "A multiallelic record is cut short"); return false;}
         for(uint32_t j = 0; j < J; j++) if((p[j >> 3] >> (j & 7)) & 1) patched.push_back(ids2[j]);
         p += bytes;
      }else if(f2 == 1){
         if(!ParseDifflist(p, e, N, false, patched, none, err)) return false;
      }else{Fail(err, "A multiallelic record has an unknown patch format (%u)", f2); return false;}
      const uint64_t K = patched.size();
      if(altCt == 2){
         const uint64_t bytes = (K + 7) / 8;
         if((uint64_t)(e - p) < bytes){Fail(err, "A multiallelic record is cut short"); return false;}
         for(size_t k = 0; k < K; k++){
            const bool homAlt2 = (p[k >> 3] >> (k & 7)) & 1;
            a1[patched[k]] = homAlt2? 2: 1;
            a2[patched[k]] = 2;
         }
         p += bytes;
      }else{
         const unsigned w = (altCt <= 4)? 2: (altCt <= 16)? 4: 8;
         const uint64_t bytes = (K * 2 * w + 7) / 8;
         if((uint64_t)(e - p) < bytes){Fail(err, "A multiallelic record is cut short"); return false;}
         const unsigned mask = (1u << w) - 1;
         for(size_t k = 0; k < K; k++){
            const unsigned v = ReadBits(p, (uint64_t) k * 2 * w, 2 * w);
            a1[patched[k]] = (unsigned char)((v & mask) + 1);
            a2[patched[k]] = (unsigned char)(((v >> w) & mask) + 1);
         }
         p += bytes;
      }
   }
   memset(out, 0, sampleBytes);
   for(uint32_t i = 0; i < N; i++){
      unsigned code;
      if(a1[i] == 255) code = 1;
      else{
         const unsigned copies = (a1[i] == alt) + (a2[i] == alt);
         code = (copies == 2)? 0: (copies == 1)? 2: 3;
      }
      out[i >> 2] |= (unsigned char)(code << ((i & 3) * 2));
   }
   return true;
}

bool PgenReader::ReadRow(uint32_t v, unsigned altCt, unsigned alt, unsigned char *out, std::string & err)
{
   if(v >= M){Fail(err, "Variant %u is beyond the end of the .pgen file", v); return false;}
   const PgenTables & tb = Tables();
   if(fixedWidth){
      std::vector<unsigned char> & rec = work;
      rec.resize(sampleBytes);
      if(PG_FSEEK(ff, (long long)(fixedStart + (uint64_t) v * fixedRecord), SEEK_SET) != 0
         || (sampleBytes > 0 && fread(&rec[0], 1, sampleBytes, ff) < sampleBytes)){
         Fail(err, "Not enough genotypes at variant %u of the .pgen file", v);
         return false;
      }
      if(mode == 1){
         if(alt == 1) memcpy(out, &rec[0], sampleBytes);
         else for(uint32_t b = 0; b < sampleBytes; b++){       // a PLINK 1 code: only the missing ones stay missing
            unsigned char x = 0;
            for(int j = 0; j < 4; j++) x |= (unsigned char)((((rec[b] >> (2 * j)) & 3) == 1? 1: 3) << (2 * j));
            out[b] = x;
         }
      }else{
         const unsigned char *table = (alt == 1)? tb.toPlink: tb.toPlinkNone;
         for(uint32_t b = 0; b < sampleBytes; b++) out[b] = table[rec[b]];
      }
      return true;
   }
   const uint32_t block = v / VBLOCK, i = v % VBLOCK;
   if(loadedBlock != (long long) block && !LoadBlock(block, err)) return false;
   const unsigned type = types[i];
   if(altCountBytes && altCts[i] != altCt + 1){
      Fail(err, "Variant %u has %u alleles in the .pgen file but %u in the .pvar file", v, (unsigned) altCts[i], altCt + 1);
      return false;
   }
   record.resize(lengths[i]);
   if(lengths[i] > 0 && (PG_FSEEK(ff, (long long) positions[i], SEEK_SET) != 0
      || fread(&record[0], 1, lengths[i], ff) < lengths[i])){
      Fail(err, "Not enough bytes for variant %u of the .pgen file", v);
      return false;
   }
   const unsigned char *p = lengths[i] > 0? &record[0]: NULL;
   const unsigned char *e = p + lengths[i];
   const unsigned t = type & 7;
   const unsigned char *mainEnd = p;
   const unsigned char *cats = NULL;
   if(t == 2 || t == 3){
      // LD-compressed: the differences from the latest variant that is not LD-compressed
      uint32_t j = i;
      do{
         if(j == 0){Fail(err, "Variant %u of the .pgen file is LD-compressed with no variant to refer to", v); return false;}
         j--;
      }while((types[j] & 7) == 2 || (types[j] & 7) == 3);
      const uint32_t baseV = block * VBLOCK + j;
      if(baseVariant != (long long) baseV){
         baseVariant = -1;
         work.resize(lengths[j]);
         if(lengths[j] > 0 && (PG_FSEEK(ff, (long long) positions[j], SEEK_SET) != 0
            || fread(&work[0], 1, lengths[j], ff) < lengths[j])){
            Fail(err, "Not enough bytes for variant %u of the .pgen file", baseV);
            return false;
         }
         const unsigned char *bp = lengths[j] > 0? &work[0]: NULL;
         const unsigned char *be0 = NULL;
         if(!DecodeMainTrack(types[j], bp, bp + lengths[j], NULL, baseCats, &be0, err)) return false;
         baseVariant = baseV;
      }
      std::vector<uint32_t> ids;
      std::vector<unsigned char> gens;
      const unsigned char *q = p;
      if(!ParseDifflist(q, e, N, true, ids, gens, err)) return false;
      baseCatsCopy = baseCats;
      for(size_t k = 0; k < ids.size(); k++) SetCategory(&baseCatsCopy[0], ids[k], gens[k]);
      if(t == 3) for(uint32_t b = 0; b < sampleBytes; b++) baseCatsCopy[b] = tb.invert[baseCatsCopy[b]];
      mainEnd = q;
      cats = &baseCatsCopy[0];
   }else{
      baseVariant = -1;
      if(!DecodeMainTrack(type, p, e, NULL, baseCats, &mainEnd, err)) return false;
      baseVariant = v;
      cats = &baseCats[0];
   }
   if(type & 8){
      return DecodeMultiallelic(cats, mainEnd, e, altCt, alt, out, err);
   }
   const unsigned char *table = (alt == 1)? tb.toPlink: tb.toPlinkNone;
   for(uint32_t b = 0; b < sampleBytes; b++) out[b] = table[cats[b]];
   return true;
}
