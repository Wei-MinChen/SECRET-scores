////////////////////////////////////////////////////////////////////////
// ScoreEngine.h
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
// secretscores calculates polygenic scores from PLINK binary genotypes
// (.bed/.bim/.fam) and score models, either PGS Catalog scoring files (text)
// or SECRETS binary models (.bin).
////////////////////////////////////////////////////////////////////////

#ifndef __SCOREENGINE_H__
#define __SCOREENGINE_H__

#include <stdio.h>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <zlib.h>

// 64-bit file positions (a .bed file of a biobank is far larger than 2 GB)
#ifdef _WIN32
  #define FSEEK64 _fseeki64
  #define FTELL64 _ftelli64
#else
  #define FSEEK64 fseeko
  #define FTELL64 ftello
#endif

// chromosome codes of the PLINK .bim file (human): 1-22 autosomes,
// 23 = X, 24 = Y, 25 = XY (pseudo-autosomal), 26 = MT
#define SEXCHR 23

// the chromosomes that are scored: autosomes and XY
inline bool ScoredChromosome(int chr) { return (chr > 0 && chr < SEXCHR) || chr == SEXCHR + 2; }

// chromosome code of a label ("1".."26", "X", "Y", "XY", "MT", optional "chr" prefix); 0 when
// it is not recognized
int ChromosomeCode(const char *label, size_t length);

// A code for an allele label that is compared instead of the label: the (upper case) letter
// for a one-letter allele, a 31-bit hash (high bit set) for longer alleles, so that the
// alleles of a whole .bim file take 4 bytes each.  0xFFFFFFFF means no allele.
const unsigned int NO_ALLELE = 0xFFFFFFFFu;
unsigned int AlleleCode(const char *label, size_t length);

// ---- models ----

// true for a model file name that ends in .bin (a SECRETS binary model; see BinaryInput.cpp)
bool IsBinaryModelName(const std::string & file);

// score name from a model file name: no directory, no .gz, no final extension
std::string ModelNameFromFile(const std::string & file);

// What the lines of a text model came to.
struct TextModelStats
{
   int fileSNPs;                 // lines with a usable effect allele and weight
   int malformed;                // lines without
   int noCoord;                  // of fileSNPs: no usable chromosome (1-26) or position
   int noOtherAllele;            // of fileSNPs: neither other_allele nor hm_inferOtherAllele has a value
   int sameAllele;               // of fileSNPs: the other allele is the effect allele
};

// Parse a PGS Catalog text model (plain or gzipped) without any genotype data.  For every line
// with a usable effect allele, weight, chromosome (1-26), position and a different other allele
// (other_allele, else hm_inferOtherAllele), onSNP is called with (chromosome code, position,
// effect allele, other allele, weight); the alleles are in upper case.  The other lines are
// dropped and counted in stats.
typedef std::function<void(int, long, const std::string &, const std::string &, double)> SNPCallback;
void ParseTextModel(const std::string & file, TextModelStats & stats, const SNPCallback & onSNP);

// --makebin: convert text models to SECRETS binary models, <name>.bin for every model
// (names[i] if given, else the model file name without directory and extension).  Needs no
// genotypes.
void MakeBinaryModels(const std::vector<std::string> & files, const std::vector<std::string> & names);

// Reads a text file (plain or gzipped) line by line, in large blocks.
class LineReader
{
   private:
      gzFile gz;
      std::vector<char> buffer;
      size_t begin, end;            // unread data: buffer[begin, end)
      bool finished;
   public:
      LineReader();
      ~LineReader();
      bool Open(const char *filename);
      // The next line, without its line break and ending in '\0'; the pointer is valid until the
      // next call.  False at the end of the file.
      bool Next(const char *& line, size_t & length);
      void Close();
};

// One score model: the SNPs of a model file that are in the genotype data.
struct ScoreModel
{
   std::string name;             // name of the score (column of the score file)
   std::string file;             // model file
   bool binary;                  // SECRETS binary model (.bin), else PGS Catalog text
   // the SNPs of the model that match the data, in the order of the model file
   std::vector<int> rowIdx;      // .bim (and .bed) row
   std::vector<double> weights;  // effect weight
   std::vector<char> switched;   // 1: the effect allele is allele A2 of the data, 0: it is A1
   // what happened to the SNPs of the model file
   int fileSNPs;                 // SNPs with a usable effect allele and weight
   int malformed;                // lines (SNPs) without a usable effect allele or weight
   int noCoord;                  // no (usable) chromosome / position
   int noOtherAllele;            // no other allele (other_allele / hm_inferOtherAllele): dropped
   int sameAllele;               // other allele is the effect allele: dropped
   int nonAutosomal;             // on a chromosome that is not scored (X, Y, MT)
   int notInData;                // chromosome and position are not in the data
   int alleleMismatch;           // position is in the data, the alleles are not
   bool scored;                  // at least one SNP was scored
   int used;                     // SNPs scored
   // the SNPs that are scored, prepared by PrepareModel (in the order of the model file)
   std::vector<int> vslot;       // SNP index in the batch
   std::vector<double> vval;     // 4 per SNP: score contribution for genotype code 0,1,2,3
};

class ScoreEngine
{
   protected:
      // ---- samples, from the PLINK .fam file (in file order) ----
      std::vector<std::string> FID, IID;
      int idCount;

      // ---- SNPs, from the PLINK .bim file: every row, in file order, so that
      //      row r of the .bim is row r of the .bed.  Kept compactly (13 bytes per
      //      SNP): the SNP names are not needed. ----
      std::string plinkinput;
      std::vector<unsigned char> bimChr;     // chromosome code, 0 when not recognized
      std::vector<int> bimPos;               // base-pair position
      std::vector<unsigned int> bimA1, bimA2;// allele 1 (A1) and allele 2 (A2), as AlleleCode
      int bimCount;
      bool bimSorted;                        // rows are in (chromosome, position) order
      // when they are not: (chromosome << 32 | position, row) of the SNPs that can be scored, sorted
      std::vector<std::pair<long long, int> > posIndex;

      // ---- the genotypes of the SNPs being scored, INDIVIDUAL-major ----
      // G_IND[i] holds person i's genotypes at the loaded SNPs, 2 bits each (PLINK
      // coding: 00 = 2 copies of A1, 01 = missing, 10 = one copy, 11 = no copy),
      // SNP k of the batch in byte k/4, bits 2*(k%4).
      unsigned char **G_IND;
      int loadedCount;              // SNPs loaded (the batch)
      int genoBytes;                // bytes per person: (loadedCount+3)/4
      int bedMode;                  // storage mode of the .bed file: 1 = SNP-major, 0 = individual-major
      std::vector<int> slotOfRow;   // .bim row -> SNP index in the batch, -1 when not loaded
      int threadCount;              // threads for the people-level loops

      // The SNP of the data at chromosome chr and position pos whose alleles are the effect
      // allele ea and the other allele oa (in either order).  Returns its .bim row, or -1; sw
      // tells whether the effect allele is A2; atPosition tells whether the data have any SNP
      // at the position.
      int MatchSNP(int chr, long pos, unsigned int ea, unsigned int oa, bool & sw, bool & atPosition) const;
      void ReadModelFile(ScoreModel & model);          // PGS Catalog text model
      void ReadBinaryModelFile(ScoreModel & model);    // SECRETS binary model (BinaryInput.cpp)
      FILE *OpenPLINKBinary();
      void ReadPLINKBinary(FILE *fp, std::vector<int> & rows);
      void ReadPLINKBinarySNPMajor(FILE *fp, std::vector<int> & rows);
      void FreeGenotypes(std::vector<int> & rows);
      // Tabulate what each genotype code of each matched SNP contributes.  Without SNPs to
      // score model.scored is false.
      void PrepareModel(ScoreModel & model, bool verbose);
      // Score the models first..last-1 (prepared) for every person.  scores[k][i] gets the
      // score; missSNP (optional, first model only) gets the fraction of the SNPs that was
      // missing.
      void ScoreBatch(ScoreModel *models, int first, int last, std::vector<std::vector<double> > & scores,
         std::vector<double> * missSNP);

   public:
      ScoreEngine();
      ~ScoreEngine();

      std::string prefix;           // output file prefix
      int CoreCount;                // --cpus, 0 = automatic
      std::string optionReadModel;  // --readmodel and --readlist as given, for the "Options in effect"
      std::string optionReadList;

      void ReadPlinkFam(const std::string & plinkname);
      void ReadPlinkBim(const std::string & plinkname);

      // Read an individual-major PLINK binary file (.bed starting 0x6c 0x1b 0x00: the genotypes of
      // one person at all .bim SNPs, then the next person) into the individual-major storage
      // G_IND.  rows lists the .bim rows to keep, ascending without repeats; SNP k of the
      // storage is rows[k].  Persons are read in large blocks and unpacked in parallel.
      void ReadPLINKBinaryIndividualMajor(FILE *fp, std::vector<int> & rows);

      // Score every person with every model.  Writes <prefix>_scores.txt (the score
      // file SECRET reads: one row per person, one column per model) and
      // <prefix>_scores.coverage; with a single model also <prefix>grs.txt.
      void PrintGeneticRiskScore(const std::string & binfile, std::vector<std::string> & modelFiles,
         std::vector<std::string> & modelNames);
      char *currentTime();
};

#endif
