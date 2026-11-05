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
#include <stdexcept>
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

// How many models go into a batch, besides the memory limit (GENOTYPE_MEMORY_LIMIT): a batch loads
// the union of the SNPs of its models ("touched" SNPs), so the .bed file is read once per batch,
// but thousands of models in one batch make the loaded set huge and the scoring of each person
// slow.  Models are added to a batch until it has MAX_MODELS_PER_BATCH models, or it has at least
// MIN_MODELS_PER_BATCH models and the touched SNPs have reached TARGET_TOUCHED_SNPS_PER_BATCH (the
// model that reaches the target is kept in the batch).  The memory limit comes first: a model that
// would exceed it goes to the next batch, even if the batch then has fewer than the minimum.
// --batchmodels (debug build) changes the maximum (and the minimum, if it is lower).
#define MIN_MODELS_PER_BATCH 100
#define MAX_MODELS_PER_BATCH 1000
#define TARGET_TOUCHED_SNPS_PER_BATCH 9000000

// chromosome codes of the PLINK .bim file (human): 1-22 autosomes,
// 23 = X, 24 = Y, 25 = XY (pseudo-autosomal), 26 = MT
#define SEXCHR 23

// How the X chromosome (code 23) is scored, set by --allowX (the numbers are PLINK 2's):
//   0  X is not scored (the default)
//   1  females 0/1/2 copies of the effect allele, males 0/1 (hemizygous)
//   2  females and males both 0/1/2: a male's hemizygous call counts as two copies (as in
//      PLINK 2)
// In both 1 and 2 a heterozygous genotype of a male is impossible outside the PAR and is treated
// as missing (as PLINK 2 does).  The PAR must be coded XY (25), as PLINK does (--split-x), and
// is scored as an autosome.
extern int XChrModel;

// Mean imputation of missing genotypes (on by default, off with --no-mean-imputation): a missing
// genotype counts as 2 x the frequency of the effect allele (the mean dosage) instead of 0.
// The frequencies come from --readfreq (a SNP that the file lacks gets no imputation, with a warning),
// else from the allelefrequency_effect column of the model.
extern bool MeanImputation;

// the chromosomes that are scored: autosomes, XY, and X unless --allowX 0
inline bool ScoredChromosome(int chr)
{
   return (chr > 0 && chr < SEXCHR) || chr == SEXCHR + 2 || (chr == SEXCHR && XChrModel != 0);
}

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
   bool hasFreq;                 // the model has an allelefrequency_effect column
   // trait annotation from the "#key=value" header lines of a PGS Catalog scoring file (empty
   // when the file has no such line)
   std::string pgsId;            // #pgs_id
   std::string pgsName;          // #pgs_name
   std::string traitReported;    // #trait_reported
   std::string traitMapped;      // #trait_mapped
   std::string traitEfo;         // #trait_efo
};

// The trait of a model for the "trait" column of the coverage file: trait_reported, else
// trait_mapped, else trait_efo, else empty.
std::string ModelTrait(const TextModelStats & stats);

// One line of the coverage file (<prefix>_scores.coverage), one per model.  A count that is -1
// is not known: with --makebin there are no genotypes, so Found, Used, NotInData, AlleleMismatch
// and the coverage are not known and are not written.
struct CoverageRow
{
   std::string name;             // score name
   std::string file;             // model file
   std::string trait;            // trait annotation, empty = none
   std::string traitEfo;         // EFO id(s) of the trait (#trait_efo), empty = none
   long long fileSNPs;
   long long found, used;
   double percentImputed;        // % of the genotypes (Used x individuals) that were imputed, < 0: not known
   long long noCoord, noOtherAllele, sameAllele, nonAutosomal, notInData, alleleMismatch;
   CoverageRow() : fileSNPs(0), found(-1), used(-1), percentImputed(-1), noCoord(0), noOtherAllele(0), sameAllele(0),
      nonAutosomal(0), notInData(-1), alleleMismatch(-1) {}
};

// Write the coverage file: comma-delimited with a header row (as the output of the Michigan
// Imputation Server PGS calculation), one row per model, fields with a comma, quote or line
// break in double quotes.  Columns:
//   SCORE,SNPs,Found,Used,coverage,NoCoord,NoOtherAllele,SameAllele,NonAutosomal,NotInData,AlleleMismatch,trait,trait_efo
// and PercentImputed (before trait) unless --no-mean-imputation: the percentage of the genotypes
// (Used SNPs x individuals) that were imputed, i.e. missing at a SNP that has a frequency.
// NonAutosomal (the name is kept for the programs that read the file) is the number of SNPs on
// chromosomes that are not scored: Y and MT, and X with --allowX 0.
// coverage = Used / SNPs.  With withGenotypes false (--makebin) the columns that need
// genotypes (Found, Used, coverage, NotInData, AlleleMismatch) are left out.
void WriteCoverageFile(const std::string & filename, const std::vector<CoverageRow> & rows, bool withGenotypes);

// Parse a PGS Catalog text model (plain or gzipped) without any genotype data.  For every line
// with a usable effect allele, weight, chromosome (1-26), position and a different other allele
// (other_allele, else hm_inferOtherAllele), onSNP is called with (chromosome code, position,
// effect allele, other allele, weight, effect allele frequency); the alleles are in upper case.
// The frequency (column allelefrequency_effect) is NaN when the model has no such column or
// the cell is missing or not in [0, 1].  The other lines are dropped and counted in stats.
typedef std::function<void(int, long, const std::string &, const std::string &, double, double)> SNPCallback;
void ParseTextModel(const std::string & file, TextModelStats & stats, const SNPCallback & onSNP);

// A failure that concerns one model only: the file cannot be opened, its header lacks a
// required column, or a binary model is damaged.  The model is skipped with a note, and the
// other models go on (see WriteSkippedFile).  Thrown by ThrowModelError, caught per model.
class ModelError : public std::runtime_error
{
   public:
      explicit ModelError(const std::string & message) : std::runtime_error(message) {}
};
void ThrowModelError(const char *format, ...);

// A message on a single line (line breaks and runs of blanks become one blank), for the screen
// and for the skipped-models file.
std::string OneLine(const std::string & message);

// <prefix>_scores.skipped: tab-delimited FILE and REASON of the models that were skipped.
void WriteSkippedFile(const std::string & filename, const std::vector<std::string> & files,
   const std::vector<std::string> & reasons);

// --makebin: convert text models to SECRETS binary models, <name>.bin for every model
// (names[i] if given, else the model file name without directory and extension).  Needs no
// genotypes.  The trait annotation of the text models is kept in the binary models, and one
// file <prefix>_scores.coverage describes all the models that were made.
// The models are converted by `threads` threads (one model per thread at a time).
void MakeBinaryModels(const std::vector<std::string> & files, const std::vector<std::string> & names,
   const std::string & prefix, int threads);

// The number of threads for --cpus n (n > 0: n; 0: half of the available cores, at least 1;
// 1 when the program is built without OpenMP).
int ResolveThreads(int cpus);

// --viewbin: print a binary model as a text model (tab-delimited, with the trait as PGS Catalog
// header lines) to the screen, the reverse of --makebin.  Only the model is written to stdout.
void ViewBinaryModel(const std::string & file);

// Reads a text file (plain or gzipped) line by line, in large blocks.
class LineReader
{
   private:
      gzFile gz;
      std::vector<char> buffer;
      size_t begin, end;            // unread data: buffer[begin, end)
      bool finished;
      bool bad;                     // a read error (a truncated or damaged .gz file), not a clean end
      std::string badMessage;
   public:
      LineReader();
      ~LineReader();
      bool Open(const char *filename);
      // The next line, without its line break and ending in '\0'; the pointer is valid until the
      // next call.  False at the end of the file.
      bool Next(const char *& line, size_t & length);
      void Close();
      // True when Next() stopped because of a read error (e.g. a truncated .gz file) and not at the
      // end of the file; message says what zlib reported.  Check it after the loop over the lines.
      bool Failed(std::string & message) const { message = badMessage; return bad; }
};

// One score model: the SNPs of a model file that are in the genotype data.
struct ScoreModel
{
   std::string failed;           // why the model was skipped (cannot be read); empty = read
   std::string name;             // name of the score (column of the score file)
   std::string file;             // model file
   std::string trait;            // trait annotation (text model: header lines, binary model: kept by --makebin)
   std::string traitEfo;         // trait_efo, the same way
   bool binary;                  // SECRETS binary model (.bin), else PGS Catalog text
   // the SNPs of the model that match the data, in the order of the model file
   std::vector<int> rowIdx;      // .bim (and .bed) row
   std::vector<double> weights;  // effect weight
   std::vector<char> switched;   // 1: the effect allele is allele A2 of the data, 0: it is A1
   std::vector<double> freq;     // frequency of the effect allele in the model (NaN: none)
   // what happened to the SNPs of the model file
   int fileSNPs;                 // SNPs with a usable effect allele and weight
   int malformed;                // lines (SNPs) without a usable effect allele or weight
   int noCoord;                  // no (usable) chromosome / position
   int noOtherAllele;            // no other allele (other_allele / hm_inferOtherAllele): dropped
   int sameAllele;               // other allele is the effect allele: dropped
   int nonAutosomal;             // on a chromosome that is not scored (Y, MT; also X with --allowX 0)
   int notInData;                // chromosome and position are not in the data
   int alleleMismatch;           // position is in the data, the alleles are not
   bool scored;                  // at least one SNP was scored
   int used;                     // SNPs scored (autosomes, XY and X)
   // where the frequency for mean imputation of the scored SNPs came from (PrepareModel)
   int freqFromFile;             // --readfreq
   int freqFromModel;            // allelefrequency_effect of the model
   int freqNone;                 // none: a missing genotype counts as 0
   int freqSkipped;              // SNPs scored without imputation: --readfreq has no frequency for them
   // which SNPs can be imputed: 0 none have a frequency, 1 all have, 2 some (vhas/xhas flag them)
   int impMode;
   std::vector<char> vhas, xhas; // only for impMode 2: 1 = the SNP of vslot/xslot has a frequency
   long long imputedGenotypes;   // genotypes that were imputed (missing at a SNP with a frequency; ScoreBatch)
   // the SNPs that are scored, prepared by PrepareModel (in the order of the model file)
   std::vector<int> vslot;       // SNP index in the batch (not X)
   std::vector<double> vval;     // 4 per SNP: score contribution for genotype code 0,1,2,3
   // the X-chromosome SNPs that are scored (in the order of the model file), added to the score
   // after the other SNPs; the contribution of a male differs from that of the others
   std::vector<int> xslot;       // SNP index in the batch
   std::vector<double> xval;     // 4 per SNP, females and people of unknown sex
   std::vector<double> xvalM;    // 4 per SNP, males
};

class ScoreEngine
{
   protected:
      // ---- samples, from the PLINK .fam file (in file order) ----
      std::vector<std::string> FID, IID;
      std::vector<unsigned char> sex;        // from the .fam file: 1 male, 2 female, 0 unknown
      int idCount;
      long long maleHetX;                    // heterozygous genotypes of males at X SNPs, taken as missing

      // ---- SNPs, from the PLINK .bim file: every row, in file order, so that
      //      row r of the .bim is row r of the .bed.  Kept compactly (13 bytes per
      //      SNP): the SNP names are not needed. ----
      std::string plinkinput;
      std::vector<unsigned char> bimChr;     // chromosome code, 0 when not recognized
      std::vector<int> bimPos;               // base-pair position
      std::vector<unsigned int> bimA1, bimA2;// allele 1 (A1) and allele 2 (A2), as AlleleCode
      // --readfreq: the frequency of allele A1 of each .bim row (NaN: not given), matched by SNP ID
      std::vector<unsigned long long> bimIdHash;   // hash of the SNP ID, only kept for --readfreq
      std::vector<float> bimFreqA1;
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
      long long readBlockBytes;     // size of one read of the individual-major .bed file (two are held)

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
      bool haveReadFreq;            // --readfreq was read
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
      int readBlockMB;              // --blocksize: MB of the .bed file read at a time (0 = default)
      std::string optionReadModel;  // --readmodel and --readlist as given, for the "Options in effect"
      std::string optionReadList;
      int maxBatchModels;           // the most models in one batch (default MAX_MODELS_PER_BATCH)
      std::string readFreqFile;     // --readfreq (empty = none); read after the .bim file

      void ReadPlinkFam(const std::string & plinkname);
      void ReadPlinkBim(const std::string & plinkname);
      // --readfreq: allele frequencies of the genotype data, from a PLINK --freq file (.frq:
      // columns SNP, A1, MAF = frequency of A1) or a PLINK 2 --freq file (.afreq: ID, ALT,
      // ALT_FREQS = frequency of ALT), matched to the .bim rows by SNP ID.  Plain or gzipped.
      void ReadFreqFile(const std::string & filename);

      // Read an individual-major PLINK binary file (.bed starting 0x6c 0x1b 0x00: the genotypes of
      // one person at all .bim SNPs, then the next person) into the individual-major storage
      // G_IND.  rows lists the .bim rows to keep, ascending without repeats; SNP k of the
      // storage is rows[k].  Persons are read in large blocks (readBlockBytes, --blocksize); the next
      // block is read while the current one is unpacked, and the unpacking is split into
      // (person, SNPs) tasks among the threads, so a block of a few persons keeps them all busy.
      void ReadPLINKBinaryIndividualMajor(FILE *fp, std::vector<int> & rows);

      // Score every person with every model.  Writes <prefix>_scores.txt (the score
      // file SECRET reads: one row per person, one column per model) and
      // <prefix>_scores.coverage (comma-delimited, with the trait of each model); with a single
      // model also <prefix>grs.txt.
      void PrintGeneticRiskScore(const std::string & binfile, std::vector<std::string> & modelFiles,
         std::vector<std::string> & modelNames);
      char *currentTime();
};

#endif
