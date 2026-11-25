////////////////////////////////////////////////////////////////////////
// Main.cpp
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
// secretscores: polygenic score calculation from PLINK binary genotypes and
// PGS Catalog scoring files.
////////////////////////////////////////////////////////////////////////

// release: VERSION is defined; debug mode: comment it out (the banner then shows the date instead of
// the version, and the debug-only options --blocksize, --batchmodels and --sex-maxsnps are available)
#define VERSION "1.0.4"

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include "Parameters.h"
#include "ScoreEngine.h"
#include <zlib.h>
#include "Error.h"

void ShowBanner(bool toStderr)
{
#ifdef VERSION
   fprintf(toStderr? stderr: stdout, "secretscores " VERSION " - (c) 2026-2026 Wei-Min Chen\n");
#else
   // debug mode: the current date in place of the version number
   char today[32];
   time_t now = time(NULL);
   strftime(today, sizeof(today), "%Y-%m-%d", localtime(&now));
   fprintf(toStderr? stderr: stdout, "secretscores %s - (c) 2026-2026 Wei-Min Chen\n", today);
#endif
}

// --readlist: one model per line.  The first column is the model file name (a path is allowed,
// relative to the current directory), and the optional second column is the name of the score
// (a column of the score file; without it the name comes from the model file name).  Columns
// are separated by white space or, if the line has a tab, by tabs only (so that a file name or
// a score name may have spaces).  Empty lines and lines that start with # are skipped.
static std::string TrimBlanks(const std::string & s)
{
   size_t b = 0, e = s.size();
   while(b < e && (s[b] == ' ' || s[b] == '\t' || s[b] == '\r' || s[b] == '\n')) b++;
   while(e > b && (s[e-1] == ' ' || s[e-1] == '\t' || s[e-1] == '\r' || s[e-1] == '\n')) e--;
   return s.substr(b, e - b);
}

// "a,b,,c" -> a, b, c (empty items are dropped)
static std::vector<std::string> SplitCommas(const char *text)
{
   std::vector<std::string> items;
   std::string item;
   for(const char *p = text; ; p++){
      if(*p == ',' || *p == 0){
         if(!item.empty()) items.push_back(item);
         item.clear();
         if(*p == 0) break;
      }else item += *p;
   }
   return items;
}

// the next line of a text file (plain or .gz) without the end of line; false at the end of the file
static bool ReadTextLine(gzFile input, std::string & line)
{
   line.clear();
   char buffer[4096];
   bool any = false;
   while(gzgets(input, buffer, sizeof(buffer)) != NULL){
      any = true;
      line += buffer;
      if(!line.empty() && line[line.size() - 1] == '\n') break;
   }
   while(!line.empty() && (line[line.size() - 1] == '\n' || line[line.size() - 1] == '\r')) line.erase(line.size() - 1);
   return any;
}

static void ReadModelList(const char *listfile, std::vector<std::string> & files, std::vector<std::string> & names)
{
   gzFile input = gzopen(listfile, "rb");
   if(input == NULL) error("Cannot open %s to read", listfile);
   std::string text;
   int lineNumber = 0, listed = 0;
   while(ReadTextLine(input, text)){
      lineNumber++;
      if(lineNumber == 1 && text.compare(0, 3, "\xEF\xBB\xBF") == 0) text.erase(0, 3);    // a UTF-8 byte order mark
      // blank lines and comments are skipped; trailing blanks and leading spaces go, but a leading
      // tab is an empty first column (an error, as the first column is always the model file)
      const size_t firstChar = text.find_first_not_of(" \t\r\n");
      if(firstChar == std::string::npos || text[firstChar] == '#') continue;
      text.erase(text.find_last_not_of(" \t\r\n") + 1);
      text.erase(0, text.find_first_not_of(' '));
      std::vector<std::string> column;
      if(text.find('\t') != std::string::npos){
         size_t from = 0;
         for(;;){
            size_t to = text.find('\t', from);
            column.push_back(TrimBlanks(text.substr(from, to == std::string::npos? std::string::npos: to - from)));
            if(to == std::string::npos) break;
            from = to + 1;
         }
         while(column.size() > 1 && column.back().empty()) column.pop_back();
      }else{
         size_t from = 0;
         while(from < text.size()){
            size_t to = text.find_first_of(" \t", from);
            if(to == std::string::npos) to = text.size();
            column.push_back(text.substr(from, to - from));
            from = text.find_first_not_of(" \t", to);
            if(from == std::string::npos) break;
         }
      }
      if(column.empty() || column[0].empty())
         error("Line %d of %s has no model file name in its first column.", lineNumber, listfile);
      if(column.size() > 2)
         error("Line %d of %s has %d columns: %s\n"
            "  A --readlist file has the model file name in the first column and, optionally, the name of the\n"
            "  score in the second column. If a file name or a score name has spaces, separate the two columns\n"
            "  with a tab.", lineNumber, listfile, (int) column.size(), text.c_str());
      files.push_back(column[0]);
      names.push_back(column.size() > 1? column[1]: std::string());
      listed++;
   }
   gzclose(input);
   if(listed == 0) error("No model files are listed in %s", listfile);
}

// --viewbin writes the text model to stdout, so its messages go to stderr
static void ViewbinFail(const char *message)
{
   fprintf(stderr, "\nFATAL ERROR - \n%s\n\n", message);
   exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
   // line-buffered output: the messages reach the screen, a log file or a pipe as they are
   // printed, not when a buffer fills (a long run to a redirected file shows nothing otherwise)
   // --viewbin: stdout is the text model only (the banner goes to stderr), with a large buffer
   bool viewbinArg = false;
   for(int i = 1; i < argc; i++)
      if(strcmp(argv[i], "--viewbin") == 0) viewbinArg = true;
   setvbuf(stdout, NULL, viewbinArg? _IOFBF: _IOLBF, viewbinArg? 1 << 16: 0);
   ShowBanner(viewbinArg);

   OptionString binfile;
   OptionString vcffile("");
   OptionString pfilePrefix("");
   ParameterList pl;
   bool makebinFlag = false;
   bool viewbinFlag = false;
   bool scoreFlag = false;
   bool inferABOFlag = false;
   bool inferAPOEFlag = false;
   bool apoePartialFlag = false;
   bool inferSexFlag = false;
   int sexMaxSnps = 50000;
   OptionString weightfile("");
   OptionString listfile("");
   OptionString freqfile("");
   OptionString keepfile("");
   OptionString removefile("");
   bool noMeanImputation = false;
   OptionString prefix("secretscores");
   int CoreCount = 0;
   int allowX = 0;
   int blockSize = 0;
   int batchModels = MAX_MODELS_PER_BATCH;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Input Files")
         LONG_STRINGPARAMETER("bfile", &binfile)
         LONG_STRINGPARAMETER("vcf", &vcffile)
         LONG_STRINGPARAMETER("pfile", &pfilePrefix)
         LONG_STRINGPARAMETER("readmodel", &weightfile)
         LONG_STRINGPARAMETER("readlist", &listfile)
         LONG_STRINGPARAMETER("readfreq", &freqfile)
      LONG_PARAMETER_GROUP("Sample Filter")
         LONG_STRINGPARAMETER("keep", &keepfile)
         LONG_STRINGPARAMETER("remove", &removefile)
      LONG_PARAMETER_GROUP("Score Action")
         LONG_PARAMETER("score", &scoreFlag)
         LONG_PARAMETER("makebin", &makebinFlag)
         LONG_PARAMETER("viewbin", &viewbinFlag)
      LONG_PARAMETER_GROUP("Inferences")
         LONG_PARAMETER("inferABO", &inferABOFlag)
         LONG_PARAMETER("inferAPOE", &inferAPOEFlag)
         LONG_PARAMETER("apoe-partial", &apoePartialFlag)
         LONG_PARAMETER("infersex", &inferSexFlag)
#ifndef VERSION
         LONG_INTPARAMETER("sex-maxsnps", &sexMaxSnps)
#endif
      LONG_PARAMETER_GROUP("Score Parameter")
         LONG_PARAMETER("no-mean-imputation", &noMeanImputation)
         LONG_INTPARAMETER("allowX", &allowX)
      LONG_PARAMETER_GROUP("Computing Parameter")
         LONG_INTPARAMETER("cpus", &CoreCount)
#ifndef VERSION
         LONG_INTPARAMETER("blocksize", &blockSize)
         LONG_INTPARAMETER("batchmodels", &batchModels)
#endif
      LONG_PARAMETER_GROUP("Output Parameter")
         LONG_STRINGPARAMETER("prefix", &prefix)
   END_LONG_PARAMETERS()

   pl.Add(new LongParameters("Options", longParameters));
   pl.Read(argc, argv);
   if(!viewbinFlag) pl.Status();

   // the three modes: --score (calculate scores; also the default when none is given),
   // --makebin and --viewbin; only one of them at a time
   if(scoreFlag && viewbinFlag) ViewbinFail("--score and --viewbin cannot be used together.");
   if(scoreFlag && makebinFlag)
      error("--score and --makebin cannot be used together.");

   if(viewbinFlag){
      // print one binary model as a text model: secretscores --viewbin --readmodel model.bin
      if(makebinFlag) ViewbinFail("--viewbin and --makebin cannot be used together.");
      if(!listfile.IsEmpty()) ViewbinFail("--viewbin shows a single model: use --readmodel <file.bin>, not --readlist.");
      const std::vector<std::string> files = SplitCommas(weightfile);
      const int fileCount = (int) files.size();
      if(fileCount != 1)
         ViewbinFail("--viewbin shows a single model. Please use --readmodel <file.bin> with one binary model,\n"
            "  e.g.,  secretscores --viewbin --readmodel model.bin | more\n"
            "         secretscores --viewbin --readmodel model.bin > model.txt");
      for(size_t i = 0; i < files.size(); i++){
         if(!IsBinaryModelName(files[i]))
            ViewbinFail("--viewbin needs a SECRETS binary model (file name ending in .bin, made with --makebin).");
         ViewBinaryModel(files[i]);
      }
      return 0;
   }

   {
      const int inputs = (binfile.IsEmpty()? 0: 1) + (vcffile.IsEmpty()? 0: 1) + (pfilePrefix.IsEmpty()? 0: 1);
      if(inputs > 1) error("Use only one of --bfile (PLINK binary), --vcf (VCF) and --pfile (PLINK 2 pgen) for the genotypes.");
   }
   if(binfile.IsEmpty() && vcffile.IsEmpty() && pfilePrefix.IsEmpty() && !makebinFlag)
      error("Genotype files are required. e.g.,\n  secretscores --score --bfile ex.bed --readmodel model.txt\n"
         "  secretscores --score --vcf ex.vcf.gz --readmodel model.txt\n"
         "  secretscores --score --pfile ex --readmodel model.txt   (ex.pgen, ex.pvar, ex.psam)\n"
         "  (--score is optional: without --makebin or --viewbin, scores are calculated)\n\n"
         "The model file is a PGS Catalog harmonized scoring file. Its header needs the columns\n"
         "  effect_allele, effect_weight, hm_chr, hm_pos\n"
         "  and other_allele or hm_inferOtherAllele; or a SECRETS binary model (file name ending in .bin,\n"
         "  made from a text model with --makebin).");
   if(apoePartialFlag && !inferAPOEFlag)
      error("--apoe-partial is used with --inferAPOE.");
   if((int) inferABOFlag + (int) inferAPOEFlag + (int) inferSexFlag > 1)
      error("--inferABO, --inferAPOE and --infersex are run one at a time.");
   if(sexMaxSnps < 20) error("--sex-maxsnps %d is out of range (at least 20).", sexMaxSnps);
   if((inferABOFlag || inferAPOEFlag || inferSexFlag) && (scoreFlag || makebinFlag))
      error("--inferABO, --inferAPOE and --infersex cannot be used with --score or --makebin.");
   if((inferABOFlag || inferAPOEFlag || inferSexFlag) && binfile.IsEmpty() && vcffile.IsEmpty() && pfilePrefix.IsEmpty())
      error("--inferABO, --inferAPOE and --infersex need genotypes: --bfile, --vcf or --pfile.");
   if(!inferABOFlag && !inferAPOEFlag && !inferSexFlag && weightfile.IsEmpty() && listfile.IsEmpty())
      error("Please use --readmodel <file> to specify a score model,\n  or --readlist <file> to specify a file with the names of score models.");
   if(CoreCount < 0)
      error("--cpus %d is out of range.", CoreCount);
   if(allowX < 0 || allowX > 2)
      error("--allowX %d is out of range. Use 0 (the X chromosome is not scored, the default), 1 (males 0/1,\n"
         "  females 0/1/2) or 2 (males and females 0/1/2).", allowX);
   if(blockSize < 0 || blockSize > 65536)
      error("--blocksize %d is out of range. Use the number of MB of the .bed file that are read at a time (1-65536).", blockSize);
   if(batchModels < 1 || batchModels > 1000000)
      error("--batchmodels %d is out of range. Use the most models in one batch (1-1000000).", batchModels);
   XChrModel = allowX;
   MeanImputation = !noMeanImputation;

   time_t timer;
   struct tm *tblock;
   timer = time(NULL);
   tblock = localtime(&timer);
   printf("secretscores starts at %s", asctime(tblock));

   // the models: --readmodel file1,file2,...  and  --readlist file (one model file per line,
   // with an optional second column that names the score; see ReadModelList).  A file whose name ends in .bin
   // is a SECRETS binary model, any other file is a PGS Catalog text model.
   std::vector<std::string> modelFiles, modelNames;
   if(!weightfile.IsEmpty()){
      const std::vector<std::string> files = SplitCommas(weightfile);
      for(size_t i = 0; i < files.size(); i++){
         modelFiles.push_back(files[i]);
         modelNames.push_back("");
      }
   }
   if(!listfile.IsEmpty())
      ReadModelList((const char*) listfile, modelFiles, modelNames);
   if(modelFiles.empty() && !inferABOFlag && !inferAPOEFlag && !inferSexFlag)
      error("Please use --readmodel <file> to specify a score model,\n  or --readlist <file> to specify a file with the names of score models.");

   if(makebinFlag){
      // convert text models to binary models: no genotypes are read
      if(!binfile.IsEmpty() || !vcffile.IsEmpty() || !pfilePrefix.IsEmpty())
         printf("--makebin only converts the models; the genotype data (--bfile, --vcf, --pfile) are not used.\n");
      if(!keepfile.IsEmpty() || !removefile.IsEmpty())
         printf("--makebin only converts the models; --keep and --remove are not used.\n");
      if(!freqfile.IsEmpty())
         printf("--makebin only converts the models; --readfreq is not used (the frequencies in the\n"
                "  allelefrequency_effect column of the models are stored in the binary models).\n");
      const int makeThreads = ResolveThreads(CoreCount);
      if(makeThreads > 1) printf("%d threads\n", makeThreads);
      MakeBinaryModels(modelFiles, modelNames, (const char*) prefix, makeThreads);
      timer = time(NULL);
      printf("secretscores ends at %s", asctime(localtime(&timer)));
      return 0;
   }

   ScoreEngine engine;
   engine.prefix = (const char*) prefix;
   engine.optionReadModel = (const char*) weightfile;
   engine.optionReadList = (const char*) listfile;
   engine.readFreqFile = (const char*) freqfile;
   engine.keepFile = (const char*) keepfile;
   engine.removeFile = (const char*) removefile;
   engine.CoreCount = CoreCount;
   engine.readBlockMB = blockSize;
   engine.maxBatchModels = batchModels;
   std::string genotypeName = (const char*) binfile;
   engine.genotypeFormat = GENO_BED;
   if(!vcffile.IsEmpty()){engine.genotypeFormat = GENO_VCF; genotypeName = (const char*) vcffile;}
   else if(!pfilePrefix.IsEmpty()){engine.genotypeFormat = GENO_PGEN; genotypeName = (const char*) pfilePrefix;}
   if(inferABOFlag || inferAPOEFlag || inferSexFlag){
      if(inferABOFlag) engine.InferABO(genotypeName);
      else if(inferSexFlag){engine.sexMaxSnps = sexMaxSnps; engine.InferSex(genotypeName);}
      else{engine.apoePartial = apoePartialFlag; engine.InferAPOE(genotypeName);}
      timer = time(NULL);
      printf("secretscores ends at %s", asctime(localtime(&timer)));
      return 0;
   }
   engine.PrintGeneticRiskScore(genotypeName, modelFiles, modelNames);
   return 0;
}
