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

#define VERSION "1.0.0"

#include <time.h>
#include <string>
#include <vector>
#include "Parameters.h"
#include "StringArray.h"
#include "ScoreEngine.h"
#include "InputFile.h"
#include "Error.h"

void ShowBanner()
{
   printf("secretscores " VERSION " - (c) 2026-2026 Wei-Min Chen\n");
}

int main(int argc, char **argv)
{
   // line-buffered output: the messages reach the screen, a log file or a pipe as they are
   // printed, not when a buffer fills (a long run to a redirected file shows nothing otherwise)
   setvbuf(stdout, NULL, _IOLBF, 0);
   ShowBanner();

   String binfile;
   ParameterList pl;
   bool makebinFlag = false;
   String weightfile("");
   String listfile("");
   String prefix("secretscores");
   int CoreCount = 0;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Polygenic Score")
         LONG_STRINGPARAMETER("readmodel", &weightfile)
         LONG_STRINGPARAMETER("readlist", &listfile)
         LONG_PARAMETER("makebin", &makebinFlag)
      LONG_PARAMETER_GROUP("Computing Parameter")
         LONG_INTPARAMETER("cpus", &CoreCount)
      LONG_PARAMETER_GROUP("Output Parameter")
         LONG_STRINGPARAMETER("prefix", &prefix)
   END_LONG_PARAMETERS()

   pl.Add(new StringParameter('b', "Binary File", binfile));
   pl.Add(new LongParameters("Additional Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

   if(binfile.IsEmpty() && !makebinFlag)
      error("Genotype files are required. e.g.,\n  secretscores -b ex.bed --readmodel model.txt\n\n"
         "The model file is a PGS Catalog harmonized scoring file. Its header needs the columns\n"
         "  effect_allele, effect_weight, hm_chr, hm_pos\n"
         "  and other_allele or hm_inferOtherAllele; or a SECRETS binary model (file name ending in .bin,\n"
         "  made from a text model with --makebin).");
   if(weightfile.IsEmpty() && listfile.IsEmpty())
      error("Please use --readmodel <file> to specify a score model,\n  or --readlist <file> to specify a file with the names of score models.");
   if(CoreCount < 0)
      error("--cpus %d is out of range.", CoreCount);

   time_t timer;
   struct tm *tblock;
   timer = time(NULL);
   tblock = localtime(&timer);
   printf("secretscores starts at %s", asctime(tblock));

   // the models: --readmodel file1,file2,...  and  --readlist file (one model file per line,
   // with an optional second column that names the score).  A file whose name ends in .bin
   // is a SECRETS binary model, any other file is a PGS Catalog text model.
   std::vector<std::string> modelFiles, modelNames;
   if(!weightfile.IsEmpty()){
      StringArray files;
      files.AddTokens(weightfile, ',');
      for(int i = 0; i < files.Length(); i++)
         if(!files[i].IsEmpty()){
            modelFiles.push_back((const char*) files[i]);
            modelNames.push_back("");
         }
   }
   if(!listfile.IsEmpty()){
      IFILE input = ifopen((const char*)listfile, "rt");
      if(input == NULL) error("Cannot open %s to read", (const char*)listfile);
      String line;
      StringArray tokens;
      int fromList = 0;
      while(!ifeof(input)){
         tokens.Clear();
         line.ReadLine(input);
         tokens.AddTokens(line);
         if(tokens.Length() < 1 || tokens[0][0] == '#') continue;
         modelFiles.push_back((const char*) tokens[0]);
         modelNames.push_back(tokens.Length() > 1? (const char*) tokens[1]: "");
         fromList++;
      }
      ifclose(input);
      if(fromList == 0) error("No model files are listed in %s", (const char*)listfile);
   }
   if(modelFiles.empty())
      error("Please use --readmodel <file> to specify a score model,\n  or --readlist <file> to specify a file with the names of score models.");

   if(makebinFlag){
      // convert text models to binary models: no genotypes are read
      if(!binfile.IsEmpty())
         printf("--makebin only converts the models; the genotype data (-b) are not used.\n");
      MakeBinaryModels(modelFiles, modelNames);
      timer = time(NULL);
      printf("secretscores ends at %s", asctime(localtime(&timer)));
      return 0;
   }

   ScoreEngine engine;
   engine.prefix = (const char*) prefix;
   engine.optionReadModel = (const char*) weightfile;
   engine.optionReadList = (const char*) listfile;
   engine.CoreCount = CoreCount;
   engine.PrintGeneticRiskScore((const char*) binfile, modelFiles, modelNames);
   return 0;
}
