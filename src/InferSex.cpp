////////////////////////////////////////////////////////////////////////
// InferSex.cpp
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
// --infersex: the sex of every sample from the X and Y chromosomes of the genotype data (--bfile, --vcf or --pfile).
// See ALGORITHMS.md.
//   X: the inbreeding coefficient F of the X chromosome outside the pseudoautosomal regions, as in PLINK --check-sex:
//      F = 1 - (observed heterozygous) / (expected heterozygous), summed over the SNPs of a sample; males have F near 1
//      (one X chromosome), females near 0.  F > 0.8 male, F < 0.2 female, else ambiguous.
//   Y: the proportion of the Y SNPs outside the pseudoautosomal regions with a genotype; the samples are split into two
//      groups at the largest gap of these proportions (when the gap is at least 0.2): the group with more genotypes is male.
//   The call is that of X when the two agree or Y is not available; a disagreement gives NA (flag X_Y_conflict).
// Output: <prefix>_sex.txt
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "ScoreEngine.h"
#include "Error.h"

namespace {

// pseudoautosomal regions, GRCh37 and GRCh38 together (the regions of the two builds do not overlap, except PAR1)
bool InParX(int pos)
{
   return (pos >= 10001 && pos <= 2781479) || (pos >= 154931044 && pos <= 155260560) || (pos >= 155701383 && pos <= 156030895);
}
bool InParY(int pos)
{
   return (pos >= 10001 && pos <= 2781479) || (pos >= 56887903 && pos <= 57217415) || (pos >= 59034050 && pos <= 59363566);
}

// at most max rows, evenly spaced
void Thin(std::vector<int> & v, int max)
{
   if(max <= 0 || (int) v.size() <= max) return;
   std::vector<int> t;
   for(int j = 0; j < max; j++) t.push_back(v[(size_t) j * v.size() / max]);
   v.swap(t);
}

}

void ScoreEngine::InferSex(const std::string & name)
{
   printf("Loading genotype data in %s...\n", genotypeFormat == GENO_VCF? "VCF format":
      genotypeFormat == GENO_PGEN? "PLINK 2 format (pgen)": "PLINK binary format");
   LoadSamplesAndVariants(name);
   if(bimCount == 0) error("No variants were found in %s", name.c_str());

   std::vector<int> xRows, yRows;
   long long xAll = 0, yAll = 0;
   for(int r = 0; r < bimCount; r++){
      if(bimChr[r] == SEXCHR){xAll++; if(bimPos[r] > 0 && !InParX(bimPos[r])) xRows.push_back(r);}
      else if(bimChr[r] == SEXCHR + 1){yAll++; if(bimPos[r] > 0 && !InParY(bimPos[r])) yRows.push_back(r);}
   }
   printf("  X chromosome SNPs: %lld (%d outside the pseudoautosomal regions); Y chromosome SNPs: %lld (%d outside)\n",
      xAll, (int) xRows.size(), yAll, (int) yRows.size());
   const bool useX = xRows.size() >= 20, useY = yRows.size() >= 20;
   if(!useX && !useY)
      error("Fewer than 20 SNPs of the X or Y chromosome outside the pseudoautosomal regions are in the data:\n"
         "  the sex cannot be inferred. (The chromosome has to be X or Y, or 23 or 24; the positions are GRCh37 or GRCh38.)");
   if(!useX) printf("  Warning: too few X SNPs: the sex is inferred from Y only.\n");
   if(!useY) printf("  Note: too few Y SNPs: the sex is inferred from X only.\n");
   Thin(xRows, sexMaxSnps);
   Thin(yRows, sexMaxSnps);

   std::vector<int> rows;
   if(useX) rows.insert(rows.end(), xRows.begin(), xRows.end());
   if(useY) rows.insert(rows.end(), yRows.begin(), yRows.end());
   threadCount = ResolveThreads(CoreCount);
   slotOfRow.assign(bimCount, -1);
   for(size_t t = 0; t < rows.size(); t++) slotOfRow[rows[t]] = (int) t;
   OpenGenotypes();
   ReadGenotypes(rows);
   CloseGenotypes();

   struct Geno{
      static int Code(unsigned char **g, int i, int slot){return (g[i][slot >> 2] >> (2 * (slot & 3))) & 3;}   // 0 A1/A1, 1 missing, 2 het, 3 A2/A2
   };
   const int n = idCount;

   // X: the F coefficient
   std::vector<double> obsHet(n, 0.0), expHet(n, 0.0), fX(n, 0.0);
   std::vector<int> usedX(n, 0);
   int xKept = 0;
   if(useX){
      for(size_t k = 0; k < xRows.size(); k++){
         const int slot = slotOfRow[xRows[k]];
         int called = 0, a1 = 0;
         for(int i = 0; i < n; i++){
            const int c = Geno::Code(G_IND, i, slot);
            if(c == 1) continue;
            called++;
            a1 += c == 0? 2: c == 2? 1: 0;
         }
         if(called < 2 || called < 0.8 * n) continue;
         const double p = a1 / (2.0 * called);
         if(p < 0.05 || p > 0.95) continue;
         xKept++;
         const double e = 2.0 * p * (1.0 - p) * (2.0 * called) / (2.0 * called - 1.0);
         for(int i = 0; i < n; i++){
            const int c = Geno::Code(G_IND, i, slot);
            if(c == 1) continue;
            usedX[i]++;
            expHet[i] += e;
            if(c == 2) obsHet[i] += 1.0;
         }
      }
      printf("  X SNPs used for F (genotyped for at least 80%% of the samples, allele frequency 0.05 to 0.95): %d\n", xKept);
   }
   // Y: the proportion of genotypes
   std::vector<double> yRate(n, 0.0);
   std::vector<int> usedY(n, 0);
   if(useY){
      for(size_t k = 0; k < yRows.size(); k++){
         const int slot = slotOfRow[yRows[k]];
         for(int i = 0; i < n; i++) if(Geno::Code(G_IND, i, slot) != 1) usedY[i]++;
      }
      for(int i = 0; i < n; i++) yRate[i] = (double) usedY[i] / yRows.size();
   }
   double yCut = -1.0;      // Y proportion above it: male; none when the samples do not split in two groups
   if(useY && n >= 2){
      std::vector<double> s(yRate);
      std::sort(s.begin(), s.end());
      double gap = 0.0;
      for(int i = 1; i < n; i++) if(s[i] - s[i - 1] > gap){gap = s[i] - s[i - 1]; yCut = (s[i] + s[i - 1]) / 2;}
      if(gap < 0.2) yCut = -1.0;
   }
   if(useY && yCut < 0) printf("  The proportions of Y genotypes do not split the samples in two groups (no gap of 0.2): Y is not used for the call.\n");

   const std::string outFile = prefix + "_sex.txt";
   FILE *fp = fopen(outFile.c_str(), "wt");
   if(fp == NULL) error("Cannot open %s to write", outFile.c_str());
   fprintf(fp, "IID\tRecordedSex\tInferredSex\tStatus\tF_X\tX_het\tX_snps\tY_rate\tY_snps\tFlag\n");
   long long nM = 0, nF = 0, nNA = 0, problems = 0, conflicts = 0, ambiguous = 0;
   for(int i = 0; i < n; i++){
      char xc = '?', yc = '?';      // M, F, A (ambiguous) or ? (none)
      double f = 0.0;
      const bool fOK = useX && usedX[i] >= 20 && expHet[i] > 0;
      if(fOK){
         f = 1.0 - obsHet[i] / expHet[i];
         xc = f > 0.8? 'M': f < 0.2? 'F': 'A';
      }
      if(useY && yCut >= 0) yc = yRate[i] > yCut? 'M': 'F';
      char call = '?';
      std::string flag;
      if(xc == 'M' || xc == 'F'){
         call = xc;
         if(yc != '?' && yc != xc){call = '?'; flag = "X_Y_conflict"; conflicts++;}
      }
      else if(yc != '?'){call = yc; flag = "Y_only";}
      if(xc == 'A'){ambiguous++; if(flag.empty() || flag == "Y_only") flag = flag.empty()? "X_ambiguous": "X_ambiguous_Y_used";}
      const char rec = sex[i] == 1? 'M': sex[i] == 2? 'F': '?';
      const char *status = (call == '?' || rec == '?')? "NA": call == rec? "OK": "PROBLEM";
      if(!strcmp(status, "PROBLEM")) problems++;
      if(call == 'M') nM++; else if(call == 'F') nF++; else nNA++;
      fprintf(fp, "%s\t%s\t%s\t%s\t", IID[i].c_str(), rec == '?'? "NA": rec == 'M'? "M": "F", call == '?'? "NA": call == 'M'? "M": "F", status);
      if(fOK) fprintf(fp, "%.4f\t%.4f\t", f, obsHet[i] / usedX[i]); else fprintf(fp, "NA\tNA\t");
      fprintf(fp, "%d\t", usedX[i]);
      if(useY) fprintf(fp, "%.4f\t%d\t", yRate[i], (int) yRows.size()); else fprintf(fp, "NA\tNA\t");
      fprintf(fp, "%s\n", flag.empty()? "NA": flag.c_str());
   }
   fclose(fp);
   FreeGenotypes(rows);

   printf("Sex of %d samples: %lld male, %lld female, %lld not inferred\n", n, nM, nF, nNA);
   if(ambiguous > 0) printf("  %lld samples have an F between 0.2 and 0.8 (ambiguous on X)\n", ambiguous);
   if(conflicts > 0) printf("  %lld samples have X and Y that disagree: not inferred (flag X_Y_conflict)\n", conflicts);
   bool anyRecorded = false;
   for(int i = 0; i < n; i++) if(sex[i] == 1 || sex[i] == 2) anyRecorded = true;
   if(!anyRecorded) printf("  No sex is recorded in the data: Status is NA.\n");
   else printf("  %lld samples have a recorded sex that differs from the inferred sex (Status PROBLEM)\n", problems);
   printf("Sex is saved in file %s\n", outFile.c_str());
}
