////////////////////////////////////////////////////////////////////////
// InferAPOE.cpp
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
// --inferAPOE: the APOE genotype (e2, e3, e4 alleles) of every sample from two SNPs of the genotype data
// (--bfile, --vcf or --pfile).  See ALGORITHMS.md.
//
//   rs429358 T>C (chr19:44908684 GRCh38, 45411941 GRCh37, plus strand)   C = e4 and e1 side
//   rs7412   C>T (chr19:44908822 GRCh38, 45412079 GRCh37, plus strand)   T = e2 and e1 side
// Haplotypes (rs429358, rs7412):  e2 = T,T   e3 = T,C   e4 = C,C   (e1 = C,T, very rare).
// With c the copies of C at rs429358 and t the copies of T at rs7412 (the data are unphased):
//     c=0: t=0 e3/e3, t=1 e2/e3, t=2 e2/e2        c=1: t=0 e3/e4, t=1 e2/e4 (or e1/e3, "ambiguous"), t=2 e1/e2 ("rare")
//     c=2: t=0 e4/e4, t=1 e1/e4 ("rare"), t=2 e1/e1 ("rare")
// A SNP is found by the variant ID, else by chromosome 19 and position (GRCh38 or GRCh37, whichever finds more);
// the alleles have to be the two of the SNP or their complement.  Both SNPs are needed; a missing genotype gives NA.
// --apoe-partial: with only one of the two SNPs in the data, the copies of e2 (rs7412 T: no_e2, e2_het, e2_hom) or of e4
// (rs429358 C: no_e4, e4_het, e4_hom) are given, flagged "partial"; without it, a missing SNP is an error.
// Output: <prefix>_APOE.txt  (IID, APOE, rs429358_C, rs7412_T, E4_copies, Flag; Flag is NA when empty).
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <map>
#include "ScoreEngine.h"
#include "Error.h"

namespace {

struct ApoeSNP
{
   const char *rsid;
   char focal;         // the allele counted
   char other;
   int pos38, pos37;
};

const ApoeSNP SNPS[2] = {
   {"rs429358", 'C', 'T', 44908684, 45411941},
   {"rs7412", 'T', 'C', 44908822, 45412079}
};

char Complement(char c)
{
   switch(c){case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C';}
   return 0;
}

// the focal allele in the coding of the data, 0 when the alleles are not those of the SNP (or their complement)
unsigned int FocalAllele(const ApoeSNP & s, unsigned int a1, unsigned int a2)
{
   const unsigned int x = (unsigned char) s.focal, y = (unsigned char) s.other;
   if((a1 == x && a2 == y) || (a1 == y && a2 == x)) return x;
   const unsigned int cx = (unsigned char) Complement(s.focal), cy = (unsigned char) Complement(s.other);
   if((a1 == cx && a2 == cy) || (a1 == cy && a2 == cx)) return cx;
   return 0;
}

}

void ScoreEngine::InferAPOE(const std::string & name)
{
   wantVariantIds = true;
   printf("Loading genotype data in %s...\n", genotypeFormat == GENO_VCF? "VCF format":
      genotypeFormat == GENO_PGEN? "PLINK 2 format (pgen)": "PLINK binary format");
   LoadSamplesAndVariants(name);
   if(bimCount == 0) error("No variants were found in %s", name.c_str());
   if(bimIdHash.size() != (size_t) bimCount) error("The variant IDs were not kept");

   int idRow[2] = {-1, -1};
   std::vector<int> posRow[2][2];    // [snp][build 0 = GRCh37, 1 = GRCh38]: rows at the position with the right alleles
   for(int r = 0; r < bimCount; r++){
      for(int k = 0; k < 2; k++){
         if(bimIdHash[r] == IdHash(SNPS[k].rsid, strlen(SNPS[k].rsid)) && idRow[k] < 0 &&
            FocalAllele(SNPS[k], bimA1[r], bimA2[r]) != 0) idRow[k] = r;
         if(bimChr[r] == 19 && FocalAllele(SNPS[k], bimA1[r], bimA2[r]) != 0){
            if(bimPos[r] == SNPS[k].pos37) posRow[k][0].push_back(r);
            if(bimPos[r] == SNPS[k].pos38) posRow[k][1].push_back(r);
         }
      }
   }
   int build = 0;
   {
      int hits[2] = {0, 0};
      for(int k = 0; k < 2; k++) if(idRow[k] < 0) for(int b = 0; b < 2; b++) if(!posRow[k][b].empty()) hits[b]++;
      if(hits[0] > 0 || hits[1] > 0) build = hits[1] >= hits[0]? 38: 37;
   }
   int row[2] = {-1, -1};
   for(int k = 0; k < 2; k++){
      const char *how = "by ID";
      row[k] = idRow[k];
      if(row[k] < 0 && build != 0 && !posRow[k][build == 38? 1: 0].empty()){
         row[k] = posRow[k][build == 38? 1: 0][0];
         how = build == 38? "by position (GRCh38)": "by position (GRCh37)";
      }
      printf("  %-9s %s\n", SNPS[k].rsid, row[k] >= 0? (std::string("found ") + how).c_str(): "not in the data (by ID or by position, with the expected alleles)");
   }
   if(row[0] < 0 && row[1] < 0)
      error("--inferAPOE needs rs429358 and rs7412; neither is in the data.");
   if((row[0] < 0 || row[1] < 0) && !apoePartial)
      error("--inferAPOE needs both SNPs, rs429358 and rs7412; %s is not in the data.\n"
         "  (--apoe-partial reports what the other SNP alone tells: e2 from rs7412, or e4 from rs429358.)",
         row[0] < 0? "rs429358": "rs7412");
   if(row[0] < 0 || row[1] < 0)
      printf("  Warning: %s is not in the data: partial results (see ALGORITHMS.md); e3 and e4 cannot be told apart without rs429358.\n",
         row[0] < 0? "rs429358": "rs7412");

   std::vector<int> rows;
   for(int k = 0; k < 2; k++) if(row[k] >= 0) rows.push_back(row[k]);
   std::sort(rows.begin(), rows.end());
   rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
   threadCount = ResolveThreads(CoreCount);
   slotOfRow.assign(bimCount, -1);
   for(size_t t = 0; t < rows.size(); t++) slotOfRow[rows[t]] = (int) t;
   OpenGenotypes();
   ReadGenotypes(rows);
   CloseGenotypes();

   unsigned int focal[2];
   for(int k = 0; k < 2; k++) focal[k] = row[k] >= 0? FocalAllele(SNPS[k], bimA1[row[k]], bimA2[row[k]]): 0;

   std::map<std::string, long long> counts;
   const std::string outFile = prefix + "_APOE.txt";
   FILE *fp = fopen(outFile.c_str(), "wt");
   if(fp == NULL) error("Cannot open %s to write", outFile.c_str());
   fprintf(fp, "IID\tAPOE\trs429358_C\trs7412_T\tE4_copies\tFlag\n");
   for(int i = 0; i < idCount; i++){
      int c[2];
      for(int k = 0; k < 2; k++){
         if(row[k] < 0){c[k] = -2; continue;}                 // the SNP is not in the data (--apoe-partial)
         const int slot = slotOfRow[row[k]];
         const int code = (G_IND[i][slot >> 2] >> (2 * (slot & 3))) & 3;
         if(code == 1) c[k] = -1;
         else{
            const int copiesA1 = code == 0? 2: code == 2? 1: 0;
            c[k] = focal[k] == bimA1[row[k]]? copiesA1: 2 - copiesA1;
         }
      }
      std::string type, flag;
      int e4 = -1;
      if(c[0] == -1 || c[1] == -1) type = "NA";
      else if(c[0] == -2){      // only rs7412: the copies of e2 (e1 is extremely rare)
         type = c[1] == 0? "no_e2": c[1] == 1? "e2_het": "e2_hom";
         flag = "partial";
      }
      else if(c[1] == -2){      // only rs429358: the copies of e4 (e1 is extremely rare)
         type = c[0] == 0? "no_e4": c[0] == 1? "e4_het": "e4_hom";
         e4 = c[0];
         flag = "partial";
      }
      else{
         static const char *T[3][3] = {{"e3/e3", "e2/e3", "e2/e2"}, {"e3/e4", "e2/e4", "e1/e2"}, {"e4/e4", "e1/e4", "e1/e1"}};
         type = T[c[0]][c[1]];
         if(c[0] == 1 && c[1] == 1) flag = "ambiguous";          // e2/e4 or e1/e3 (e2/e4 is called: e1 is extremely rare)
         else if(c[0] >= 1 && c[1] >= 1 && !(c[0] == 1 && c[1] == 1)) flag = "rare";
         e4 = c[0] == 0? 0: c[0] == 1? (c[1] == 2? 0: 1): (c[1] == 0? 2: c[1] == 1? 1: 0);
      }
      counts[type]++;
      fprintf(fp, "%s\t%s\t", IID[i].c_str(), type.c_str());
      if(c[0] < 0) fprintf(fp, "NA\t"); else fprintf(fp, "%d\t", c[0]);
      if(c[1] < 0) fprintf(fp, "NA\t"); else fprintf(fp, "%d\t", c[1]);
      if(e4 < 0) fprintf(fp, "NA\t"); else fprintf(fp, "%d\t", e4);
      fprintf(fp, "%s\n", flag.empty()? "NA": flag.c_str());
      if(!flag.empty()) counts["flag " + flag]++;
   }
   fclose(fp);
   FreeGenotypes(rows);

   printf("APOE genotypes of %d samples:", idCount);
   for(std::map<std::string, long long>::const_iterator it = counts.begin(); it != counts.end(); ++it)
      if(it->first.compare(0, 5, "flag ") != 0) printf("  %s %lld", it->first.c_str(), it->second);
   printf("\n");
   if(counts["flag partial"] > 0)
      printf("  Partial results: only %s was used (e1, extremely rare, is read as e2 or e4).\n", row[0] < 0? "rs7412": "rs429358");
   if(counts["flag ambiguous"] > 0)
      printf("  %lld samples are heterozygous at both SNPs: called e2/e4 (the other reading, e1/e3, is extremely rare); marked \"ambiguous\"\n", counts["flag ambiguous"]);
   if(counts["flag rare"] > 0)
      printf("  %lld samples have a genotype with e1, which is extremely rare (check the data); marked \"rare\"\n", counts["flag rare"]);
   printf("APOE genotypes are saved in file %s\n", outFile.c_str());
}
