////////////////////////////////////////////////////////////////////////
// InferABO.cpp
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
// --inferABO: the ABO blood type of every sample from a few SNPs of the genotype data (--bfile, --vcf or --pfile).
//
// The SNPs and the rules are those of the review of Disparities in ABO blood type determination across diverse
// ancestries (medRxiv 2024.02.21.24302372): 
//   O or not O:  rs8176719 (the deletion of a G in exon 6: allele T = O, TC = not O), the functional variant;
//                else a tag SNP (reference = not O > alternative = O): rs505922 C>T, rs657152 A>C, rs687289 A>G,
//                rs612169 G>A, rs529565 C>T, rs514659 C>A.  (rs8176704 and rs8176693 are not used: the most
//                discordant with serology in that study; rs8176645 is A/T and cannot be oriented.)
//   A or B:      rs8176746 (G = A, T = B), the functional variant; else a tag SNP (reference = A > alternative = B):
//                rs8176749 C>T, rs8176672 C>T, rs8176741 G>A, rs8176722 C>A, rs8176720 T>C.
// The SNPs are found by the variant ID (rsID) of the genotype data.  A SNP not found by its ID is looked for by
// chromosome 9 and position, in GRCh38 or GRCh37 (the build with more SNPs found; the alleles have to be those of the
// SNP, or their complement), for the SNPs whose positions are in the table (all but rs612169, rs529565, rs514659,
// rs8176672, and rs687289 in GRCh37).
// Alleles are as above or their complement (the other strand).  The SNPs are tried in the order above, for each
// sample: a missing genotype at a SNP is replaced by the next SNP.
//
// Type of a sample, with nO the copies of the O allele and a the copies of the A allele at the A/B SNP:
//     nO = 2                          O
//     nO = 0:  a = 2 A, a = 0 B, a = 1 AB
//     nO = 1:  a = 2 A, a = 0 B, a = 1 B (flagged: the O allele is assumed to be on an A background, as O alleles are
//              descended from A; the other reading is A/O with the O allele on a B background, which is rare)
//     no SNP for A/B, or no genotype there: nO < 2 gives "nonO";  no genotype for O: NA.
// The tag SNPs are much less accurate than rs8176719 and rs8176746, above all in African ancestry
// (median r2 with the deletion 0.44 in AFR, 0.87 in EUR, 0.95 in EAS).  Only ABO is inferred; the Rh type is not.
// Output: <prefix>_ABO.txt  (IID, ABO, O_copies, A_copies, Flag; Flag is NA when empty).
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <map>
#include "ScoreEngine.h"
#include "Error.h"

namespace {

struct AboSNP
{
   const char *rsid;
   char focal;         // the O allele (O SNPs) or the A allele (A/B SNPs); 0 for the deletion rs8176719
   char other;
   bool indel;
   int pos38, pos37;   // chromosome 9 positions in GRCh38 and GRCh37 (0 = not known)
};

const AboSNP O_SNPS[] = {
   {"rs8176719", 'T', 'C', true, 133257521, 136132908},
   {"rs505922", 'T', 'C', false, 133273813, 136149229}, {"rs657152", 'C', 'A', false, 133263862, 136139265},
   {"rs687289", 'G', 'A', false, 133261703, 0},
   {"rs612169", 'A', 'G', false, 0, 0}, {"rs529565", 'T', 'C', false, 0, 0}, {"rs514659", 'A', 'C', false, 0, 0}
};
const AboSNP AB_SNPS[] = {
   {"rs8176746", 'G', 'T', false, 133255935, 136131322},
   {"rs8176749", 'C', 'T', false, 133255801, 136131188}, {"rs8176672", 'C', 'T', false, 0, 0},
   {"rs8176741", 'G', 'A', false, 133256074, 136131461}, {"rs8176722", 'C', 'A', false, 133257367, 136132754},
   {"rs8176720", 'T', 'C', false, 133257486, 136132873}
};

char Complement(char c)
{
   switch(c){case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C';}
   return 0;
}

// a SNP of the data that is the SNP of the table: the focal allele in the coding of the data, 0 when not usable
unsigned int FocalAllele(const AboSNP & s, unsigned int a1, unsigned int a2)
{
   if(s.indel){
      const bool l1 = a1 < 128, l2 = a2 < 128;         // a letter, or a hash of a longer allele
      if(l1 != l2) return l1? a1: a2;                  // T against TC: the shorter allele is the deletion
      if(l1 && l2 && ((a1 == 'D' && a2 == 'I') || (a1 == 'I' && a2 == 'D'))) return 'D';
      return 0;
   }
   const unsigned int x = (unsigned char) s.focal, y = (unsigned char) s.other;
   if((a1 == x && a2 == y) || (a1 == y && a2 == x)) return x;
   const unsigned int cx = (unsigned char) Complement(s.focal), cy = (unsigned char) Complement(s.other);
   if((a1 == cx && a2 == cy) || (a1 == cy && a2 == cx)) return cx;
   return 0;
}

struct Site
{
   const AboSNP *snp;
   int row;
   unsigned int focal;
   int slot;
};

}

void ScoreEngine::InferABO(const std::string & name)
{
   wantVariantIds = true;
   printf("Loading genotype data in %s...\n", genotypeFormat == GENO_VCF? "VCF format":
      genotypeFormat == GENO_PGEN? "PLINK 2 format (pgen)": "PLINK binary format");
   LoadSamplesAndVariants(name);
   if(bimCount == 0) error("No variants were found in %s", name.c_str());
   if(bimIdHash.size() != (size_t) bimCount) error("The variant IDs were not kept");

   // the rows of the SNPs of the table
   std::map<unsigned long long, std::vector<int> > byHash;
   const int nO = (int)(sizeof(O_SNPS) / sizeof(O_SNPS[0])), nAB = (int)(sizeof(AB_SNPS) / sizeof(AB_SNPS[0]));
   for(int k = 0; k < nO; k++) byHash[IdHash(O_SNPS[k].rsid, strlen(O_SNPS[k].rsid))];
   for(int k = 0; k < nAB; k++) byHash[IdHash(AB_SNPS[k].rsid, strlen(AB_SNPS[k].rsid))];
   for(int r = 0; r < bimCount; r++){
      std::map<unsigned long long, std::vector<int> >::iterator it = byHash.find(bimIdHash[r]);
      if(it != byHash.end()) it->second.push_back(r);
   }
   // found by the ID: a row of the right alleles
   std::vector<Site> oSites, abSites;
   std::vector<int> idRow[2];   // per table, per SNP: the row found by ID (-1 none), -2 in the data with other alleles
   for(int pass = 0; pass < 2; pass++){
      const AboSNP *table = pass == 0? O_SNPS: AB_SNPS;
      const int n = pass == 0? nO: nAB;
      idRow[pass].assign(n, -1);
      for(int k = 0; k < n; k++){
         const std::vector<int> & rows = byHash[IdHash(table[k].rsid, strlen(table[k].rsid))];
         if(!rows.empty()) idRow[pass][k] = -2;
         for(size_t j = 0; j < rows.size(); j++)
            if(FocalAllele(table[k], bimA1[rows[j]], bimA2[rows[j]]) != 0){idRow[pass][k] = rows[j]; break;}
      }
   }
   // the position fallback: for the SNPs not found by ID, chromosome 9 and the position of GRCh38 or of GRCh37.
   // The build is the one with more SNPs found of the right alleles.
   int hits[2] = {0, 0};      // GRCh37, GRCh38
   std::vector<std::pair<int, int> > posRows;   // (position, row) on chromosome 9 at a position of the table
   {
      std::map<int, int> wanted;
      for(int k = 0; k < nO; k++){wanted[O_SNPS[k].pos38]; wanted[O_SNPS[k].pos37];}
      for(int k = 0; k < nAB; k++){wanted[AB_SNPS[k].pos38]; wanted[AB_SNPS[k].pos37];}
      wanted.erase(0);
      for(int r = 0; r < bimCount; r++)
         if(bimChr[r] == 9 && wanted.find(bimPos[r]) != wanted.end()) posRows.push_back(std::make_pair(bimPos[r], r));
   }
   struct ByPos{
      static int Find(const std::vector<std::pair<int, int> > & v, const AboSNP & s, int build, const std::vector<unsigned int> & a1,
         const std::vector<unsigned int> & a2)
      {
         const int pos = build == 38? s.pos38: s.pos37;
         if(pos == 0) return -1;
         for(size_t j = 0; j < v.size(); j++)
            if(v[j].first == pos && FocalAllele(s, a1[v[j].second], a2[v[j].second]) != 0) return v[j].second;
         return -1;
      }
   };
   int posBuild = 0;
   if(!posRows.empty()){
      for(int pass = 0; pass < 2; pass++){
         const AboSNP *table = pass == 0? O_SNPS: AB_SNPS;
         for(int k = 0; k < (pass == 0? nO: nAB); k++)
            if(idRow[pass][k] < 0)
               for(int b = 0; b < 2; b++)
                  if(ByPos::Find(posRows, table[k], b == 0? 37: 38, bimA1, bimA2) >= 0) hits[b]++;
      }
      if(hits[0] > 0 || hits[1] > 0) posBuild = hits[1] >= hits[0]? 38: 37;
   }
   for(int pass = 0; pass < 2; pass++){
      const AboSNP *table = pass == 0? O_SNPS: AB_SNPS;
      const int n = pass == 0? nO: nAB;
      for(int k = 0; k < n; k++){
         int row = idRow[pass][k];
         const char *how = "by ID";
         if(row < 0 && posBuild != 0){
            const int r = ByPos::Find(posRows, table[k], posBuild, bimA1, bimA2);
            if(r >= 0){row = r; how = posBuild == 38? "by position (GRCh38)": "by position (GRCh37)";}
         }
         if(row >= 0){
            Site s;
            s.snp = &table[k]; s.row = row; s.focal = FocalAllele(table[k], bimA1[row], bimA2[row]); s.slot = -1;
            (pass == 0? oSites: abSites).push_back(s);
         }
         printf("  %-10s %s: %s\n", table[k].rsid, pass == 0? "O / not O": "A / B",
            row >= 0? (std::string("found ") + how).c_str(): idRow[pass][k] == -2? "in the data, but its alleles are not the expected ones":
            (table[k].pos38 == 0 && table[k].pos37 == 0)? "not in the data (by ID; its position is not in the table)": "not in the data");
      }
   }
   if(posBuild != 0)
      printf("  Some SNPs were found by position, as GRCh%d (chromosome 9, the alleles have to match too).\n", posBuild);
   if(oSites.empty())
      error("None of the SNPs for O versus not O is in the data, by variant ID or by position (GRCh37 or GRCh38, chromosome 9,\n"
         "  with the expected alleles): rs8176719, rs505922, rs657152, rs687289, rs612169, rs529565, rs514659.");
   if(abSites.empty())
      printf("  Warning: none of the SNPs for A versus B is in the data: the samples that are not O are typed \"nonO\".\n");
   if(oSites[0].snp->indel == false)
      printf("  Warning: rs8176719, the functional variant for O, is not in the data: tag SNPs are used (see the notes in the README).\n");

   // genotypes of the sites
   std::vector<int> rows;
   for(size_t k = 0; k < oSites.size(); k++) rows.push_back(oSites[k].row);
   for(size_t k = 0; k < abSites.size(); k++) rows.push_back(abSites[k].row);
   std::sort(rows.begin(), rows.end());
   rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
   threadCount = ResolveThreads(CoreCount);
   slotOfRow.assign(bimCount, -1);
   for(size_t t = 0; t < rows.size(); t++) slotOfRow[rows[t]] = (int) t;
   OpenGenotypes();
   ReadGenotypes(rows);
   CloseGenotypes();

   // copies of the focal allele of a site: -1 when missing
   struct Copies{
      static int Of(const ScoreEngine *e, unsigned char **g, int person, const Site & s, const std::vector<unsigned int> & a1,
         const std::vector<unsigned int> & a2, int slot)
      {
         (void) e;
         const int code = (g[person][slot >> 2] >> (2 * (slot & 3))) & 3;
         if(code == 1) return -1;
         const int copiesA1 = code == 0? 2: code == 2? 1: 0;
         return s.focal == a1[s.row]? copiesA1: 2 - copiesA1;
         (void) a2;
      }
   };
   std::vector<long long> usedO(oSites.size(), 0), usedAB(abSites.size(), 0);
   std::map<std::string, long long> counts;
   long long assumed = 0;
   const std::string outFile = prefix + "_ABO.txt";
   FILE *fp = fopen(outFile.c_str(), "wt");
   if(fp == NULL) error("Cannot open %s to write", outFile.c_str());
   fprintf(fp, "IID\tABO\tO_copies\tA_copies\tFlag\n");
   for(int i = 0; i < idCount; i++){
      int o = -1, a = -1;
      for(size_t k = 0; k < oSites.size() && o < 0; k++){
         const int c = Copies::Of(this, G_IND, i, oSites[k], bimA1, bimA2, slotOfRow[oSites[k].row]);
         if(c >= 0){o = c; usedO[k]++;}
      }
      for(size_t k = 0; k < abSites.size() && a < 0; k++){
         const int c = Copies::Of(this, G_IND, i, abSites[k], bimA1, bimA2, slotOfRow[abSites[k].row]);
         if(c >= 0){a = c; usedAB[k]++;}
      }
      std::string type, flag;
      if(o < 0) type = "NA";
      else if(o == 2) type = "O";
      else if(a < 0) type = "nonO";
      else if(o == 0) type = a == 2? "A": a == 0? "B": "AB";
      else{
         type = a == 2? "A": "B";
         if(a == 1){flag = "assumed"; assumed++;}
      }
      counts[type]++;
      fprintf(fp, "%s\t%s\t", IID[i].c_str(), type.c_str());
      if(o < 0) fprintf(fp, "NA\t"); else fprintf(fp, "%d\t", o);
      if(a < 0) fprintf(fp, "NA\t"); else fprintf(fp, "%d\t", a);
      fprintf(fp, "%s\n", flag.empty()? "NA": flag.c_str());
   }
   fclose(fp);
   FreeGenotypes(rows);

   printf("ABO blood types of %d samples:", idCount);
   const char *order[] = {"A", "B", "AB", "O", "nonO", "NA"};
   for(int k = 0; k < 6; k++) if(counts[order[k]] > 0) printf("  %s %lld", order[k], counts[order[k]]);
   printf("\n");
   for(size_t k = 0; k < oSites.size(); k++)
      if(usedO[k] > 0) printf("  O / not O from %s for %lld samples\n", oSites[k].snp->rsid, usedO[k]);
   for(size_t k = 0; k < abSites.size(); k++)
      if(usedAB[k] > 0) printf("  A / B from %s for %lld samples\n", abSites[k].snp->rsid, usedAB[k]);
   if(assumed > 0)
      printf("  %lld samples with one O allele are heterozygous at the A/B SNP: called B (O assumed on the A background); marked \"assumed\"\n", assumed);
   printf("ABO types are saved in file %s\n", outFile.c_str());
}
