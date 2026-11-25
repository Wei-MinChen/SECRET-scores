////////////////////////////////////////////////////////////////////////
// Parameters.cpp
////////////////////////////////////////////////////////////////////////
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "Parameters.h"
#include "Error.h"

namespace {

std::string Lower(const char *s)
   {
   std::string r(s);
   for (size_t i = 0; i < r.size(); i++)
      r[i] = (char) tolower((unsigned char) r[i]);
   return r;
   }

// "--name": an option word (a lone "--", and words with one dash such as -0.5, are not)
bool IsOptionWord(const char *s)
   {
   return s[0] == '-' && s[1] == '-' && s[2] != 0;
   }

// "-b" or "-bvalue": one dash and a letter
bool IsShortWord(const char *s)
   {
   return s[0] == '-' && isalpha((unsigned char) s[1]);
   }

// the value of a number option must be the whole word
bool ParseInt(const char *s, int &v)
   {
   char *end = NULL;
   errno = 0;
   long x = strtol(s, &end, 10);
   if (end == s || *end != 0 || errno == ERANGE || x > 2147483647L || x < -2147483647L)
      return false;
   v = (int) x;
   return true;
   }

bool ParseDouble(const char *s, double &v)
   {
   char *end = NULL;
   errno = 0;
   double x = strtod(s, &end);
   if (end == s || *end != 0 || errno == ERANGE)
      return false;
   v = x;
   return true;
   }

// how the settings table shows a number: 0.50, 1.0e-03
std::string ShowDouble(double v)
   {
   char buf[64];
   if (v == 0.0 || v >= 0.01)
      snprintf(buf, sizeof(buf), "%.2f", v);
   else
      snprintf(buf, sizeof(buf), "%.1e", v);
   return buf;
   }

// "--name [value]" for the settings table; a switch shows [ON] when it is set, a zero integer shows nothing
std::string ShowOption(const OptionDef &d)
   {
   std::string s = std::string("--") + d.name;
   char buf[64];
   switch (d.kind)
      {
      case OPTION_FLAG:
         if (*(bool *) d.var)
            s += " [ON]";
         break;
      case OPTION_SHORT_STRING:
         s = std::string("-") + d.name;
         s += std::string(" [") + (const char *) *(OptionString *) d.var + "]";
         break;
      case OPTION_STRING:
         s += std::string(" [") + (const char *) *(OptionString *) d.var + "]";
         break;
      case OPTION_INT:
         if (*(int *) d.var != 0)
            {
            snprintf(buf, sizeof(buf), " [%d]", *(int *) d.var);
            s += buf;
            }
         break;
      case OPTION_DOUBLE:
         s += " [" + ShowDouble(*(double *) d.var) + "]";
         break;
      default:
         break;
      }
   return s;
   }

const int kLabelWidth = 22;      // the group titles are right-aligned in this width
const int kLineWidth = 79;       // the settings table wraps before this many characters

}

ParameterList::~ParameterList()
   {
   for (size_t i = 0; i < lists_.size(); i++)
      delete lists_[i];
   }

void ParameterList::Read(int argc, char **argv)
   {
   warnings_.clear();
   // all options of all lists
   std::vector<const OptionDef *> all;
   for (size_t k = 0; k < lists_.size(); k++)
      for (const OptionDef *d = lists_[k]->Defs(); d->kind != OPTION_END; d++)
         if (d->kind != OPTION_GROUP)
            all.push_back(d);

   for (int i = 1; i < argc; i++)
      {
      const char *word = argv[i];
      char note[512];
      if (IsShortWord(word))
         {
         const OptionDef *shortHit = NULL;
         for (size_t k = 0; k < all.size(); k++)
            if (all[k]->kind == OPTION_SHORT_STRING && all[k]->name[0] == word[1])
               shortHit = all[k];
         if (shortHit == NULL)
            {
            snprintf(note, sizeof(note), "Command line parameter %s is undefined\n", word);
            warnings_ += note;
            continue;
            }
         if (word[2] != 0)                                   // -bvalue
            {
            *(OptionString *) shortHit->var = word + 2;
            continue;
            }
         if (i + 1 >= argc || IsOptionWord(argv[i + 1]))
            error("-%s needs a value.", shortHit->name);
         *(OptionString *) shortHit->var = argv[++i];
         continue;
         }
      if (!IsOptionWord(word))
         {
         snprintf(note, sizeof(note), "Command line parameter %s (#%d) ignored\n", word, i);
         warnings_ += note;
         continue;
         }
      // an exact name wins; otherwise the only option that starts with the word
      const std::string want = Lower(word + 2);
      const OptionDef *hit = NULL;
      int prefixHits = 0;
      for (size_t k = 0; k < all.size(); k++)
         {
         if (all[k]->kind == OPTION_SHORT_STRING) continue;
         const std::string name = Lower(all[k]->name);
         if (name == want)
            {
            hit = all[k];
            prefixHits = 1;
            break;
            }
         if (name.compare(0, want.size(), want) == 0)
            {
            hit = all[k];
            prefixHits++;
            }
         }
      if (prefixHits == 0)
         {
         snprintf(note, sizeof(note), "Command line parameter %s is undefined\n", word);
         warnings_ += note;
         continue;
         }
      if (prefixHits > 1)
         {
         snprintf(note, sizeof(note), "Command line parameter %s is ambiguous\n", word);
         warnings_ += note;
         continue;
         }
      if (hit->kind == OPTION_FLAG)
         {
         *(bool *) hit->var = true;
         continue;
         }
      // the other options take the next word as their value
      if (i + 1 >= argc || IsOptionWord(argv[i + 1]))
         error("--%s needs a value.", hit->name);
      const char *value = argv[++i];
      if (hit->kind == OPTION_STRING)
         *(OptionString *) hit->var = value;
      else if (hit->kind == OPTION_INT)
         {
         int v = 0;
         if (!ParseInt(value, v))
            error("--%s needs an integer, but \"%s\" was given.", hit->name, value);
         *(int *) hit->var = v;
         }
      else
         {
         double v = 0.0;
         if (!ParseDouble(value, v))
            error("--%s needs a number, but \"%s\" was given.", hit->name, value);
         *(double *) hit->var = v;
         }
      }
   }

void ParameterList::Status()
   {
   printf("\nThe following parameters are in effect:\n");
   for (size_t k = 0; k < lists_.size(); k++)
      {
      printf("\n%s\n", lists_[k]->Title().c_str());
      const OptionDef *d = lists_[k]->Defs();
      while (d->kind != OPTION_END)
         {
         // one group: its title, then its options separated by commas, wrapped
         const char *title = d->name;
         std::vector<std::string> items;
         for (d++; d->kind != OPTION_END && d->kind != OPTION_GROUP; d++)
            items.push_back(ShowOption(*d));
         std::string line;
         char head[64];
         snprintf(head, sizeof(head), "%*s :", kLabelWidth, title);
         line = head;
         const std::string blank(line.size(), ' ');
         for (size_t n = 0; n < items.size(); n++)
            {
            const std::string item = items[n] + (n + 1 < items.size() ? "," : "");
            if (line.size() > blank.size() && line.size() + 1 + item.size() > (size_t) kLineWidth)
               {
               printf("%s\n", line.c_str());
               line = blank;
               }
            line += " " + item;
            }
         printf("%s\n", line.c_str());
         }
      }
   printf("\n");
   if (!warnings_.empty())
      printf("\nWARNING - \nProblems encountered parsing command line:\n\n%s\n", warnings_.c_str());
   }
