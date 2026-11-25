////////////////////////////////////////////////////////////////////////
// Parameters.h
//
// Command line options of the form --name [value].
//
//   OptionString idDelim("_"); bool quiet = false; int cpus = 0; double maf = -1.0;
//   BEGIN_LONG_PARAMETERS(options)
//      LONG_PARAMETER_GROUP("Input")
//         LONG_STRINGPARAMETER("id-delim", &idDelim)
//      LONG_PARAMETER_GROUP("Output")
//         LONG_PARAMETER("quiet", &quiet)
//         LONG_INTPARAMETER("cpus", &cpus)
//         LONG_DOUBLEPARAMETER("maf", &maf)
//   END_LONG_PARAMETERS()
//   ParameterList pl;
//   pl.Add(new LongParameters("Options", options));
//   pl.Read(argc, argv);
//   pl.Status();
//
//      SHORT_STRINGPARAMETER("b", &binfile)      gives  -b <value> or -b<value>
//
// Option names are not case sensitive and may be shortened to any prefix that
// is unique (an exact name always wins over longer names that start with it).
// A value option needs a value, and a number must be a number: otherwise the
// program stops with an error.  Unknown or ambiguous options and loose words
// are reported as warnings.
////////////////////////////////////////////////////////////////////////
#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <string>
#include <vector>

// The value of a string option: the text, read as a C string
class OptionString
   {
   public:
      OptionString() {}
      OptionString(const char *text) : text_(text ? text : "") {}
      OptionString &operator=(const char *text) { text_ = text ? text : ""; return *this; }
      operator const char *() const { return text_.c_str(); }
      bool IsEmpty() const { return text_.empty(); }
      int Length() const { return (int) text_.size(); }
   private:
      std::string text_;
   };

enum OptionKind { OPTION_END, OPTION_GROUP, OPTION_FLAG, OPTION_STRING, OPTION_INT, OPTION_DOUBLE, OPTION_SHORT_STRING };

struct OptionDef
   {
   OptionKind kind;
   const char *name;       // option name without the dashes (for a group: its title)
   void *var;              // bool, OptionString, int or double (an OptionString for a short option)
   };

#define BEGIN_LONG_PARAMETERS(list)         OptionDef list[] = {
#define LONG_PARAMETER_GROUP(title)         { OPTION_GROUP, title, 0 },
#define LONG_PARAMETER(name, var)           { OPTION_FLAG, name, var },
#define LONG_STRINGPARAMETER(name, var)     { OPTION_STRING, name, var },
#define LONG_INTPARAMETER(name, var)        { OPTION_INT, name, var },
#define LONG_DOUBLEPARAMETER(name, var)     { OPTION_DOUBLE, name, var },
// a one-letter option with a value, written -b <value> or -b<value> (the letter is case sensitive)
#define SHORT_STRINGPARAMETER(letter, var)  { OPTION_SHORT_STRING, letter, var },
#define END_LONG_PARAMETERS()               { OPTION_END, 0, 0 } };

class LongParameters
   {
   public:
      LongParameters(const char *title, OptionDef *defs) : title_(title), defs_(defs) {}
      const std::string &Title() const { return title_; }
      OptionDef *Defs() const { return defs_; }
   private:
      std::string title_;
      OptionDef *defs_;
   };

class ParameterList
   {
   public:
      ~ParameterList();
      void Add(LongParameters *p) { lists_.push_back(p); }       // the list takes ownership
      void Read(int argc, char **argv);                          // sets the variables
      void Status();                                             // prints the settings and the warnings to the screen
      const std::string &Warnings() const { return warnings_; }  // problems found by Read (empty = none)
   private:
      std::vector<LongParameters *> lists_;
      std::string warnings_;
   };

#endif
