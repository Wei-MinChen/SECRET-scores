#include "Error.h"

#include "stdlib.h"
#include "stdarg.h"
#include "stdio.h"

// Declare a dummy class to ensure that compilers recognize this as C++ code
class String;

void error ( const char * msg, ... )
   {
   va_list  ap;

   va_start(ap, msg);

   // models are read by several threads: the first error prints its message and ends the
   // program; another thread that fails at the same time waits here and never prints
   #pragma omp critical(secretscores_fatal)
   {
   printf("\nFATAL ERROR - \n");
   vprintf(msg, ap);
   printf("\n\n");
   fflush(stdout);

   va_end(ap);

   exit(EXIT_FAILURE);
   }
   }

void warning ( const char * msg, ... )
   {
   va_list  ap;

   va_start(ap, msg);

   printf("\n\aWARNING - \n");
   vprintf(msg, ap);
   printf("\n");

   va_end(ap);
   }

void numerror ( const char * msg , ... )
   {
   va_list  ap;

   va_start(ap, msg);

   printf("\nFATAL NUMERIC ERROR - ");
   vprintf(msg, ap);
   printf("\n\n");

   va_end(ap);

   exit(EXIT_FAILURE);
   }
