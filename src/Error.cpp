#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "Error.h"

void error(const char * msg, ...)
   {
   va_list ap;
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
