#include "MemoryInfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(_WIN32)
  #include <windows.h>
#elif defined(__APPLE__)
  #include <mach/mach.h>
  #include <unistd.h>
#else
  #include <unistd.h>
#endif

String & MemoryInfo(double bytes)
   {
   static String info;

   if (bytes < 1024)
      return info = "<1.0 kb";

   if (bytes < 1024. * 1024.)
      info.printf("%.1f kb", (bytes + 1023) / 1024.);
   else if (bytes < 1024. * 1024. * 1024.)
      info.printf("%.1f mb", (bytes + 1024. * 1024. - 1) / (1024. * 1024.));
   else if (bytes < 1024. * 1024. * 1024. * 1024.)
      info.printf("%.1f gb", bytes / (1024. * 1024. * 1024.));
   else
      info.printf("%.1f tb", bytes / (1024. * 1024. * 1024. * 1024.));

   return info;
   }

#if !defined(_WIN32) && !defined(__APPLE__)
// the first number of a small text file (cgroup files), -1 if there is none; "max" is no limit
static double FirstNumber(const char *path)
   {
   FILE *f = fopen(path, "r");
   if (f == NULL) return -1;
   char buf[64];
   double value = -1;
   if (fgets(buf, sizeof(buf), f) != NULL && buf[0] >= '0' && buf[0] <= '9')
      value = strtod(buf, NULL);
   fclose(f);
   return value;
   }
#endif

double AvailableMemory()
   {
#if defined(_WIN32)
   MEMORYSTATUSEX status;
   status.dwLength = sizeof(status);
   if (!GlobalMemoryStatusEx(&status)) return 0;
   return (double) status.ullAvailPhys;
#elif defined(__APPLE__)
   vm_statistics64_data_t vm;
   mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
   if (host_statistics64(mach_host_self(), HOST_VM_INFO64, (host_info64_t) &vm, &count) != KERN_SUCCESS)
      return 0;
   return ((double) vm.free_count + (double) vm.inactive_count) * (double) sysconf(_SC_PAGESIZE);
#else
   double available = 0;
   // MemAvailable: what can be used without swapping (free + reclaimable cache)
   FILE *f = fopen("/proc/meminfo", "r");
   if (f != NULL)
      {
      char line[256];
      while (fgets(line, sizeof(line), f) != NULL)
         if (strncmp(line, "MemAvailable:", 13) == 0)
            {
            available = strtod(line + 13, NULL) * 1024.0;   // kB
            break;
            }
      fclose(f);
      }
   if (available <= 0)
      {
      long pages = sysconf(_SC_AVPHYS_PAGES), size = sysconf(_SC_PAGESIZE);
      if (pages > 0 && size > 0) available = (double) pages * (double) size;
      }
   // a container (cgroup v2, then v1) can allow less than the machine has
   double limit = FirstNumber("/sys/fs/cgroup/memory.max");
   double used = FirstNumber("/sys/fs/cgroup/memory.current");
   if (limit < 0)
      {
      limit = FirstNumber("/sys/fs/cgroup/memory/memory.limit_in_bytes");
      used = FirstNumber("/sys/fs/cgroup/memory/memory.usage_in_bytes");
      }
   if (limit > 0 && limit < 1e15)
      {
      double room = limit - (used > 0 ? used : 0);
      if (room < 0) room = 0;
      if (available <= 0 || room < available) available = room;
      }
   return available > 0 ? available : 0;
#endif
   }
