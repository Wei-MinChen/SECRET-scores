#ifndef __MEMORYINFO_H__
#define __MEMORYINFO_H__

#include "StringBasics.h"

String & MemoryInfo(double bytes);

// RAM that is available to this process right now, in bytes (free + reclaimable
// memory, and no more than a container / cgroup limit allows); 0 when it cannot be found.
double AvailableMemory();

#endif
