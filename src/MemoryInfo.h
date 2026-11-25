#ifndef __MEMORYINFO_H__
#define __MEMORYINFO_H__

#include <string>

// a size in bytes as text: "<1.0 kb", "12.5 mb", "3.2 gb"
std::string MemoryInfo(double bytes);

// RAM that is available to this process right now, in bytes (free + reclaimable
// memory, and no more than a container / cgroup limit allows); 0 when it cannot be found.
double AvailableMemory();

#endif
