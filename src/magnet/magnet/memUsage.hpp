/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <unistd.h>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>

namespace magnet {
  /*! \brief Attempts to read the system-dependent data for a process'
   * resident set size (actual used memory), and return the results in
   * KB.
   */
  inline double process_mem_usage()
  {
    double resident_set(0);
    {//Try the getrusage method
      ::rusage ru;
      getrusage(RUSAGE_SELF, &ru);
      resident_set = ru.ru_maxrss;
    }
    
    if (resident_set != 0) return resident_set;
  
    //Try just parsing the proc file system
    std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);
  
    if (!stat_stream.is_open()) return 0;

    // dummy vars for leading entries in stat
    std::string pid, comm, state, ppid, pgrp, session, tty_nr;
    std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    std::string utime, stime, cutime, cstime, priority, nice;
    std::string O, itrealvalue, starttime;
  
    // the two fields we want
    unsigned long vsize;
    long rss;
  
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		>> utime >> stime >> cutime >> cstime >> priority >> nice
		>> O >> itrealvalue >> starttime >> vsize >> rss;
  
    // in case x86-64 is configured to use 2MB pages
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; 
    //vm_usage     = vsize / 1024.0;
    resident_set = rss * page_size_kb;
  
    return resident_set;
  }
}
