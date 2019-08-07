/*
 * \file HostName.cc
 * \author tsharpe
 * \date Jan 9, 2012
 *
 * \brief
 */
#include "system/HostName.h"
#include "system/SysConf.h"
#include "system/System.h"

std::string getHostName()
{
    size_t len = maxHostNameLen()+1; // 1 extra for null

    char* buf = new char[len+1]; // 1 extra beyond that
    if ( gethostname(buf,len) == -1 )
        FatalErr("Unable to determine host name.");
    buf[len] = 0; // make completely, absolutely sure we're null terminated

    std::string result(buf);
    delete [] buf;
    return result;
}
