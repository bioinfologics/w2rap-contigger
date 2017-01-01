//
// Created by Bernardo Clavijo (TGAC) on 31/12/2016.
//

#ifndef W2RAP_CONTIGGER_OUTPUTLOG_H
#define W2RAP_CONTIGGER_OUTPUTLOG_H

#include <iostream>
#include <system/System.h>

extern unsigned OutputLogLevel;

std::ostream & OutputLog(const unsigned level,const bool include_date = true);

#endif //W2RAP_CONTIGGER_OUTPUTLOG_H
