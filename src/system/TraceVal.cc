#include <stdlib.h>
#include <iostream>
#include "system/TraceVal.h"

TraceValCommon::timestamp_t TraceValCommon::nextTimeStamp_ = 1;
Bool TraceValCommon::tracingCopiesStopped_ = False;
TraceValCommon::timestamp_t TraceValCommon::stopAtStamp_ = 0;

template class TraceVal<int>;
