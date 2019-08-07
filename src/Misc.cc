#include "Misc.h"

char PrintBool( Bool b ) {
  if ( b != False )
    return 'R';
  else
    return 'F';
}

int BoolToInt( Bool b ) {
  return ( b ? 1 : 0 );
}
