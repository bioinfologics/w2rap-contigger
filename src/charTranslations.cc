// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "charTranslations.h"
#include "system/Crash.h"

int carraySz(int strl) {
  return strl/3 + (strl%3 ? 1 : 0) + 1;
}

int textToCharSimple(char *s, unsigned char *carr) {

  int i, ptr=0, cl = strlen(s);
  for (i=0; i<cl-2; i+=3)
    carr[ptr++] = 25*char2num[(int)s[i]] + 5*char2num[(int)s[i+1]] + char2num[(int)s[i+2]];

  switch (cl - i) {
  case 2:
    carr[ptr++] = 125 + 5*char2num[(int)s[i]] + char2num[(int)s[i+1]];
    carr[ptr++] = 255;
    return ptr;
  case 1:
    carr[ptr++] = 150 + char2num[(int)s[i]];
    carr[ptr++] = 255;
    return ptr;
  case 0:
    carr[ptr++] = 255;
    return ptr;
  default:
    Assert(false);
  }
Assert( false ) ;
return ptr;
}

int textToCharCompact(const char *s, unsigned char *carr) {

  int i, ptr=0, cl = strlen(s);
  for (i=0; i<cl-2; i+=3) {
    carr[ptr++] = simple2compact[ 25*char2num[(int)s[i]] + 5*char2num[(int)s[i+1]] + char2num[(int)s[i+2]] ];
  }

  switch (cl - i) {
  case 2:
    carr[ptr++] = simple2compact[ 125 + 5*char2num[(int)s[i]] + char2num[(int)s[i+1]] ];
    carr[ptr++] = 255;
    return ptr;
  case 1:
    carr[ptr++] = simple2compact[ 150 + char2num[(int)s[i]] ];
    carr[ptr++] = 255;
    return ptr;
  case 0:
    carr[ptr++] = 255;
    return ptr;
  default:
    Assert(false);
  }
Assert( false ) ;
return ptr;
}


int textToCharSimple(int textl, char *s, unsigned char *carr) {

  int i, ptr=0;
  for (i=0; i<textl-2; i+=3)
    carr[ptr++] = 25*char2num[(int)s[i]] + 5*char2num[(int)s[i+1]] + char2num[(int)s[i+2]];

  switch (textl - i) {
  case 2:
    carr[ptr++] = 125 + 5*char2num[(int)s[i]] + char2num[(int)s[i+1]];
    carr[ptr++] = 255;
    return ptr;
  case 1:
    carr[ptr++] = 150 + char2num[(int)s[i]];
    carr[ptr++] = 255;
    return ptr;
  case 0:
    carr[ptr++] = 255;
    return ptr;
  default:
    Assert(false);
  }
Assert( false ) ;
return ptr;
}

int textToCharCompact(int textl, char *s, unsigned char *carr) {

  int i, ptr=0;
  for (i=0; i<textl-2; i+=3) {
    carr[ptr++] = simple2compact[ 25*char2num[(int)s[i]] + 5*char2num[(int)s[i+1]] + char2num[(int)s[i+2]] ];
  }

  switch (textl - i) {
  case 2:
    carr[ptr++] = simple2compact[ 125 + 5*char2num[(int)s[i]] + char2num[(int)s[i+1]] ];
    carr[ptr++] = 255;
    return ptr;
  case 1:
    carr[ptr++] = simple2compact[ 150 + char2num[(int)s[i]] ];
    carr[ptr++] = 255;
    return ptr;
  case 0:
    carr[ptr++] = 255;
    return ptr;
  default:
    Assert(false);
  }
Assert( false ) ;
return ptr;
}

int charCompactToText(unsigned char *carr, char *ctext) {

  int ptr = 0, txtptr = 0;
  const char *tempbuf;
  while (carr[ptr] != EOCH) {
    tempbuf = char2string_compact[carr[ptr++]];
    ctext[txtptr++] = tempbuf[0];
    ctext[txtptr++] = tempbuf[1];
    ctext[txtptr++] = tempbuf[2];
  }
  ctext[txtptr++] = '\0';
  return txtptr;
}
int charSimpleToText(unsigned char *carr, char *ctext) {

  int ptr = 0, txtptr = 0;
  const char *tempbuf;
  while (carr[ptr] != EOCH) {
    tempbuf = char2string_simple[carr[ptr++]];
    ctext[txtptr++] = tempbuf[0];
    ctext[txtptr++] = tempbuf[1];
    ctext[txtptr++] = tempbuf[2];
  }
  ctext[txtptr++] = '\0';
  return txtptr;
}
