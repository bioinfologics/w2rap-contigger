#include "dvString.h"
#include "Vec.h"
#include "ParseSet.h"
#include "MemberOf.h"


bool IsMemberOfSet(const String& set, const String& value) {
  vec<String> validValues;
  ParseStringSet(set, validValues, false);
  return Member(validValues, value);
}

bool IsMemberOfSet(const String& set, const int value) {
  vec<int> validValues;
  ParseIntSet(set, validValues, false);
  return Member(validValues, value);
}

bool IsMemberOfSet(const String& set, const longlong value) {
  vec<longlong> validValues;
  ParseLongLongSet(set, validValues, false);
  return Member(validValues, value);
}

bool IsMemberOfSet(const String& set, const unsigned int value) {
  int valueLongLong = static_cast<longlong> (value);
  return IsMemberOfSet(set, valueLongLong);
}

bool IsMemberOfSet(const String& set, const double value) {
  vec<double> validValues;
  ParseDoubleSet(set, validValues, false);
  return Member(validValues, value);
}
