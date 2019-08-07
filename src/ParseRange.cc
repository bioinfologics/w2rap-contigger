#include "dvString.h"
#include "ParseRange.h"

template<> int ConvertFromString<int>(const String& value) {
  return value.Int();
}

template<> unsigned int ConvertFromString<unsigned int>(const String& value) {
  return value.Int();
}

template<> long ConvertFromString<long>(const String& value) {
  return value.Int(); // Int() returns a long
}

template<> double ConvertFromString<double>(const String& value) {
  return value.Double();
}

template<> String ConvertFromString<String>(const String& value) {
  return value;
}
