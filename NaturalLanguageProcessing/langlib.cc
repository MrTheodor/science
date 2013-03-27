// Copyright 2013 jkrajniak@gmail.com (Jakub Krajniak)
// Use several groups for rare words (<5)
//  - NUMERIC for numeric values
//  - ALL_CAPITALS for all capital letters
//  - LAST_CAPITAL ends with capital letter
//  - RARE for all other cases
//

#include "langlib.h"
#include <string>

using namespace std;

string getRareToken(string token) {
  bool allCapitals = true;
  for (int i = 0; i < token.length(); i++) {
    // at least one numeric character
    if (token[i] >= 48 && token[i] <= 57) {
      return "NUMERIC";
    }
    // Check if all letters are capital.
    if (token[i] < 65 || token[i] > 90)
      allCapitals = false;
  }
  if (allCapitals)
    return "ALL_CAPITALS";
  
  char lastChar = token[token.length()-1];
  if (lastChar >= 65 && lastChar <= 90)
    return "LAST_CAPITAL";

  return "_RARE_";
}
