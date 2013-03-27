// Copyright 2013 jkrajniak@gmail.com (Jakub Krajniak)
//
// Calculates the frequency of words.
// Use several groups for rare words (<5)
//  - NUMERIC for numeric values
//  - ALL_CAPITALS for all capital letters
//  - LAST_CAPITAL ends with capital letter
//  - RARE for all other cases
//

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cstdio>
#include <sstream>
#include <vector>
#include "langlib.h"

using namespace std;

int main(int argc, char **argv) {
  int rareCount = 5;

  fstream wordsFile(argv[1], fstream::in);
  map<string, int> wordFreq;
  
  string line;
  string token;
  string tag;

  while (wordsFile.good()) {
    getline(wordsFile, line);
    stringstream lineStream(line);
    getline(lineStream, token, ' ');
    getline(lineStream, tag, ' ');
    if (wordFreq.count(token) == 0)
      wordFreq[token] = 1;
    else
      wordFreq[token]++;
  }
  wordsFile.clear();
  wordsFile.seekg(0);

  while(wordsFile.good()) {
    getline(wordsFile, line);
    stringstream lineStream(line);
    getline(lineStream, token, ' ');
    getline(lineStream, tag, ' ');
    if (wordFreq[token] < rareCount) {
      cout << getRareToken(token) << " " << tag << endl; 
    } else {
      cout << line << endl;
    }
  }
  wordsFile.close();
  return 0;
}
