// Copyright 2013 jkrajniak@gmail.com (Jakub Krajniak)
// Calculates emission parameters for words and tags.
//
// Usage:
//   emission word.count test_case
//
// Returns:
//  For each word in the test_case return the word and the tag.
//

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
  fstream wordCount(argv[1], fstream::in);
  map<string, map<string, int> > wordTags; // Count(y->x)
  map<string, int> unigrams; // Count(y)

  string line, token, word, event;
  int wCount;

  while(getline(wordCount, line)) {
    stringstream lineStream(line);
    lineStream >> wCount;
    lineStream >> event;
    if (event == "WORDTAG") {
      lineStream >> token;
      lineStream >> word;
      wordTags[token][word] = wCount;
    } else if (event == "1-GRAM") {
      lineStream >> token;
      unigrams[token] = wCount;
    }
  }

  wordCount.close();
  
  map<string, map<string, double> > epsilon;
  map<string, string> bestTag;

  map<string, map<string, int> >::iterator it;
  map<string, int>::iterator it2;
  double value = 0.0;
  for (it = wordTags.begin(); it != wordTags.end(); ++it) {
    for (it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      value = (double) it2->second / (double) unigrams[it->first];
      epsilon[it2->first][it->first] = value;
    }
  }

  fstream inputDev(argv[2], fstream::in);
  map<string, double>::iterator it3;
  while(getline(inputDev, word)) {
    token = "";
    value = 0.0;
    if(word == "\n" || word == "" || word == " "){
      cout << word << endl;
    } else if (epsilon.count(word) > 0) {
      for (it3 = epsilon[word].begin(); it3 != epsilon[word].end(); ++it3) {
        if (it3->second > value) {
          token = it3->first;
          value = it3->second;
        }
      }
      cout << word << " " << token << endl;
    } else {
      for (it3 = epsilon["_RARE_"].begin(); it3 != epsilon["_RARE_"].end(); ++it3) {
        if (it3->second > value) {
          token = it3->first;
          value = it3->second;
        }
      }
      cout << word << " " << token << endl;
    }
  }

  return 0;
}
