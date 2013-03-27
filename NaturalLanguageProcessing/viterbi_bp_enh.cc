// Copyright 2013 jkrajniak@gmail.com (Jakub Krajniak)
// Draft implementation of Viterbi algorithm with backpointers
// and the four categories for rare word.
//
// Usage:
//    viterbi_bp_enh word.count test_case
//
// Returns:
//  For each word in the test_case return the word and the tag.
//

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <set>
#include <deque>
#include <algorithm>
#include <cstdio>
#include "langlib.h"

using namespace std;

int main(int argc, char **argv) {
  fstream wordCount(argv[1], fstream::in);
  map<string, map<string, int> > wordTags;  // Count(y->x)
  map<string, int> unigrams;
  map<string, map<string, int> > bigrams;
  map<string, map<string, map<string, int> > > trigrams;
  set<string> tags;
  tags.insert("I-GENE");
  tags.insert("O");

  string line, token, token1, token2, token3, word, event;
  int wCount;

  while (getline(wordCount, line)) {
    stringstream lineStream(line);
    lineStream >> wCount;
    lineStream >> event;
    if (event == "WORDTAG") {
      lineStream >> token;
      lineStream >> word;
      wordTags[token][word] = wCount;
      //cout << "\t" << line << endl;
      //cout << "wordTags[" << token << "][" << word << "] " << wCount << endl;
    } else if (event == "1-GRAM") {
      lineStream >> token;
      unigrams[token] = wCount;
    } else if (event == "2-GRAM") {
      lineStream >> token;
      lineStream >> token1;
      bigrams[token][token1] = wCount;
    } else if (event == "3-GRAM") {
      lineStream >> token;
      lineStream >> token1;
      lineStream >> token2;
      trigrams[token][token1][token2] = wCount;
    }
  }

  wordCount.close();

  // Calculates emission parameters.
  map<string, map<string, double> > epsilon;

  map<string, map<string, int> >::iterator it;
  map<string, int>::iterator it2;
  double value = 0.0;
  for (it = wordTags.begin(); it != wordTags.end(); ++it) {
    for (it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      value = static_cast<double>(it2->second) /
        static_cast<double>(unigrams[it->first]);
      epsilon[it2->first][it->first] = value;
      //cout << "e[" << it2->first << "][" << it->first << "] " << value << endl;
    }
  }

  // Calculates q parameters.
  map<string, map<string, map<string, double> > > q;
  map<string, map<string, map<string, double> > >::iterator qIt;
  map<string, map<string, map<string, int> > >::iterator triIt;

  string yn2, yn1, yn;
  for (triIt = trigrams.begin(); triIt != trigrams.end(); triIt++) {
    yn2 = triIt->first;
    for (it = triIt->second.begin(); it != triIt->second.end(); it++) {
      yn1 = it->first;
      for (it2 = it->second.begin(); it2 != it->second.end(); it2++) {
        yn = it2->first;
        double c = it2->second;
        q[yn][yn2][yn1] = c / static_cast<double>(bigrams[yn2][yn1]);
      }
    }
  }

  // Reads input sentences.
  vector< vector<string>* > sentences;
  vector<string> *tmpSentence = new vector<string>;
  fstream inputDev(argv[2], fstream::in);
  while (getline(inputDev, word)) {
    if (word == "") {
      sentences.push_back(tmpSentence);
      tmpSentence = new vector<string>;
    } else {
      tmpSentence->push_back(word);
    }
  }

  // Calculates pi and backpointers
  // Viterbi algorithm implementation
  map<int, map<string, map<string, double> > > pi;  // pi(k, u, v);
  map<int, map<string, map<string, string> > > bp;  // bp(k, u, v);
  vector<string> resultsTag;

  set<string> *uSet = NULL;
  set<string> *vSet = NULL;
  set<string> *wSet = NULL;
  set<string> singleTag;
  singleTag.insert("*");

  set<string>::iterator uSetIt, vSetIt, wSetIt;
  string maxToken, xk;
  double maxValue;

  map<string, double> *emissionParam;
  double piScale = 1000.0;

  pi[0]["*"]["*"] = piScale*1.0;

  for (int sIdx = 0; sIdx < sentences.size(); sIdx++) {
    vector<string> *sentence = sentences[sIdx];
    resultsTag.clear();
    resultsTag.resize(sentence->size()+1);
    fill(resultsTag.begin(), resultsTag.begin()+sentence->size(), "");
    // Iterates over sentence.
    for (int k = 1; k <= sentence->size(); k++) {
      // Word.
      xk = sentence->at(k-1);
      // Define sets.
      uSet = (k < 2) ? &singleTag : &tags;
      vSet = &tags;
      wSet = (k < 3) ? &singleTag : &tags;

      // cout << "----DEBUG k=" << k << ": \"" << xk << "\" --------" << endl;
      // Emission parameter.
      if (epsilon.count(xk) > 0) {
        emissionParam = &epsilon[xk];
      } else { // RARE
        string rareToken = getRareToken(xk);
        //cout << rareToken << " for " << xk << endl;
        emissionParam = &epsilon[rareToken];
      }
      double emissionValue = 0.0;
      for (uSetIt = uSet->begin(); uSetIt != uSet->end(); uSetIt++) {
        for (vSetIt = vSet->begin(); vSetIt != vSet->end(); vSetIt++) {
          wSetIt = wSet->begin();
          maxToken = "O";
          if (emissionParam->count(*vSetIt) > 0)
            emissionValue = emissionParam->at(*vSetIt);
          else
            emissionValue = 0.0;
          maxValue = piScale*pi[k-1][maxToken][*uSetIt]
            * q[*vSetIt][maxToken][*uSetIt] * emissionValue;
          value = 0.0;
          for (; wSetIt != wSet->end(); wSetIt++) {
             value = piScale*pi[k-1][*wSetIt][*uSetIt]
               * q[*vSetIt][*wSetIt][*uSetIt]
               * emissionValue;
             if (value > maxValue) {
               maxValue = value;
               maxToken = *wSetIt;
             }
             
             /*
             printf("### pi cal for (%d, %s, %s) ###\n", k, (*uSetIt).c_str(), (*vSetIt).c_str());
             printf("pi[%d][%s][%s] = %e\n", k-1, (*wSetIt).c_str(), (*uSetIt).c_str(), (double) pi[k-1][*wSetIt][*uSetIt]);
             printf("q[%s][%s][%s] = %e\n", (*vSetIt).c_str(), (*wSetIt).c_str(), (*uSetIt).c_str(), q[*vSetIt][*wSetIt][*uSetIt]);
             printf("em[%s][%s] = %e\n", xk.c_str(), (*vSetIt).c_str(), emissionValue);
             printf("maxValue: %e, value: %e, maxToken: %s\n\n", maxValue, value, maxToken.c_str()); 
             */
             
          }
          pi[k][*uSetIt][*vSetIt] = maxValue;
          bp[k][*uSetIt][*vSetIt] = maxToken;
        }
      }
    }
    // Backtrace for sentence.
    // Fill last two tags.
    maxValue = 0.0;
    yn1 = "";  // tag n-1
    yn = "";  // last tag
    int n = sentence->size();
    for (uSetIt = tags.begin(); uSetIt != tags.end(); uSetIt++) {
      for (vSetIt = tags.begin(); vSetIt != tags.end(); vSetIt++) {
        value = pi[n][*uSetIt][*vSetIt] * q["STOP"][*uSetIt][*vSetIt];
        if (value > maxValue) {
          maxValue = value;
          yn = *vSetIt;
          yn1 = *uSetIt;
        }
      }
    }
    resultsTag[n] = yn;
    resultsTag[n-1] = yn1;
    // Fill the rest tags.
    for (int k = n-2; k >= 1; k--) {
      resultsTag[k] = bp[k+2][resultsTag[k+1]][resultsTag[k+2]];
    }

    for (int k = 0; k < n; ++k) {
      cout << sentence->at(k) << " " << resultsTag[k+1] << endl;
    }
    cout << "\n";
  }  // Finished of processing sentence

  return 0;
}
