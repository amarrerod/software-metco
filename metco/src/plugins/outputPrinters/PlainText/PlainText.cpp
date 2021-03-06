#include "PlainText.h"

#include <iostream>
#include <new>

#include "MOFrontVector.h"

//#define __MOEAD_MPP_ARCHIVE__

void PlainText::printSolution(EA *ga, bool end) {
  outputFile.precision(18);
  if (end) {
    ga->printInfo(outputFile);
    outputFile << endl << "Time =  " << ga->getElapsedTime() << endl << endl;
  }

  if (ga->getCritStop() == ga->getTypeStoppingCriterion("TIME")) {
    if (end) {
      outputFile << "Current Evaluation = " << ga->getPerformedEvaluations()
                 << endl;
    } else {
      outputFile << "Current Evaluation = " << ga->getNextPrint() << "("
                 << ga->getNextPrint() << ")" << endl;
    }
    if (ga->isGeneratingArchive()) {
      outputFile << ga->getLocalSolution();
    } else {
      try {
        MOFront *p = new MOFrontVector(ga->getSampleInd(), false, false);
        ga->getSolution(p);
        outputFile << p;
        delete p;
        p = nullptr;
#ifdef __MOEAD_MPP_ARCHIVE__
        p = ga->getLocalSolution();
        if (p != nullptr) {
          outputFile << p;
          delete p;
          p = nullptr;
        }
#endif
      } catch (const std::bad_alloc &except) {
        std::cerr << "Allocation failed at PlainText" << endl;
      }
    }
  } else if (ga->getCritStop() == ga->getTypeStoppingCriterion("EVALUATIONS")) {
    if (end) {
      outputFile << "Current Evaluation = " << ga->getPerformedEvaluations()
                 << endl;
    } else {
      outputFile << "Current Evaluation = " << ga->getNextPrint() << "("
                 << ga->getPerformedEvaluations() << ")" << endl;
    }
    if (ga->isGeneratingArchive()) {
      outputFile << ga->getLocalSolution();
    } else {
      try {
        MOFront *p = new MOFrontVector(ga->getSampleInd(), false, false);
        ga->getSolution(p);
        outputFile << p;
        delete (p);
      } catch (const std::bad_alloc &except) {
        std::cerr << "Allocation failed at PlainText" << endl;
      }
    }
  }
}

bool PlainText::init(const vector<string> &params) {
  if (params.size() != 1) {
    cerr << "Incorrect parameters to module PlainText: outputFile" << endl;
    return false;
  }
  outputFile.open(params[0].c_str(), ios::in | ios::out | ios::trunc);
  if (!outputFile.is_open()) {
    cerr << "No pudo abrir fichero: " << params[0] << endl;
  }
  return true;
}

void PlainText::finish() { outputFile.close(); }
