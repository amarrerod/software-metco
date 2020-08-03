#ifndef __MPP_CHECKER_H__
#define __MPP_CHECKER_H__

#include <math.h>

#include <iostream>
#include <string>
#include <vector>

#include "EA.h"
#include "Individual.h"
#include "MOFront.h"

class MPP_Checker : public EA {
 public:
  // Constructor
  MPP_Checker(){};

  // Define una generaci�n de b�squeda del algoritmo
  void runGeneration();

  // Inicializa los par�metros iniciales del algoritmo
  bool init(const vector<string> &params);

  // Rellena un frente con las soluciones actuales
  void getSolution(MOFront *p);

  // Muestra la informaci�n relativa al algoritmo
  void printInfo(ostream &os) const;

  string getName() const { return "MPP_Checker"; }

 private:
  std::string fileToParse;

 private:
  static const int N_PARAMS;
};

#endif
