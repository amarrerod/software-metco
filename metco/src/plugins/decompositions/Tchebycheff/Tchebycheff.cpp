
#include "Tchebycheff.h"

#include <cmath>
#include <cstdlib>
#include <limits>

double Tchebycheff::evaluate(Individual* ind, Individual* reference,
                             std::vector<double>& weights) {
  const int objs = ind->getNumberOfObj();
  double objective = numeric_limits<double>::min();
  double fitness = std::numeric_limits<double>::min();
  for (int i = 0; i < objs; i++) {
    double eval = (abs(ind->getObj(i) - reference->getObj(i))) * weights[i];
    if (eval > fitness) {
      fitness = eval;
    }
  }
  return fitness;
}
