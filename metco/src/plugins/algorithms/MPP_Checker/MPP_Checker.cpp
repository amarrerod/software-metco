#include "MPP_Checker.h"

#include <float.h>

#include <algorithm>
#include <fstream>
#include <iterator>

const int MPP_Checker::N_PARAMS = 1;

void MPP_Checker::runGeneration() {
  std::ifstream file(fileToParse, ios::in);
  if (file.is_open()) {
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
      lines.push_back(line);
    }
    // Despues de abrir miramos los individuos finales
    std::reverse(std::begin(lines), std::end(lines));
    for (std::string line : lines) {
      if (line.rfind("Front Size", 0) == 0) {
        break;
      }
      if (line.empty()) {
        continue;
      }
      // Iteramos cada int
      Individual *ind = getSampleInd()->internalClone();
      std::istringstream iss(line);
      // Guardamos las variables
      std::vector<std::string> variables{
          std::istream_iterator<std::string>(iss), {}};
      for (int i = 0; i < variables.size() - 2; i++) {
        ind->setVar(i, stoi(variables[i]));
        std::cout << variables[i] << " ";
      }
      std::cout << std::endl;
      ind->evaluate();
    }
    exit(0);
  } else {
    std::cerr << "Error abriendo el fichero: " << fileToParse << std::endl;
    exit(-1);
  }
  file.close();
}

// Inicializa los par�metros iniciales del algoritmo
bool MPP_Checker::init(const vector<string> &params) {
  // Check for the algorithm parameters
  if (params.size() < MPP_Checker::N_PARAMS) {
    cerr << "Error MPP_Checker: Incorrect parameters (file_to_parse)" << endl;
    cerr << "Size: " << params.size() << std::endl;
    return false;
  }
  fileToParse = params[0];
  return true;
}

// Rellena un frente con las soluciones actuales
void MPP_Checker::getSolution(MOFront *p) {
  for (unsigned int i = 0; i < population->size(); i++) {
    p->insert((*population)[i]);
  }
}

// Muestra la informaci�n relativa al algoritmo
void MPP_Checker::printInfo(ostream &os) const { os << "MPP_Checker" << endl; }