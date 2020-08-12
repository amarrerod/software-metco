/***********************************************************************************
 * AUTHORS
 *   - Alejandro Marrero
 *   - Eduardo Segredo
 *
 * DATE
 *   August 2020
 *
 * DESCRIPTION
 *
 * Implementation of the Multi-objective Evolutionary Algorithm based on
 * Decomposition (MOEA/D) given in: Qingfu Zhang, Hui Li (2007). "MOEA/D: A
 * Multiobjective Evolutionary Algorithm Based on Decomposition". IEEE
 * Transactions on Evolutionary Computation.
 *
 * **********************************************************************************/

#ifndef __MOEAD_MPP_H__
#define __MOEAD_MPP_H__

#include <float.h>
#include <math.h>

#include <iostream>
#include <string>
#include <vector>

#include "EA.h"
#include "Individual.h"
#include "MOFront.h"

#define __MOEAD_MPP_DEBUG__

using namespace std;

class MenuPlanning;

class MOEAD_MPP : public EA {
 public:
  // Constructor
  MOEAD_MPP();

  // Destructor
  virtual ~MOEAD_MPP(void);

  // Describes one iteration of the algorithm
  void runGeneration();

  // Initialises all data structures and parameters required by the approach
  bool init(const vector<string> &params);

  // Front with the current non-dominated solutions (external population or EP
  // in the corresponding paper)
  void getSolution(MOFront *p);

  // Shows information of the algorithm
  void printInfo(ostream &os) const;

  // Getters and setters
  inline string getName() const { return "MOEA/D"; }
  inline int getNumberOfSubProblems() const { return numberSubProblems; }
  inline int getNeighbourhoodSize() const { return neighbourhoodSize; }
  inline int getNumberOfObj() const { return numberOfObj; }
  inline string getFileName() const { return fileName; }
  inline double getCrossoverRate() const { return pc; }
  inline double getMutationRate() const { return pm; }

  inline void setNumberOfSubProblems(const int &numberSubProblems) {
    this->numberSubProblems = numberSubProblems;
  }
  inline void setNeighbourhoodSize(const int &neighbourhoodSize) {
    this->neighbourhoodSize = neighbourhoodSize;
  }
  inline void setNumberOfObj(const int &numberOfObj) {
    this->numberOfObj = numberOfObj;
  }
  inline void setFileName(const string &fileName) { this->fileName = fileName; }
  inline void setCrossoverRate(const double &pc) { this->pc = pc; }
  inline void setMutationRate(const double &pm) { this->pm = pm; }

 private:
  // Crossover and mutation rates
  double pc, pm;

  // Number of subproblems, neighbourhood size, number of objective functions
  unsigned int numberSubProblems;
  unsigned int neighbourhoodSize = 10;
  unsigned int numberOfObj;

  // Weight vectors input file name
  string fileName;

  // Weight vectors
  vector<vector<double>> weightVector;

  // Neighbourhoods
  vector<vector<int>> neighbourhood;

  // Reference point
  Individual *referencePoint;

  // External population
  vector<Individual *> *secondPopulation;

  // Offspring population generated in each iteration
  vector<Individual *> offspringPop;

  // Private methods:

  // Initialises the reference point
  void initialiseReferencePoint();

  // Initialises the weight vectors
  void initialiseWeightVector();

  // Initialises the different neighbourhoods
  void initialiseNeighbourhood();

  // Initialises the population
  void initialisePopulation();

  // Generate an offsprings
  Individual *createOffspring(const int &i, const bool &);

  // Updates the reference point
  void updateReferencePoint(Individual *ind);

  // Updates the external population with non-dominated solutions
  void updateSecondPopulation(Individual *ind);

  // Compares a novel offspring to its neighboring solutions in order to update
  // the neighbourhood
  void updateNeighbouringSolution(Individual *ind, const int &i, const double &,
                                  const bool &);

  // Computes the fitness value of a particular individual by considering the
  // Tchebycheff approach
  double computingFitnessValue(Individual *ind, vector<double> &lambda,
                               const double &);

  void computeIDRanges(double &minIDS, double &maxIDS);

  void BNPSurvivalSelection(vector<Individual *> &offspring);

  // Survival Selection Technique Based on BNP - MPP CEC 2019
  void survivorSelectionBNP();

  double computeClosesDistance(const unsigned int &i);

  // Sorts neighbour weight vectors in terms of the Euclidean distance
  // between each of them and a particular weight vector in ascending order
  void minFastSort(vector<double> &dist, vector<int> &indx,
                   int numberSubProblems, int neighbourhoodSize);

  // Calculates the Euclidean distance between two weight vectors
  double getEuclideanDistance(vector<double> weightA, vector<double> weightB);

 private:
  const static int INITIAL_GENERATION;
  const static int NUM_PARAMS;
};

#endif
