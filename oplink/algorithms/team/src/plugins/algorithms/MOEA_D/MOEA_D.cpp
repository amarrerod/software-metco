/***********************************************************************************
 * AUTHORS
 *   - Alejandro Marrero
 *   - Eduardo Segredo
 *
 * DATE
 *   November 2018
 * *********************************************************************************/
#include "MOEA_D.h"
#include <float.h>
#include <random>
#include <math.h>
#include <limits>
#include <cstring>

const int MOEA_D::PARAMS = 6;
const int MOEA_D::RANDOM = 0;
const int MOEA_D::FILE = 1;
int MOEA_D::DIRECTION;

MOEA_D::MOEA_D() {

}

MOEA_D::~MOEA_D() {
    //
}

// Inicializa los parámetros iniciales del algoritmo
bool MOEA_D::init(const vector<string> &params) {
    if(params.size() < PARAMS){
        cerr << "Error\nParameters for MOEAD_D are:";
        cerr << " <population_size>";
        cerr << " <neighborhood_size>";
        cerr << " <mutationProb>";
        cerr << " <crossoverProb>";
        cerr << " <0 = MINIMIZE | 1 = MAXIMIZE>";
        cerr << " <weightGeneration> || 0 = Random || [1 = file <filename>]" << endl;
        return false;
    } else {
        setPopulationSize(stoi(params[0].c_str()));
        neighSize = stoi(params[1].c_str());
        mutationProb = stof(params[2].c_str());
        crossoverProb = stof(params[3].c_str());
        neighborhood.resize(getPopulationSize(), std::vector<int>());
        weights.resize(getPopulationSize(), std::vector<double>(getSampleInd()->getNumberOfObj(), 0));
        refPoint = unique_ptr<Individual>(getSampleInd()->internalClone());
        int direction = stoi(params[4].c_str());
        if(direction != MINIMIZE && direction != MAXIMIZE) {
            cerr << "Error\n Direction must be MINIZE or MAXIMIZE" << endl;
            return false;
        } else {
            MOEA_D::DIRECTION = direction;
        }
        if(stoi(params[5]) == FILE && params.size() == PARAMS + 1) {
            filename = params[6];
            initWeights(true);
        } else {
            initWeights(false);
        }
        initNeighborhood();
        initReferencePoint();
        updateReferencePoint((*population));
        return true;
    }
}

// Define una generación de búsqueda del algoritmo
void MOEA_D::runGeneration() {
    const int size = getPopulationSize();
    random_device rd;
    mt19937 range(rd());
    default_random_engine generator;
    uniform_int_distribution<int> intDist(0, neighSize - 1);
    uniform_real_distribution<double> realDist(0.0, 1.0);
    for(int i = 0; i < size; i++) { // Update loop
        int idxL = 0;
        int idxK = 0;
        while(idxL == idxK) {
            idxL = intDist(range);
            idxK = intDist(range);
        }
        Individual* indL = (*population)              [neighborhood[i][idxL]]->internalClone();
        Individual* indK = (*population)                        [neighborhood[i][idxK]]->internalClone();
        double genXProb = realDist(generator);
        if(genXProb < crossoverProb) {
            indL->crossover(indK);
        }
        indL->mutation(mutationProb);
        improvement(indL);                  // Keep Individual in ranges
        evaluate(indL);               // Evaluating the new Individual
        updateReferencePoint(indL);   // Updating the reference point
        updateNeighbours(i, indL); // Update Neighbors
        updateExternalPopulation(indL);
    }
}

void MOEA_D::initWeights(bool withFile) {
    const int objs = getSampleInd()->getNumberOfObj();
    const int size = getPopulationSize();
    if(withFile) {
        ifstream inputFile;
        try {
            inputFile.open(filename, ifstream::in);
            if(inputFile.is_open()) {
                int popSize;
                inputFile >> popSize;
                if (popSize != size) {
                  cerr << "Error MOEA/D: population size and file size are different" << endl;
                  inputFile.close();
                  exit(1);
                }
                for (int i = 0; i < size; i++) {
                  for (int j = 0; j < objs; j++) {
                    inputFile >> weights[i][j];
                  }
                }
                inputFile.close();
            } else {
                cerr << "Error MOEA/D: file containing weight vectors could not be opened" << endl;
                exit(1);
            }
        } catch (std::ifstream::failure e) {
            std::cerr << "Exception opening file: " << std::strerror(errno) << endl;
        }
    } else {
        default_random_engine generator;
        uniform_real_distribution<double> distribution(0.01, 1.0);
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < objs; j++){
                weights[i][j] = distribution(generator);
            }
        }
    }
}

void MOEA_D::initNeighborhood() {
    const int size = getPopulationSize();
    vector<double> distance(size);
    vector<int> idx(size);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            distance[j] = getEuclideanDistanceFromCoeff(weights[i], weights[j]);
            idx[j] = j;
        }
        minFastSort(distance, idx);
        neighborhood[i].resize(neighSize);
        for(int k = 0; k < neighSize; k++){
            neighborhood[i][k] = idx[k];
        }
    }
}

double MOEA_D::getEuclideanDistanceFromCoeff(const vector<double>& coeff1,
                                            const vector<double>& coeff2)
                                            const {
	double dist = 0;
    const int size = coeff1.size();
	for (int i = 0; i < size; i++){
		dist += (coeff1[i] - coeff2[i]) * (coeff1[i] - coeff2[i]);
	}
	return sqrt(dist);
}

void MOEA_D::minFastSort(vector<double>& distances, vector<int>& idx) {
    const int size = getPopulationSize();
    for (int i = 0; i < neighSize; i++) {
        for (int j = i + 1; j < size; j++) {
            if (distances[i] > distances[j]) {
                double temp = distances[i];
                distances[i] = distances[j];
                distances[j] = temp;
                int id = idx[i];
                idx[i] = idx[j];
                idx[j] = id;
            }
        }
    }
}

void MOEA_D::initReferencePoint() {
    const int objs = getSampleInd()->getNumberOfObj();
    for(int i = 0; i < objs; i++) {
        if(getSampleInd()->getInternalOptDirection(i) == MAXIMIZE) {
            refPoint->setVar(i, numeric_limits<double>::min());
        } else if (getSampleInd()->getInternalOptDirection(i) == MINIMIZE) {
            refPoint->setVar(i, numeric_limits<double>::max());
        }
    }
}

void MOEA_D::updateReferencePoint(Individual* ind) {
    const int objs = ind->getNumberOfObj();
    for(int i = 0; i < objs; i++) {
            if(ind->getInternalOptDirection(i) == MINIMIZE) {
                if(isless(ind->getObj(i), refPoint->getObj(i))) {
                    refPoint->setObj(i, ind->getObj(i));
                }
            } else {
                if(isgreater(ind->getObj(i), refPoint->getObj(i))) {
                    refPoint->setObj(i, ind->getObj(i));
                }
            }
    }
}

void MOEA_D::updateReferencePoint(vector<Individual*>& popu) {
    const int objs = getSampleInd()->getNumberOfObj();
    for(const Individual* ind : popu) {
        for(int j = 0; j < objs; j++) {
            if(ind->getInternalOptDirection(j) == MINIMIZE) {
                if(isless(ind->getObj(j), refPoint->getObj(j))) {
                    refPoint->setObj(j, ind->getObj(j));
                }
            } else {
                if(isgreater(ind->getObj(j), refPoint->getObj(j))) {
                    refPoint->setObj(j, ind->getObj(j));
                }
            }
        }
    }
}

void MOEA_D::updateNeighbours(const int ind, Individual* child) {
    for(int j = 0; j < neighSize; j++) {
        double childEv = decompose(child, refPoint.get(), weights.at(j));
        double neighbourdEv = decompose((*population)[neighborhood[ind][j]], refPoint.get(), weights.at(j));
        if (childEv < neighbourdEv){
            (*population)[neighborhood[ind][j]] = child->internalClone();
        }
    }

}

void MOEA_D::updateExternalPopulation(Individual* child) {
    int dominance = 0;
    const int NO_DOMINATED = 0;
    int i = 0;
    while(i < exPopulation.size()) {
        int result = dominanceTest(exPopulation[i], child);
        if(result == SECOND_DOMINATES) {
            delete(exPopulation[i]);
            exPopulation[i] = exPopulation[exPopulation.size() - 1];
            exPopulation.pop_back();
        } else if(result == FIRST_DOMINATES) {
            dominance++;
        }
        i++;
    }
    if(dominance == NO_DOMINATED) {
        exPopulation.push_back(child);
    }
}


// Mantenemos las variables dentro de los límites
void MOEA_D::improvement(Individual* ind) {
    const int vars = ind->getNumberOfVar();
    default_random_engine generator;
    for(int i = 0; i < vars; i++) {
        if(ind->getVar(i) > ind->getMaximum(i)
        || ind->getVar(i) < ind->getMinimum(i)){
            uniform_real_distribution<double> distribution(ind->getMinimum(i), ind->getMaximum(i));
            ind->setVar(i, distribution(generator));
        }
    }
}

// Rellena un frente con las soluciones actuales
void MOEA_D::getSolution(MOFront *p) {
   for (Individual* ind : exPopulation)  {
      p->insert(ind);
   }
}

// Muestra la información relativa al algoritmo
void MOEA_D::printInfo(ostream &os) const {
    os << "MOEA/D: Multiobjective Evolutionary Algorithm Based on Decomposition" << endl;
    os << "Number of Evaluations = " << getEvaluations() << endl;
    os << "Mutation Probability = " << mutationProb << endl;
    os << "Population Size = " << getPopulationSize() << endl;
    os << "Neighborhood size = " << neighSize << endl;
    os << "GLOBAL DIRECTION = " << DIRECTION << endl;
}
