#include "MultiObjBestHammingDist.h"

#include <stdlib.h>
#include <math.h>

void MultiObjBestHammingDist::convertToMultiObjective(vector<Individual *> *population, vector<Individual *> *offSpring, vector<Individual *> *reference,  int ordObj, int numberGeneration){
	//Look for best individual (in obj0)
	Individual *best = NULL;
	if (population->size() > 0){
		best = (*population)[0];
	} else if (offSpring->size() > 0){
		best = (*offSpring)[0];
	} else {
		cerr << "Error interno: offSpring y population de tam. 0" << endl;
		exit(-1);
	}
	int direction = best->getInternalOptDirection(0);
	for (int i = 0; i < population->size(); i++){
		if ( ( (direction == MINIMIZE) && ((*population)[i]->getObj(0) < best->getObj(0)) ) ||
				 ( (direction == MAXIMIZE) && ((*population)[i]->getObj(0) > best->getObj(0))) ){
			best = (*population)[i];
		}
	}

	for (int i = 0; i < offSpring->size(); i++){
		if ( ( (direction == MINIMIZE) && ((*offSpring)[i]->getObj(0) < best->getObj(0)) ) ||
				 ( (direction == MAXIMIZE) && ((*offSpring)[i]->getObj(0) > best->getObj(0))) ){
			best = (*offSpring)[i];
		}
	}

	//Calculates the Hamming distance to the best individual
	for (int i = 0; i < population->size(); i++){
		double dist = 0;
		for (int j = 0; j < best->getNumberOfVar(); j++){
			dist += ((best->getVar(j) - (*population)[i]->getVar(j)) * (best->getVar(j) - (*population)[i]->getVar(j)));
		}
		(*population)[i]->setObj(ordObj, dist);
	}
	for (int i = 0; i < offSpring->size(); i++){
		double dist = 0;
		for (int j = 0; j < best->getNumberOfVar(); j++){
			dist += ((best->getVar(j) - (*offSpring)[i]->getVar(j)) * (best->getVar(j) - (*offSpring)[i]->getVar(j)));
		}
		(*offSpring)[i]->setObj(ordObj, dist);
	}
}
