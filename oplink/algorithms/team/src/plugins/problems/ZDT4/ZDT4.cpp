/***********************************************************************************
 * AUTORES
 *   - Ofelia Gonz�lez P�rez
 *   - Gara Miranda Valladares
 *   - Carlos Segura Gonz�lez
 * 
 * FECHA
 *    Noviembre 2007
 ************************************************************************************/

#include "ZDT4.h"
#include <vector>
#include <time.h>
#include <sys/time.h>

// Valores minimos de cada variable
const double ZDT4::MINIMUM[] = {0.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0};
// Valores maximos de cada variable
const double ZDT4::MAXIMUM[] = {1.0,  5.0,  5.0,  5.0,  5.0,  5.0, 5.0,  5.0,  5.0,  5.0};

// Inicializaci�n. 
// Fija:
//   - Numero de variables
//   - Numero de objetivos
//   - Direccion de optimizacion de los objetivos
//   - Rangos de las variables (m�nimos y m�ximos)
bool ZDT4::init (const vector<string> &params){ 
	if (params.size() != 0){
		cerr << "ZDT4 called with parameters" << endl;
		return false;
	}
	setNumberOfVar(NPARAM);
	setNumberOfObj(NOBJ);
	return true;
}

// Evaluacion
void ZDT4::evaluate (void) {
	struct timeval ini;
	if (delay != 0)
		gettimeofday(&ini, NULL);

	setObj(0, getVar(0));
	double sum = 0;
	for (int i = 1; i < NPARAM; i++)
		sum += (getVar(i)*getVar(i) - 10.0 * cos(4.0 * M_PI * getVar(i)));
	double g = 1.0 + 10.0*(NPARAM - 1) + sum ;
	double tmp = getVar(0)/ g;
	setObj(1, g*(1.0 - sqrt(tmp)));

	if (delay != 0){
		struct timeval fin;
		double at;
		do {
			gettimeofday(&fin, NULL);
			at = (double) (fin.tv_sec) * 1.0e6 + (double) (fin.tv_usec) - (double) (ini.tv_sec) * 1.0e6 - (double) (ini.tv_usec);
		} while(at < delay);
	}
}

// Clonacion
Individual* ZDT4::clone () const {
	return new ZDT4();
}

long ZDT4::delay;
