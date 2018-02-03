 /***********************************************************************************
 * AUTORES
 *    Jesica de Armas Adrian 
 *    
 * FECHA
 *    Mayo 2010
 *
 *    Enero 2011
 *    Modificacion de la mutacion, para todos los elementos
 *    Gen doble tamano con relleno en 2da parte
 ************************************************************************************/

#include "TDCSP.h"
#include <vector>
#include <stack>
#include <set>
#include <algorithm>

// Direcciones de optimizacion
double* TDCSP::len;
double* TDCSP::wid;
double* TDCSP::prof;
double TDCSP::widSheet;
double TDCSP::lenSheet;

int TDCSP::totalSize;
int TDCSP::numPat;
int TDCSP::totalNumPat;
bool TDCSP::firstin;

//-------------------------------------------------------------------------------
// Funcion auxiliar para la ordenacion de pieces segun su largo
//------------------------------------------------------------------------------
bool cmpLength(const P_Order &p1, const P_Order &p2) {
                return (p1.len < p2.len);
}
//-------------------------------------------------------------------------------
// Funcion auxiliar para la ordenacion de pieces segun su ancho
//------------------------------------------------------------------------------
bool cmpWidth(const P_Order &p1, const P_Order &p2) {
                return (p1.wid < p2.wid);
}
//------------------------------------------------------------------------------
// Funcion auxiliar para la ordenacion de pieces segun su relacion prof/area
//------------------------------------------------------------------------------
bool cmpProfArea(const ProfArea &p1, const ProfArea &p2) {
                return (p1.val > p2.val);
}                
//--------------------------------------------------------------------------------
// Devuelve un entero entre 0 y maxValue - 1
// Para m�s informaci�n mirar en: man 3 rand
//--------------------------------------------------------------------------------
int TDCSP::uniformRandom(int maxValue) const {
        return (int) (rand() / ((RAND_MAX + 1.0) / maxValue));
}
//-------------------------------------------------------------------------------
// Fijar si el objetivo se debe maximizar o minimizar (indice i)
//-------------------------------------------------------------------------------
unsigned int TDCSP::getOptDirection(const int i) const{ 
	if (i == 0) 
		return MAXIMIZE;
	else	
		return MINIMIZE;     
}
// ---------------------------------------------------------------------------------
// Inicializaci�n. 
// Fija:
//   - Numero de variables
//   - Numero de objetivos
//   - Direccion de optimizacion de los objetivos
//   - Rangos de las variables (m�nimos y m�ximos)
// ---------------------------------------------------------------------------------
bool TDCSP::init (const vector<string> &params) {
	if (params.size() != 1) {
		cerr << "Error: TDCSP called with error number of parameters" << endl;
		return false;
	}
	if (!readFile(params[0])) 
		return false;
	
	setNumberOfVar(2*totalSize); // vars number 
	setNumberOfObj(NOBJ);

	return true;
}

//----------------------------------------------------------------------------------
// Reads a problem from file
// ---------------------------------------------------------------------------------
bool TDCSP::readFile(const string filename) {
	double* length;
	double* width;
	double* profit;
	int* number;
	int total;

	ifstream input(filename.c_str());
	if (!input.is_open()) {
		cerr << "Error: TDCSP called with error name of file" << endl;
		return false;
	}
	input >> numPat;
	length = new double[numPat];
  	width = new double[numPat];
	profit = new double[numPat];
	number = new int[numPat];
	

	input >> lenSheet;
	input >> widSheet;
	total = 0;
	for (int i = 0; i < numPat; i++) {
		input >> length[i];
		input >> width[i];
		input >> profit[i];
		input >> number[i];

		total += number[i];

		if (width[i] > widSheet) {
			cerr << " Error: width of pattern bigger than width of mother sheet "<< endl;
			return false;
		}
		if (length[i] > lenSheet) {
			cerr << " Error: length of pattern bigger than length of mother sheet "<< endl;
			return false;
		}
	}
	input.close();

	totalNumPat = total;
	
	totalSize = totalNumPat * 2 - 1;
	
	len = new double[totalNumPat];
	wid = new double[totalNumPat];
	prof = new double[totalNumPat];

	total = 0;
	for (int i = 0; i < numPat; i++) {
		for (int j = total; j < (total + number[i]); j++) {
			len[j] =  length[i];
			wid[j] = width[i];
			prof[j] = profit[i];
		}
		total += number[i];		
	}

	delete(length);
	delete(width);
	delete(profit);
	delete(number);

	return true;
}

//----------------------------------------------------------------------------------
// Check condition 1 <= no <= np - 1 and no = np - 1, for an operator:
// no = number of operators to its left plus itself
// np = number of patterns to its left
// ---------------------------------------------------------------------------------
bool TDCSP::checkCondition() {
	int variable;
	int numP = 0;
	int numO = 0;

	for (int i = 0; i < totalSize; i++) {
		variable = (int) getVar(i);
		if (variable != H && variable != V)
			numP++;
		else {
			numO++;
			if ((numO < 1) || (numO > (numP - 1)))
				return false;
		}
	}
	if (numO != numP - 1)
		return false;
	
	return true;
}
// -------------------------------------------------------------------------------
//  Calcula l y w 
// devuelve true cuando las dimensiones son correctas
// devuelve false cuando se sobrepasan las dimensiones
// -------------------------------------------------------------------------------
bool TDCSP::fit() {
	int var, a, b, c;
	double lengthA, lengthB, widthA, widthB;
	double  *l, *w;
	int numNewP = size/2;
	stack<int> st;
	
 	l = new double[numNewP];
   	w = new double[numNewP];
	 
	c = totalNumPat;
	
	for (int i = 0; i < size; i++) {
		var = (int) getVar(i);

		if (var != V && var != H) 
			st.push(var);	
		else {
			a = st.top();
			st.pop();
			b = st.top();
			st.pop();
			if (a >= totalNumPat) { 		// es una pieza creada por concatenacion
				lengthA = l[a-totalNumPat];
				widthA = w[a-totalNumPat];
			}
			else {			// es una pieza original
				lengthA = len[a];
				widthA = wid[a];
			}
			if (b >= totalNumPat) {		// es una pieza creada por concatenaci�n
				lengthB = l[b-totalNumPat];
				widthB = w[b-totalNumPat];                                                                                                  
			}
			else {                         		// es una pieza original                
				lengthB = len[b];
				widthB = wid[b]; 
			}
			
			if (var == H) {
				l[c - totalNumPat] = lengthA + lengthB;
				w[c - totalNumPat] = max(widthA, widthB);
			}
			else {
				l[c - totalNumPat] = max(lengthA, lengthB);
				w[c - totalNumPat] = widthA + widthB;
			}

			if (l[c - totalNumPat] > lenSheet) {
				delete(l);
				delete(w);
				return false;
			}
			if (w[c - totalNumPat] > widSheet) {
				delete(l);
				delete(w);
				return false;
			}
			st.push(c);
			c++;
		}
	}
	a = st.top();
	st.pop();
/*	
	if (a >= totalNumPat) {
		lengthA = l[a-totalNumPat];
	 	widthA = w[a-totalNumPat];
	}
	else {
		lengthA = len[a];
		widthA = wid[a];
	}

	if (lengthA > lenSheet)
		return false;
	if (widthA > widSheet)
		return false;
*/
	delete(l);
	delete(w);

	return true;
}
// -------------------------------------------------------------------------------
// Generation of an individual
// -------------------------------------------------------------------------------
void TDCSP::restart() {
	randomRestart();
}
// -------------------------------------------------------------------------------
// Random generation of an individual
// -------------------------------------------------------------------------------
void TDCSP::randomRestart() {
	int *pat, *op, j, k, numO;
	int value = uniformRandom(totalNumPat);
	set<int> st;

	pat = new int[totalNumPat]; 		// patrones
	op = new int[totalNumPat - 1];	// operadores

	profarea.resize(totalNumPat);
        // Se genera un vector de patrones ordenados aleatoriamente
        for (int i = 0; i < totalNumPat; i++) {
                while(st.find(value) != st.end()) {
                        value = uniformRandom(totalNumPat);
                }
                pat[i] = value;
                profarea[i].pat = value;
                profarea[i].val = prof[value] * (prof[value]/(len[value] * wid[value]));
                st.insert(value);
        }

        /*if (firstin) {
                // Introducir solo un individuo con profit/area ordenado de mayor a menor
                sort (profarea.begin(), profarea.end(), cmpProfArea);
                for (int i = 0; i < totalNumPat; i++) {
                        pat[i] = profarea[i].pat;
                }
                firstin = false;
        }*/

	// se genera un vector de operadores H y V 
	for (int i = 0; i < totalNumPat - 1; i++) {		
		op[i] = (int) round((rand() / (double) RAND_MAX) - 2);
	}
	// los dos primeros elementos han de ser patrones
	setVar(0, pat[0]);
	setVar(1, pat[1]);
	
	j = 2;
	k = 0;
	
	// se escoge aleatoriamente un patron o un operador (pat = 0, op = 1)
	// si es patron se agrega si aun quedan por agregar,
	// de lo contrario se agregara un operador
	// si es operador se agrega solo si cumple con la condicion 1 <= no <= np -1,
	// de lo contrario se agrega un patron
	for (int i = 2; i < totalSize; i++) {
		value = (int) round (rand()/(double) RAND_MAX);
		if (value == 0 && j < totalNumPat) {
			setVar(i,pat[j]);
			j++;
		}
		else {
			numO = k + 1;
			if ((numO < 1) || (numO > (j - 1) )) {
		 		setVar(i,pat[j]);
				j++;
			}
			else {
				setVar(i,op[k]);
				k++;
			}
		}
	}

	size = totalSize;
	delete(pat);
	delete(op);
}
//-------------------------------------------------------------------------------
// Mutacion 
// -- Se eligen 2 elementos de la cadena, p1 y p2 (donde p1 esta mas a la izquierda)
//
//		- Si ambos son piezas u operadores, se intercambian.
//		- Si p1 es un operador y p2 pieza, se intercambian.
//		- Si p1 es pieza y p2 operador, se intercambian solo cuando se cumpla:
//					1 <= no <= np - 1
//		  Se debe de seguir cumpliendo despu�s del intercambio.
//
// -- Se elige un operador aleatorio y se cambia bas�ndose en la probabilidad de mutaci�n.
//------------------------------------------------------------------------------							
void TDCSP::dependentMutation (double pm) {
	int var1, var2, aux, point1, point2;

	for (int i = 0; i < totalSize; i++) {
		point1 = i;
		double prob = rand() / (RAND_MAX + 1.0);
	        if (prob <= pm) {
			point2 = uniformRandom(totalSize);

			// Asegurar que son diferentes
			while (point1 == point2) {
				point2 = uniformRandom(totalSize);
			}

			// El primer punto ha de estar mas a la izquierda
			if (point1 > point2) {
				aux = point1;
				point1 = point2;
				point2 = aux;
			}
			var1 = (int) getVar(point1);
			var2 = (int) getVar(point2);
			
			if ( (var1 < 0 && var2 < 0) || (var1 >= 0 && var2 >= 0) || (var1 < 0 && var2 >= 0)) {
                		setVar(point1, var2);
                		setVar(point2, var1);
        		}
        		else {
                		setVar(point1, var2);
                		setVar(point2, var1);

                		if (!checkCondition()) {
                        		setVar(point1, var1);
                         		setVar(point2, var2);
                		}
        		}
		}
	}	
}
//-------------------------------------------------------------------------------
// Crossover
//-------------------------------------------------------------------------------
void TDCSP::dependentCrossover (Individual* ind) {
	//cout << "cross" << endl;
	crossoverPMX(ind);
	//cout << "sale cross" << endl;
}
//-------------------------------------------------------------------------------
// Crossover PMX: PMX Partially Mapped Crossover
//
// 1. Elegir aleatoriamente dos puntos de cruce
// 2. Intercambiar estos 2 segmentos en los hijos que se generan (solo piezas)
// 3. El resto de las cadenas se obtienen haciendo mapeos entre los 2 padres:
// 	Si un valor no est� contenido en el segmento intercambiado, permanece igual.
// 	Si est� contenido en el segmento intercambiado, entonces se sustituye por 
// 	el valor que tenga dicho segmento en el otro.
// 	
// ------------------------------------------------------------------------------
void TDCSP::crossoverPMX (Individual* ind) {
	int *p1, *p2, *h1, *h2;
	int aux, j, k, v1, v2, point1, point2;

	p1 = new int[totalNumPat];	
	h1 = new int[totalNumPat];
	p2 = new int[totalNumPat];
	h2 = new int[totalNumPat];
	j = 0;
	k = 0;

	for (int i = 0; i < totalSize; i++) {
		v1 = (int) getVar(i);
		v2 = (int) ind->getVar(i);
		if (v1 >= 0)
			p1[j++] = v1;
		if (v2 >= 0)
			p2[k++] = v2;
	}

	
	// Selecci�n de puntos de cruce
	point1 = uniformRandom(totalNumPat); 	
	point2 = uniformRandom(totalNumPat);
	
	// El primer punto ha de estar m�s a la izquierda
	if (point1 > point2) {
		aux = point1;
		point1 = point2;
		point2 = aux;
	}

	for (int i = 0; i < totalNumPat; i++) {
		// Si no pertenece al segmento de intercambio
 		if (i < point1 || i > point2) {
			// se comprueba si el numero de patron se encuentra entre los ya intercambiados
			v1 = findpat(p1[i], p1, p2, point1, point2);
			// si pertenece se cambia por el valor que corresponde
			if (v1 != -1)
				h1[i] = v1;
			// si no simplemente se copia
			else
				h1[i] = p1[i];
			v2 = findpat(p2[i], p2, p1, point1, point2);
			if (v2 != -1)
				h2[i] = v2;
			else
				h2[i] = p2[i];

		}
		// Si pertenece al segmento de intercambio se intercambian
		else {
			h1[i] = p2[i];
			h2[i] = p1[i];
		}
	}

	// reasignaci�n a los genes correspondientes
	j = 0;
	k = 0;
	for (int i = 0; i < totalSize; i++) {
		if (getVar(i) > -1)
			setVar(i, h1[j++]);
		if (ind->getVar(i) > -1)
			ind->setVar(i,h2[k++]); 
	}


	
	//repair();
	//((TDCSP*) ind)->repair();

	delete(p1);
	delete(p2);
	delete(h1);
	delete(h2);
}

//------------------------------------------------------------------------------
// B�squeda de un patron en el segmento de intercambio, devolviendo el nuevo valor
// si ya se encuentra en el segmento intercambiado.
// Es recursiva ya que puede ocurrir que el nuevo valor tambi�n se encuentre en
// el segmento de intercambio.
// Used by crossoverPMX.
//------------------------------------------------------------------------------
int TDCSP::findpat(int elem, int *p1, int *p2, int point1, int point2) {
        int newelem;
        for(int i = point1; i <= point2; i++){
                if(p2[i] == elem) {
                        elem = p1[i];
                        newelem = findpat(elem, p1, p2, point1, point2);
                        if (newelem != -1)
                                elem = newelem;
                        return elem;
                }
        }
        return -1;
}

//-------------------------------------------------------------------------------
// Evaluacion
// ------------------------------------------------------------------------------
void TDCSP::evaluate (void) {
	//cout << "entro evaluate" << endl;
	int var, a, b, c, numOp, numPie, finalSize, finalNc;
	double lengthA, lengthB, widthA, widthB;
	double  *l, *w;
	int nc = 0;
	int numNewP = totalSize/2;      //????????????
	stack<int> st;
	double finalProfit, total_profit = 0;	

 	l = new double[numNewP];
   	w = new double[numNewP];
	 
	c = totalNumPat;

	finalSize = 1;
	numOp = numPie = 0;
	finalProfit = 0;
	finalNc = 0;
	bool change = false;

	for (int i = 0; i < totalSize; i++) {
		var = (int) getVar(i);
	
		if (var != V && var != H) {
			st.push(var);
			numPie++;
			if (i == 0)
				finalProfit = prof[var];
			total_profit += prof[var];
		}
		else {
			numOp++;
			nc++;
			b = st.top();
			st.pop();
			a = st.top();
			st.pop();

			if (a >= totalNumPat) { 		// es una pieza creada por concatenacion
				lengthA = l[a-totalNumPat];
				widthA = w[a-totalNumPat];
			}
			else {				// es una pieza original
				lengthA = len[a];
				widthA = wid[a];
			}
			if (b >= totalNumPat) {		// es una pieza creada por concatenaci�n
				lengthB = l[b-totalNumPat];
				widthB = w[b-totalNumPat];                                                                                                  
			}
			else {                         // es una pieza original                
				lengthB = len[b];
				widthB = wid[b]; 
			}
			
			if (var == H) {
				l[c - totalNumPat] = lengthA + lengthB;
				w[c - totalNumPat] = max(widthA, widthB);
				if (widthA != widthB)
		                        nc++;
				if (l[c - totalNumPat] > lenSheet) {                                                                                               
                                	if (change) {
						size = finalSize;      
						break;
					}
					else {                  
                                        	setVar(i, V);
						// cambiar H por V                                                                                        
                                        	i--;                                                     
                                      		st.push(a);
						st.push(b);                                                                                      
                                                change = true;                                                                                                                            
                                                nc--;                                                   
                                                if (widthA != widthB)
                                                	nc--;                                                                                                                                                                     numOp--;                                                                                                                                  
                                                continue;
                                	}                                                                                                                                                                         }                                                                                                                                                         
                        	change = false;
			}
			else {
				l[c - totalNumPat] = max(lengthA, lengthB);
				w[c - totalNumPat] = widthA + widthB;
				if (lengthA != lengthB)
		                        nc++;
				 if (w[c - totalNumPat] > widSheet) {
                                 	if (change) {
                                        	size = finalSize;
                                        	break;
                                 	}
                                 	else {
                                        	setVar(i, H);                   // cambiar V por H
                                        	i--;
                                        	st.push(a);
                                        	st.push(b);
                                      		change = true;
                                        	nc--;
                                        	if (lengthA !=  lengthB)
                                                	nc--;
                                        	numOp--;
                                        	continue;
                                	}
	                        }
        	                change = false;


			}
			
			if (numOp == numPie - 1) {
				finalSize = 2 * numPie  - 1;
				finalProfit = total_profit;
				finalNc = nc;
			}

			st.push(c);
			c++;
		}
	}

	while (st.size() != 0) {
		a = st.top();
		st.pop();
	}

	if (a >= totalNumPat) {
		lengthA = l[a-totalNumPat];
	 	widthA = w[a-totalNumPat];
	}
	else {
		lengthA = len[a];
		widthA = wid[a];
	}


	setObj(0, finalProfit);
	
	if (widthA != widSheet) {
		//cout << "anchos diferentes" << endl;
		finalNc += 1;
	}
	if (lengthA != lenSheet)
		//cout << "largos diferentes" << endl;
		//cout << lengthA << " " << widthA << endl;
		finalNc += 1;
	
	setObj(1, finalNc);


	// Hacer copia en la segunda parte hasta size
	for (int i = 0; i < size; i++) {
		setVar(totalSize + i, (int) getVar(i));
	}
	
	// Aplicar relleno, que cambie size!!!
	fill();
	
	// Actualizar cortes y profit

	delete(l);
	delete(w);	
	l = new double[numNewP];
        w = new double[numNewP];

        c = totalNumPat;
	total_profit  = 0;
	nc = 0;

        for (int i = totalSize; i < size + totalSize; i++) {
                var = (int) getVar(i);

                if (var != V && var != H) {
                        st.push(var);
                        total_profit += prof[var];
                }
                else {
                        nc++;
                        a = st.top();
                        st.pop();
                        b = st.top();
                        st.pop();

                        if (a >= totalNumPat) {                 // es una pieza creada por concatenacion
                                lengthA = l[a-totalNumPat];
                                widthA = w[a-totalNumPat];
                        }
                        else {                          // es una pieza original
                                lengthA = len[a];
                                widthA = wid[a];
                        }
                        if (b >= totalNumPat) {         // es una pieza creada por concatenacisn
                                lengthB = l[b-totalNumPat];
                                widthB = w[b-totalNumPat];
                        }
                        else {                         // es una pieza original
                                lengthB = len[b];
                                widthB = wid[b];
                        }

                        if (var == H) {
                                l[c - totalNumPat] = lengthA + lengthB;
                                w[c - totalNumPat] = max(widthA, widthB);
                                if (widthA != widthB)
                                        nc++;
                        }
                        else {
                                l[c - totalNumPat] = max(lengthA, lengthB);
                                w[c - totalNumPat] = widthA + widthB;
                                if (lengthA != lengthB)
                                        nc++;
                        }

                        st.push(c);
                        c++;
                }
        }
        a = st.top();
        st.pop();

        if (a >= totalNumPat) {
                lengthA = l[a-totalNumPat];
                widthA = w[a-totalNumPat];
        }
        else {
                lengthA = len[a];
                widthA = wid[a];
        }


        setObj(0, total_profit);

        if (widthA != widSheet)
                nc += 1;
        if (lengthA != lenSheet)
                nc += 1;

        setObj(1, nc);
	

	//print();
	delete(l);
	delete(w);
	//cout << "salgo evaluate" << endl;
}
//-------------------------------------------------------------------------------
// Apila en vertical u horizontal las piezas sobrantes que se van a colocar de relleno
//-------------------------------------------------------------------------------
void TDCSP::pilUp(vector<P_Order> &piecesOrd, int it, int op, int maxDim) {
        int k = 0;
        int count = 0;
        int ltotal = 0;
        int wtotal = 0;
        //cout << "entra pill" << endl;
        while (k <= it) { // se deberia de salir antes de que se cumpla la condicion
                count = 0;
                ltotal = 0;
                wtotal = 0;

                for (int j = it; j >= k; j--) {
                        setVar(size + totalSize, piecesOrd[j].pat);
                        size++;
                        count++;
                        ltotal += piecesOrd[j].len;
                        wtotal += piecesOrd[j].wid;
                }
                for (int j = 1; j < count; j++) {
                        if (op == H)
                                setVar(size + totalSize, V);
                        else
                                setVar(size + totalSize, H);
                        size++;
                }
                setVar(size + totalSize, op);

                size++;

                if ((op == H && maxDim >= wtotal) || (op == V && maxDim >= ltotal)) {
                        piecesOrd.erase(piecesOrd.begin()+k, piecesOrd.begin()+it+1);
                        break;
                }
                else {
                        size -= 2*count;
                        k++;
                }
        }
        //cout << "sale pill" << endl;
}
// -------------------------------------------------------------------------------
//  Calcula l y w sobrantes
//  CASO 0: Devuelve el largo sobrante y a = W si (dimension = 0)
//              o el ancho sobrante y L-lsobrante si (dimension = 1)
//  CASO 1: Devuelve el largo sobrante y ancho de la contruccion si (dimension 0)
//              o el ancho sobrante y largo total (dimension 1)
// -------------------------------------------------------------------------------
void TDCSP::getRemainSize(int dimension, int build, int &lgap, int & wgap) {
        int var, a, b, c;
        double lengthA, lengthB, widthA, widthB;
        double  *l, *w;
        int numNewP = size/2;
        stack<int> st;

        //cout << "entra getR" << endl;
        l = new double[numNewP];
        w = new double[numNewP];

        c = totalNumPat;

        for (int i = totalSize; i < size + totalSize; i++) {
                var = (int) getVar(i);

                if (var != V && var != H)
                        st.push(var);
                else {
                        a = st.top();
                        st.pop();
                        b = st.top();
                        st.pop();
                        if (a >= totalNumPat) {                 // es una pieza creada por concatenacion
                                lengthA = l[a-totalNumPat];
                                widthA = w[a-totalNumPat];
                        }
                        else {                  // es una pieza original
                                lengthA = len[a];
                                widthA = wid[a];
                        }
                        if (b >= totalNumPat) {         // es una pieza creada por concatenacisn
                                lengthB = l[b-totalNumPat];
                                widthB = w[b-totalNumPat];
                        }
                        else {                                  // es una pieza original
                                lengthB = len[b];
                                widthB = wid[b];
                        }

                        if (var == H) {
                                l[c - totalNumPat] = lengthA + lengthB;
                                w[c - totalNumPat] = max(widthA, widthB);
                        }
                        else {
                                l[c - totalNumPat] = max(lengthA, lengthB);
                                w[c - totalNumPat] = widthA + widthB;
                        }

                        if (l[c - totalNumPat] > lenSheet) {
                                delete(l);
                                delete(w);
                                //return 0;
                                cout << "error"<< endl;
                        }
                        if (w[c - totalNumPat] > widSheet) {
                                delete(l);
                                delete(w);
                                cout << "error"<< endl;
                                //return 0;
                        }
                        st.push(c);
                        c++;
                }
        }
        a = st.top();
        st.pop();

        if (a >= totalNumPat) {
                lengthA = l[a-totalNumPat];
                widthA = w[a-totalNumPat];
        }
        else {
                lengthA = len[a];
                widthA = wid[a];
        }

        delete(l);
        delete(w);

        if (dimension == 0) {
                lgap = lenSheet - lengthA;
                if (build == 0)
                        wgap = widSheet;
                else
                        wgap = widthA;
        }
        else {
                wgap = widSheet - widthA;
                if (build == 0)
                        lgap = lengthA;
                else
                        lgap = lenSheet;
        }
        //cout << "sale getR" << endl;
}
//-----------------------------------------------------------------------------------
// Operador de Relleno
//-----------------------------------------------------------------------------------
void TDCSP::fill(){
        int v1;
        set <int> pieces;
        vector <P_Order> piecesOrd;
        vector <P_Order> aux;

        for (int i = 0; i < size; i++) {
                v1 = (int) getVar(i);
                if (v1 >= 0)
                        pieces.insert(v1);
        }
        int f = 0;
        // Relleno con las piezas sobrantes
        piecesOrd.resize(totalNumPat - pieces.size());
        for (int i = 0; i < totalNumPat; i++) {
                if (pieces.find(i) == pieces.end()) {
                        piecesOrd[f].pat = i;
                        piecesOrd[f].len = len[i];
                        piecesOrd[f].wid = wid[i];
                        f++;
                }
        }

        // fill = 0     p x V y H       fill = 1        p y H x V
        //  ------------                 ------------
        // |  x  |      |               |  x         |
        // |-----    y  |               |----- ------|
        // |  p  |      |               |  p  |  y   |
        //  ------------                 ------------

        int lgap, wgap, it = 0;
        int fill = uniformRandom(2);
        if (fill == 0) {
                // Las ordena de menor a mayor ancho
                sort (piecesOrd.begin(), piecesOrd.end(), cmpWidth);
                lgap = 0;
                wgap = 0;
                it = 0;

                while (!piecesOrd.empty()) {
                        getRemainSize(1, 0, lgap, wgap);
                        // si cabe la pieza de ancho
                        if (wgap >= piecesOrd[it].wid) {
                                // si no cabe de largo
                                if (lgap < piecesOrd[it].len) {
                                        aux.push_back(piecesOrd[it]);
                                        piecesOrd.erase(piecesOrd.begin()+it);
                                        if (it == piecesOrd.size())
                                                it--;
                                        // si la siguiente no cabe de ancho no va a apilar, hacerlo ahora
                                        if (it > 0 && it < piecesOrd.size()) {
                                                if (wgap < piecesOrd[it].wid){
                                                         pilUp(piecesOrd, it-1, V, lgap);
                                                         it = 0;
                                                }
                                        }
                                }
                                else {
                                        // si cabe de largo
                                        if (piecesOrd.size() > 1 && it < piecesOrd.size()-1) {
                                                // si no cabe la siguiente pieza
                                                if (wgap < piecesOrd[it+1].wid) {
                                                        pilUp(piecesOrd,it,V, lgap);
                                                        it = 0;
                                                }
                                                else // si cabe la siguiente pieza
                                                        it++;
                                        }
                                        else {
                                                // si es la unica pieza
                                                if (piecesOrd.size() == 1) {
                                                        setVar(size + totalSize, piecesOrd[it].pat);
                                                        setVar(size + totalSize + 1, V);
                                                        size += 2;
                                                        piecesOrd.pop_back();
                                                }
                                                else {
                                                        // si it esta en el ultimo lugar
                                                        pilUp(piecesOrd,it,V,lgap);
                                                        it = 0;
                                                }
                                        }
                                }
                        }
                        else
                                break;
                }

                piecesOrd.insert(piecesOrd.begin(), aux.begin(), aux.end());
                aux.clear();
                if (!piecesOrd.empty()) {
                        // Las ordeno de menor a mayor largo
                        sort (piecesOrd.begin(), piecesOrd.end(), cmpLength);
                        lgap = 0;
                        wgap = 0;
                        it = 0;

                        while (!piecesOrd.empty()) {
                                // si cabe la pieza
                                getRemainSize(0, 0, lgap, wgap);
                                if (lgap >= piecesOrd[it].len) {
                                        if (wgap < piecesOrd[it].wid) {
                                                aux.push_back(piecesOrd[it]);
                                                piecesOrd.erase(piecesOrd.begin()+it);
                                                if (it == piecesOrd.size())
                                                        it--;
                                                // si la siguiente no cabe de ancho no va a apilar, hacerlo ahora
                                                if (it > 0 && it < piecesOrd.size()) {
                                                        if (lgap < piecesOrd[it].len){
                                                                pilUp(piecesOrd,it-1,H, wgap);
                                                                it = 0;
                                                        }
                                                }
                                        }
                                        else {
                                                if (piecesOrd.size() > 1 && it < piecesOrd.size()-1) {
                                                        // si no cabe la siguiente pieza
                                                        if (lgap < piecesOrd[it+1].len) {
                                                                pilUp(piecesOrd,it,H,wgap);
                                                                it = 0;
                                                        }
                                                        else // si cabe la siguiente pieza
                                                                it++;
                                                }
                                                else {
                                                        if (piecesOrd.size() == 1) {
                                                                setVar(size + totalSize, piecesOrd[it].pat);
                                                                setVar(size + totalSize + 1, H);
                                                                size += 2;
                                                                piecesOrd.pop_back();
                                                        }
                                                        else { // si it esta en el ultimo lugar,no hay mas piezas despues
                                                                pilUp(piecesOrd,it,H,wgap);
                                                                it = 0;
                                                        }
                                                }
                                        }
                                }
                                else
                                        break;
                        }
                }
        }
        else {
                // Las ordeno de menor a mayor largo
                sort (piecesOrd.begin(), piecesOrd.end(), cmpLength);
                lgap = 0;
                wgap = 0;
                it = 0;

                while (!piecesOrd.empty()) {
                        // si cabe la pieza
                        getRemainSize(0, 1, lgap, wgap);
                        if (lgap >= piecesOrd[it].len) {
                                if (wgap < piecesOrd[it].wid) {
                                        aux.push_back(piecesOrd[it]);
                                        piecesOrd.erase(piecesOrd.begin()+it);
                                        if (it == piecesOrd.size())
                                                it--;
                                        // si la siguiente no cabe de ancho no va a apilar, hacerlo ahora
                                        if (it > 0 && it < piecesOrd.size()) {
                                                if (lgap < piecesOrd[it].len){
                                                        pilUp(piecesOrd,it-1,H, wgap);
                                                        it = 0;
                                                }
                                        }
                                }
                                else {
                                        if (piecesOrd.size() > 1 && it < piecesOrd.size()-1) {
                                                // si no cabe la siguiente pieza
                                                if (lgap < piecesOrd[it+1].len) {
                                                        pilUp(piecesOrd,it,H,wgap);
                                                        it = 0;
                                                }
                                                else // si cabe la siguiente pieza
                                                        it++;
                                        }
                                        else {
                                                if (piecesOrd.size() == 1) {
                                                        setVar(size + totalSize, piecesOrd[it].pat);
                                                        setVar(size + totalSize + 1, H);
                                                        size += 2;
                                                        piecesOrd.pop_back();
                                                }
                                                else { // si it esta en el ultimo lugar,no hay mas piezas despues
                                                        pilUp(piecesOrd,it,H,wgap);
                                                        it = 0;
                                                }
                                        }
                                }
                        }
                        else
                                break;
                }

                piecesOrd.insert(piecesOrd.begin(), aux.begin(), aux.end());
                aux.clear();
                if (!piecesOrd.empty()) {
                        sort (piecesOrd.begin(), piecesOrd.end(), cmpWidth);
                        lgap = 0;
                        wgap = 0;
                        it = 0;

                        while (!piecesOrd.empty()) {
                                // si cabe la pieza
                                getRemainSize(1, 1, lgap, wgap);
                                if (wgap >= piecesOrd[it].wid) {
                                         // si no cabe de largo
                                        if (lgap < piecesOrd[it].len) {
                                                aux.push_back(piecesOrd[it]);
                                                piecesOrd.erase(piecesOrd.begin()+it);
                                                if (it == piecesOrd.size())
                                                        it--;
                                                // si la siguiente no cabe de ancho no va a apilar, hacerlo ahora
                                                if (it > 0 && it < piecesOrd.size()) {
                                                        if (wgap < piecesOrd[it].wid){
                                                                pilUp(piecesOrd,it-1,V, lgap);
                                                                it = 0;
                                                        }
                                                }
                                        }
                                        else {
                                                if (piecesOrd.size() > 1 && it < piecesOrd.size()-1) {
                                                        // si no cabe la siguiente pieza
                                                        if (wgap < piecesOrd[it+1].wid) {
                                                                pilUp(piecesOrd,it,V, lgap);
                                                                it = 0;
                                                        }
                                                        else // si cabe la siguiente pieza
                                                                it++;
                                                }
                                                else {
                                                        if (piecesOrd.size() == 1) {
                                                                setVar(size + totalSize, piecesOrd[it].pat);
                                                                setVar(size + totalSize + 1, V);
                                                                size += 2;
                                                                piecesOrd.pop_back();
                                                        }
                                                        else { // si it estC! en el ultimo lugar
                                                                pilUp(piecesOrd,it,V,lgap);
                                                                it = 0;
                                                        }
                                                }
                                        }
                                }
                                else
                                        break;
                        }
                }
        }
        aux.clear();
        pieces.clear();
        piecesOrd.clear();
}
//----------------------------------------------------------------------------
// Clonacion
// ---------------------------------------------------------------------------
Individual* TDCSP::clone (void) const {
  	TDCSP* ind = new TDCSP();
	ind->size = size;
	return ind;
}
//------------------------------------------------------------------------------
// Imprimir por pantalla un gen o individuo
//-----------------------------------------------------------------------------
void TDCSP::print(void) {
	unsigned int i;
		  
	cout << "Gene: ";
	for (i = 0; i < size; i++) {
		if (getVar(i) == -1)
			cout << "H ";
		else if (getVar(i) == -2)
			cout << "V ";
		else
			cout << getVar(i) << " ";
	}
			 
  	cout << endl << "Profit: " << getObj(0) << endl;	
  	cout << "Number of cuts: " << getObj(1) << endl << endl;
	
}
/*
//------------------------------------------------------------------------------
// Codigo auxiliar para iniciar una determinada solucion
//------------------------------------------------------------------------------
void TDCSP::fixSolution(void) {
	// Fijar el gen y su correspondiente orientacion
	
	// Test Problem 4 - Solution 1
	//string gene("7 8 10 3 H H H 2 4 9 H 5 H 1 6 H H H V");
	//string orientation("0 0 1 1 0 1 0 1 0 0");
	// Tests Problem 4 - Solution 2
	//string gene("5 7 2 H 9 3 H H H 1 6 H 8 4 H H 10 H V");
	//string orientation("0 0 0 1 0 0 1 1 1 1");	
	
	// Tests Problem 5 - Solution 1
	//string gene("5 13 2 H 15 9 1 H H 8 H 18 16 H H 3 H H H 12 7 17 14 19 H 4 20 H H H H H 11 10 H H 6 H V");
	//string orientation("0 0 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0");	
	// Tests Problem 5 - Solution 2
	//string gene("18 8 5 H H 11 19 H 6 20 H 7 4 H H H 17 10 H H H 13 2 16 H 3 1 H 9 H 15 14 H 12 H H H H V");
	//string orientation("0 0 1 1 0 1 1 1 0 0 0 1 1 0 1 0 0 1 0 1");	
	// Tests Problem 5 - Solution 3
	//string gene("13 6 9 H 20 17 H 16 15 H 8 H H H 3 H H 10 19 2 H H 18 H 1 H 12 H 5 H 7 H 11 4 H H 14 H V");
	//string orientation("0 1 1 1 0 1 1 1 1 0 0 1 0 0 0 0 0 1 0");	
	// Tests Problem 5 - Solution 4
	//string gene("5 12 19 10 17 9 H H 2 7 16 H H 15 H H H H H 3 1 18 6 H 13 H 11 20 H H H H 14 H 4 8 H H V");
	//string orientation("0 0 1 0 0 1 1 1 0 1 0 1 1 1 0 0 1 1 0 0");	
	// Tests Problem 5 - Solution 5
	//string gene("12 9 H 4 13 2 20 17 H 3 H H H H 14 19 H H H 11 7 H 5 H 18 H 6 H 16 15 8 H 1 H H 10 H H V");
	//string orientation("0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 0 1 1 0 1");	

	//string gene("7 8 9 H H 20 2 5 V H V 1 19 14 H 4 V 10 V 3 V H 17 12 18 H V 13 V H 15 11 H 6 16 H V H H");
	//string orientation("0 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 1 0 1");
	
	string gene("17 14 5 19 9 H V 13 V H 3 18 20 H V 6 7 V 1 H 11 H 2 15 H V 4 10 H V H 12 16 V H H H H");
	string orientation("0 0 0 0 1 0 1 1 1 1 1 0 0 0 1 0 0 1 1 0");
	stringstream geneStream(gene); 
	stringstream orientationStream(orientation); 
	string elem; 
	unsigned int i = 0;

	// Almacenar el gen
	while (geneStream >> elem) {
		if (elem == "H") 
			setVar(i, -1); 
		else if (elem == "V") 
			setVar(i, -2); 
		else 
			setVar(i, (atoi(elem.c_str())) - 1); 
		i++;
	}
	
	while (orientationStream >> elem) {
		setVar(i, atoi(elem.c_str())); 
		i++;
	}

	// Evalua y muestra la solucion
	evaluate();
	print();
	exit(-1);
}

*/
