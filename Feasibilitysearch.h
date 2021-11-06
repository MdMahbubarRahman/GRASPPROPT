#include <iostream>

#include "Initialsolution.h"

#ifndef FEASIBILITYSEARCH_H
#define FEASIBILITYSEARCH_H 


/*
TODO: Implement Feasibility Search algoirhtm based on "GRASP with
Path Relinking for the Two-echelon Vehicle Routing Problem" paper.
*/

class Feasibilitysearch{
private:
	ProblemParameters problemParams;
	TwoEchelonSolution currentSolution;
	bool isFeasible;
public:
	Feasibilitysearch();
	Feasibilitysearch(ProblemParameters problemParams, TwoEchelonSolution currentSolution);
	Feasibilitysearch(const Feasibilitysearch & febSearch);
	bool getFeasibilityStatus();
	TwoEchelonSolution getCurrentSolution();
	void showFeasibilityStatus();
	void showCurrentSolution();
	void runFeasibilitySearch();
};


#endif // !FEASIBILITYSEARCH_H




