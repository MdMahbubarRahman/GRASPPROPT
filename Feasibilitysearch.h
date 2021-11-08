#include <iostream>

#include "Initialsolution.h"

#ifndef FEASIBILITYSEARCH_H
#define FEASIBILITYSEARCH_H 


/*
TODO: Implement Feasibility Search algoirhtm based on "GRASP with
Path Relinking for the Two-echelon Vehicle Routing Problem" paper.
*/

//route class
class Route {
private:
	std::list<int> routeSol;
	double routeCost;
	int capacityServed;
	int sepInt;
public:
	Route();
	Route(std::list<int> routeSol,	double routeCost, int capacityServed, int sepInt);
	Route(const Route & route);
	std::list<int> getRouteSol();
	double getRouteCost();
	int getCapacityServed();
	int getSepInt();
	void showRoute();
	void updateRoute(int cus, int cap, double additionalCost);
};

//feasibility search class
class Feasibilitysearch{
private:
	ProblemParameters problemParams;
	TwoEchelonSolution currentSolution;
	TwoEchelonSolution feasibilitySearchSolution;
	bool isFeasible;
public:
	Feasibilitysearch();
	Feasibilitysearch(ProblemParameters problemParams, TwoEchelonSolution currentSolution);
	Feasibilitysearch(const Feasibilitysearch & febSearch);
	bool getFeasibilityStatus();
	TwoEchelonSolution getCurrentSolution();
	TwoEchelonSolution getFeasibilitySearchSolution();
	void showFeasibilityStatus();
	void showCurrentSolution();
	void runFeasibilitySearch();
};


#endif // !FEASIBILITYSEARCH_H




