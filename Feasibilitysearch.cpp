#include "Feasibilitysearch.h"

//Implement Feasibility search algorithm here.

//default constructor
Feasibilitysearch::Feasibilitysearch() {
	isFeasible = false;
}

//constructor
Feasibilitysearch::Feasibilitysearch(ProblemParameters probParams, TwoEchelonSolution currSol) {
	problemParams = probParams;
	currentSolution = currSol;
}

//copy constructor
Feasibilitysearch::Feasibilitysearch(const Feasibilitysearch& febSearch) {
	problemParams = febSearch.problemParams;
	currentSolution = febSearch.currentSolution;
	isFeasible = febSearch.isFeasible;
}

//returns feasibility status
bool Feasibilitysearch::getFeasibilityStatus() {
	return isFeasible;
}

//returns current solution
TwoEchelonSolution Feasibilitysearch::getCurrentSolution() {
	return currentSolution;
}

//prints feasibility status
void Feasibilitysearch::showFeasibilityStatus() {
	std::cout << "Feasibility status : " << isFeasible << std::endl;
}

//prints current solution
void Feasibilitysearch::showCurrentSolution() {
	currentSolution.showTwoEchelonSolution();
}

//runs feasibility search algorithm
void Feasibilitysearch::runFeasibilitySearch() {
	//needs implementation

}


