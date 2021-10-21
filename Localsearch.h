#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <queue>

#include "Geneticalgorithm.h"
#include "Initialsolution.h"


#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H 


/*
TODO: Implement Local Search algorithm presented in "GRASP with Path Relinking
for the Two-Echelon Vehicle Routing Problem" paper.

More insight of this algorithm can be found from the following paper:
	1. Clustering-based heuristics for the two-echelon capacitated vehicle routing problem by
	   TG Crainic, S Mancini, G Perboli, R Tadei
*/

//customer to depots potential assignment class
class CustomerDepotDifferentialCost {
private:
	int customerID;
	int currentDepot;
	int potentialDepot;
	double differentialCost;
public:
	CustomerDepotDifferentialCost();
	CustomerDepotDifferentialCost(const CustomerDepotDifferentialCost & cusDifCost);
	CustomerDepotDifferentialCost(int cusID, int curDepot, int potDepot, double diffCost); 
	int getCustomerID();
	int getCurrentDepot();
	int getPotentialDepot();
	double getDifferentialCost();
	void showCustomerDepotDifferentialCost();
};

//comparator for customer to depot distance
class Compare{
public:
	bool operator()(CustomerDepotDifferentialCost & a, CustomerDepotDifferentialCost & b);
};

//local search class
class Localsearch{
private:
	TwoEchelonSolution currentSolution;
	TwoEchelonSolution bestSolution;
	std::priority_queue<CustomerDepotDifferentialCost, std::vector<CustomerDepotDifferentialCost>, Compare> orderedCustomerList;
	CVRPSolution firstEchelonSolution;
	CVRPSolution currentSatelliteSolution;
	CVRPSolution potentialSatelliteSolution;
	bool isNodeAdditionFeasible;//needs justification?
public:
	Localsearch();
	Localsearch(const Localsearch & locsrch);
	Localsearch(TwoEchelonSolution currentSolution, TwoEchelonSolution bestSolution);
	TwoEchelonSolution getCurrentSolution();
	TwoEchelonSolution getBestSolution();
	void showCurrentSolution();
	void showBestSolution();
	void runLocalSearch();
	void customersReassignmentOrder();
	void checkNodeAdditionFeasibility();
	void addNodeToThePotentialSatellite();
	void deleteNodeFromTheCurrentSatellite();
	void updateFirstEchelonSolution();
	void computeObjectiveValueOfTheNewSolution();
};


#endif // !LOCALSEARCH_H

