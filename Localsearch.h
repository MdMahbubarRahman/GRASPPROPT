#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <queue>

#include "Geneticalgorithm.h"


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
class Comparator {
public:
	bool operator()(CustomerDepotDifferentialCost & a, CustomerDepotDifferentialCost & b);
};


//local solution class
class LocalSolution {
private:
	int numberOfSattellites;
	int depotNode;
	std::vector<int> demands;
	std::vector<std::vector<double>> distances;//cost
	std::set<int> customerNodes;
	std::set<int> sattelliteNodes;
	std::list<int> firstEchelonRoutes;
	std::list<std::list<int>> secondEchelonRoutes;
	Chromosome firstEchelonSolution;
	std::list<Chromosome> secondEchelonSolutions;

public:
	LocalSolution();
	LocalSolution(const LocalSolution & locSol);
	LocalSolution(Chromosome firstEchelonSol, std::list<Chromosome> secondEchelonSols, std::vector<int> demands, std::vector<std::vector<double>> distances, std::vector<int> customerNodes, std::vector<int> sattelliteNodes);
	int getNumberOfSattellites();
	int getDepotNode();
	std::vector<int> getDemands();
	std::vector<std::vector<double>> getDistances();//cost
	std::set<int> getCustomerNodes();
	std::set<int> getSattelliteNodes();
	std::list<int> getFirstEchelonRoutes();
	std::list<std::list<int>> getSecondEchelonRoutes();
	Chromosome getFirstEchelonSolution();
	std::list<Chromosome> getSecondEchelonSolutions();
	void updateLocalSolution(Chromosome firstEchelonSolution, std::list<Chromosome> secondEchelonSolutions, std::list<int> firstEchelonRoutes, std::list<std::list<int>> secondEchelonRoutes);
};

//local search class
class Localsearch{
private:
	LocalSolution currentSolution;
	LocalSolution bestSolution;
	std::priority_queue<CustomerDepotDifferentialCost, std::vector<CustomerDepotDifferentialCost>, Comparator> orderedCustomerList;
public:
	Localsearch();
	Localsearch(const Localsearch & locsrch);
	Localsearch(LocalSolution currentSolution, LocalSolution bestSolution);
	LocalSolution getCurrentSolution();
	LocalSolution getBestSolution();
	void showCurrentSolution();
	void showBestSolution();
	void runLocalSearch();
	void orderCustomersBasedOnReassignmentCosts();
};


#endif // !LOCALSEARCH_H





