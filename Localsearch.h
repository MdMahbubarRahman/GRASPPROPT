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

//class cvrp solution 
class CVRPSolution {
private:
	std::vector<int> solution;
	std::set<int> customers;
	std::map<int, int> customerTodemand;
	double solCost = 0;
	int satelliteNode = 0;
	int maxRouteCapacity = 0;
	int totalDemandSatisfied = 0;
	int numberOfRoutes = 0;
public:
	CVRPSolution();//default constructor
	CVRPSolution(std::vector<int> sol, std::set<int> customers, std::map<int, int> customerTodemand, double cost, int satelliteNode, int maxRouteCapacity, int totalDemandSatisfied, int numberOfRoutes);//constructor
	CVRPSolution(const CVRPSolution & sol); //copy constructor
	void showSolution() const;
	void showCustomers();
	void showCustomerTodemandMap();
	void showSolCost();
	void showSatelliteNode();
	void showMaxRouteCapacity();
	void showTotalDemandSatisfied();
	void showNumberOfRoutes();
	std::vector<int> getSolution();
	int getSatelliteNode();
	double getCost();
	int getMaxRouteCapacity();
	std::set<int> getCustomers();
	std::map<int, int> getCustomerTodemandMap();
	int getTotalDemandSatisfied();
	int getNumberOfRoutes();
};

//local solution class
class TwoEchelonSolution {
private:
	int numberOfActiveSatellites = 0;
	std::map<int, int> customerToDemandMap;
	std::map<int, int> satelliteToDemandMap;
	std::vector<std::vector<double>> distanceMatrix;
	std::set<int> customerNodes;
	std::set<int> satelliteNodes;
	CVRPSolution firstEchelonSolution;
	std::list<CVRPSolution> secondEchelonSolutions;
public:
	TwoEchelonSolution();
	TwoEchelonSolution(const TwoEchelonSolution & locSol);
	TwoEchelonSolution(int numberOfActiveSatellites, std::map<int, int> customerToDemandMap, std::map<int, int> satelliteToDemandMap, std::vector<std::vector<double>> distanceMatrix, std::set<int> customerNodes,	std::set<int> satelliteNodes, CVRPSolution firstEchelonSol, std::list<CVRPSolution> secondEchelonSols);
	int getNumberOfActiveSatellites();
	std::map<int, int> getCustomerToDemandMap();
	std::map<int, int> getSatelliteToDemandMap();
	std::vector<std::vector<double>> getDistanceMatrix();
	std::set<int> getCustomerNodes();
	std::set<int> getSatelliteNodes();
	CVRPSolution getFirstEchelonSolution();
	std::list<CVRPSolution> getSecondEchelonSolutions();
};

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

