#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <queue>

#include "Geneticalgorithm.h"

#ifndef INITIALSOLUTION_H
#define INITIALSOLUTION_H



/*
TODO: Implement the clustering constructive heuristic to get initial solution from 
the following paper:
	1. Multi-start heuristics for the two-echelon vehicle routing problem by
	TG Crainic, S Mancini, G Perboli, R Tadei

*/


//customer to satellite distance class
class CustomerToSatelliteDistance {
private:
	int customerID;
	int satelliteID;
	double distance;
public:
	CustomerToSatelliteDistance();
	CustomerToSatelliteDistance(const CustomerToSatelliteDistance& cusToSatDis);
	CustomerToSatelliteDistance(int cusID, int satelliteID, double disCost);
	int getCustomerID();
	int getSatelliteID();
	double getDistance();
	void showCustomerToSatelliteDistance();
};

//comparator for customer to satellite distance
class Distant {
public:
	bool operator()(CustomerToSatelliteDistance & a, CustomerToSatelliteDistance & b);
};

//Input parameters
class ProblemParameters {
private:
	std::map<int, int> customerToDemandMap;
	std::vector<std::vector<double>> distanceMatrix;
	std::set<int> customerNodes;
	std::set<int> satelliteNodes;
	std::set<int> customersMustServeByFirstEchelon;
	int firstEchelonVehicleCapacityLimit = 0;
	int secondEchelonVehicleCapacityLimit = 0;
	int maxNumberOfVehicleInFirstEchelon = 0;
	int maxNumberOfVehicleInSecondEchelon = 0;
public:
	ProblemParameters();
	ProblemParameters(const ProblemParameters & probParam);
	ProblemParameters(std::map<int, int> customerToDemandMap, std::vector<std::vector<double>> distanceMatrix, std::set<int> customerNodes,	std::set<int> satelliteNodes, std::set<int> customersMustServeByFirstEchelon, int firstEchelonVehicleCapacityLimit, int secondEchelonVehicleCapacityLimit, int maxNumberOfVehicleInFirstEchelon, int maxNumberOfVehicleInSecondEchelon);
	void showCustomerToDemandMap();
	void showDistanceMatrix();
	void showCustomerNodes();
	void showSatelliteNodes();
	void showCustomersMustServeByFirstEchelon();
	void showFirstEchelonVehicleCapacityLimit();
	void showSecondEchelonVehicleCapacityLimit();
	void showMaxNumberOfVehicleInFirstEchelon();
	void showMaxNumberOfVehicleInSecondEchelon();
	std::map<int, int> getCustomerToDemandMap();
	std::vector<std::vector<double>> getDistanceMatrix();
	std::set<int> getCustomerNodes();
	std::set<int> getSatelliteNodes();
	std::set<int> getCustomersMustServeByFirstEchelon();
	int getFirstEchelonVehicleCapacityLimit();
	int getSecondEchelonVehicleCapacityLimit();
	int getMaxNumberOfVehicleInFirstEchelon();
	int getMaxNumberOfVehicleInSecondEchelon();
};


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
	CVRPSolution(const CVRPSolution& sol); //copy constructor
	CVRPSolution(Chromosome & chromSol, std::set<int> customers, std::map<int, int> customerTodemand); //constructor
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
	double solutionFitness = 0;
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
	TwoEchelonSolution(const TwoEchelonSolution& locSol);
	TwoEchelonSolution(double solutionFitness, int numberOfActiveSatellites, std::map<int, int> customerToDemandMap, std::map<int, int> satelliteToDemandMap, std::vector<std::vector<double>> distanceMatrix, std::set<int> customerNodes, std::set<int> satelliteNodes, CVRPSolution firstEchelonSol, std::list<CVRPSolution> secondEchelonSols);
	int getSolutionFitness();
	int getNumberOfActiveSatellites();
	std::map<int, int> getCustomerToDemandMap();
	std::map<int, int> getSatelliteToDemandMap();
	std::vector<std::vector<double>> getDistanceMatrix();
	std::set<int> getCustomerNodes();
	std::set<int> getSatelliteNodes();
	CVRPSolution getFirstEchelonSolution();
	std::list<CVRPSolution> getSecondEchelonSolutions();
	void insertFirstEchelonSolution(CVRPSolution cvrp);
	void insertSecondEchelonSolution(CVRPSolution cvrp);
	void showTwoEchelonSolution();
	void clearSecondEchelonSolutionList();
	void populateTwoEchelonSolution(std::map<int, int> customerToDemandMap, std::map<int, int> satelliteToDemandMap, std::vector<std::vector<double>> distanceMatrix, std::set<int> customerNodes, std::set<int> satelliteNodes);
};

//initial solution
class Initialsolution {
private:
	bool isProblemFeasible;
	std::multimap<int, int> satelliteToCustomerMap;
	std::map<int, int> satelliteToDemandMap;
	std::set<int> freeCustomers;
	Chromosome chrom;
	CVRPSolution cvrpSol;
	ProblemParameters problemParameters;
	TwoEchelonSolution twoEchelonSolution;
	std::priority_queue<CustomerToSatelliteDistance, std::vector<CustomerToSatelliteDistance>, Distant> customerToSatDistList;
public:
	Initialsolution();
	Initialsolution(const Initialsolution & initSol);
	Initialsolution(Chromosome chrom, CVRPSolution cvrpSol, ProblemParameters problemParameters, TwoEchelonSolution twoEchelonSolution);
	Initialsolution(ProblemParameters problemParameters);
	Chromosome getChromosome();
	CVRPSolution getCVRPSolution();
	ProblemParameters getProblemParameters();
	TwoEchelonSolution getTwoEchelonSolution();
	bool getFeasibilityStatus();
	std::multimap<int, int> getSatelliteToCustomerMap();
	std::map<int, int> getSatelliteToDemandMap();
	std::set<int> getFreeCustomers();
	std::priority_queue<CustomerToSatelliteDistance, std::vector<CustomerToSatelliteDistance>, Distant> getCustomerToSatDistList();
	void mapCustomersToSatellites();
	void makeFeasibleCustomersCluster();
	void generateCVRPSolutions();
	void generateTwoEchelonSolution();
	void runInitialSolution();
};

#endif // !INITIALSOLUTION_H


