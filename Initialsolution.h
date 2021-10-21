#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>

#include "Geneticalgorithm.h"

#ifndef INITIALSOLUTION_H
#define INITIALSOLUTION_H



/*
TODO: Implement the clustering constructive heuristic to get initial solution from 
the following paper:
	1. Multi-start heuristics for the two-echelon vehicle routing problem by
	TG Crainic, S Mancini, G Perboli, R Tadei

*/


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
	void showTwoEchelonSolution();
};

//initial solution
class Initialsolution {
private:
	FeasibleSolution febSol;
	Chromosome chrom;
	CVRPSolution cvrpSol;
	ProblemParameters problemParameters;
	TwoEchelonSolution twoEchelonSolution;
public:
	Initialsolution();
	Initialsolution(const Initialsolution & initSol);
	Initialsolution(FeasibleSolution febSol, Chromosome chrom, CVRPSolution cvrpSol, ProblemParameters problemParameters, TwoEchelonSolution twoEchelonSolution);
	FeasibleSolution getFeasibleSolution();
	Chromosome getChromosome();
	CVRPSolution getCVRPSolution();
	ProblemParameters getProblemParameters();
	TwoEchelonSolution getTwoEchelonSolution();
};

#endif // !INITIALSOLUTION_H





