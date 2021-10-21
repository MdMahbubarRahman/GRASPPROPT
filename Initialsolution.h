#include <iostream>
#include <map>
#include <vector>
#include <set>

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





#endif // !INITIALSOLUTION_H





