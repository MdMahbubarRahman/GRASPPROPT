#include "Initialsolution.h"

//Implement Initial solution algorihtm here 

//default constructor
ProblemParameters::ProblemParameters() {
	std::cout << "The default constructor of the problem parameter class has been called!";
}

//copy constructor
ProblemParameters::ProblemParameters(const ProblemParameters& probParam) {
	customerToDemandMap = probParam.customerToDemandMap;
	distanceMatrix = probParam.distanceMatrix;
	customerNodes = probParam.customerNodes;
	satelliteNodes = probParam.satelliteNodes;
	customersMustServeByFirstEchelon = probParam.customersMustServeByFirstEchelon;
	firstEchelonVehicleCapacityLimit = probParam.firstEchelonVehicleCapacityLimit;
	secondEchelonVehicleCapacityLimit = probParam.secondEchelonVehicleCapacityLimit;
	maxNumberOfVehicleInFirstEchelon = probParam.maxNumberOfVehicleInFirstEchelon;
	maxNumberOfVehicleInSecondEchelon = probParam.maxNumberOfVehicleInSecondEchelon;
}

//constructor
ProblemParameters::ProblemParameters(std::map<int, int> customerToDemandMap, std::vector<std::vector<double>> distanceMatrix, std::set<int> customerNodes, std::set<int> satelliteNodes, std::set<int> customersMustServeByFirstEchelon, int firstEchelonVehicleCapacityLimit, int secondEchelonVehicleCapacityLimit, int maxNumberOfVehicleInFirstEchelon, int maxNumberOfVehicleInSecondEchelon) {
	customerToDemandMap = customerToDemandMap;
	distanceMatrix = distanceMatrix;
	customerNodes = customerNodes;
	satelliteNodes = satelliteNodes;
	customersMustServeByFirstEchelon = customersMustServeByFirstEchelon;
	firstEchelonVehicleCapacityLimit = firstEchelonVehicleCapacityLimit;
	secondEchelonVehicleCapacityLimit = secondEchelonVehicleCapacityLimit;
	maxNumberOfVehicleInFirstEchelon = maxNumberOfVehicleInFirstEchelon;
	maxNumberOfVehicleInSecondEchelon = maxNumberOfVehicleInSecondEchelon;
}

//prints customer to demand map
void ProblemParameters::showCustomerToDemandMap() {
	std::cout << "The customer to demand map is given below : " << std::endl;
	for (auto & it: customerToDemandMap) {
		std::cout << "customer : " << it.first << "-->  demand : " << it.second << std::endl;
	}
}

//prints distance matrix
void ProblemParameters::showDistanceMatrix() {
	int i = 0;
	/*
	for (auto & it: distanceMatrix) {
		for (auto iter: it) {
			std::cout<<"distance from node : "<<i
		}
		i++;
	}
	*/
}

//prints customer nodes
void ProblemParameters::showCustomerNodes() {
	std::cout << "The customer nodes are : " << std::endl;
	for (auto it: customerNodes) {
		std::cout << it << std::endl;
	}
}

//prints satellite nodes
void ProblemParameters::showSatelliteNodes() {
	std::cout << "The satellite nodes are : " << std::endl;
	for (auto it : satelliteNodes) {
		std::cout << it << std::endl;
	}
}

//prints customers assigned to first echelon 
void ProblemParameters::showCustomersMustServeByFirstEchelon() {
	std::cout << "The customers assigned to first echelon vehicles are : " << std::endl;
	for (auto it : customersMustServeByFirstEchelon) {
		std::cout << it << std::endl;
	}
}

//prints first echelon vehicles' capacity limit
void ProblemParameters::showFirstEchelonVehicleCapacityLimit() {
	std::cout << "The capacity limit for a first echelon vehicle is : " << firstEchelonVehicleCapacityLimit << std::endl;
}

//prints capacity limit for second echelon vehicles
void ProblemParameters::showSecondEchelonVehicleCapacityLimit() {
	std::cout << "The capacity limit for a second echelon vehicle is : " << secondEchelonVehicleCapacityLimit << std::endl;
}

//prints max number of allowable first echelon vehicles
void ProblemParameters::showMaxNumberOfVehicleInFirstEchelon() {
	std::cout << "The maximum number of vehicles in the first echelon is : " << maxNumberOfVehicleInFirstEchelon << std::endl;
}

//prints max number of allowable second echelon vehicles per satellite
void ProblemParameters::showMaxNumberOfVehicleInSecondEchelon() {
	std::cout << "The maximum number of vehicles in the second echelon is : " << maxNumberOfVehicleInSecondEchelon << std::endl;
}

//returns customer to demand map
std::map<int, int> ProblemParameters::getCustomerToDemandMap() {
	return customerToDemandMap;
}

//returns distance matrix
std::vector<std::vector<double>> ProblemParameters::getDistanceMatrix() {
	return distanceMatrix;
}

//returns customer nodes
std::set<int> ProblemParameters::getCustomerNodes() {
	return customerNodes;
}

//returns satellite nodes
std::set<int> ProblemParameters::getSatelliteNodes() {
	return satelliteNodes;
}

//returns customers assigned for first echelon vehicles
std::set<int> ProblemParameters::getCustomersMustServeByFirstEchelon() {
	return customersMustServeByFirstEchelon;
}

//returns capacity limit for first echelon vehicles
int ProblemParameters::getFirstEchelonVehicleCapacityLimit() {
	return firstEchelonVehicleCapacityLimit;
}

//returns capacity limit for second echelon vehicles
int ProblemParameters::getSecondEchelonVehicleCapacityLimit() {
	return secondEchelonVehicleCapacityLimit;
}

//returns max number of allowable first echelon vehicles
int ProblemParameters::getMaxNumberOfVehicleInFirstEchelon() {
	return maxNumberOfVehicleInFirstEchelon;
}

//returns max number of allowable second echelon vehicles
int ProblemParameters::getMaxNumberOfVehicleInSecondEchelon() {
	return maxNumberOfVehicleInSecondEchelon;
}


















