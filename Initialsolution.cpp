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


//default constructor
CVRPSolution::CVRPSolution() {
	std::cout << "The default constructor of the CVRPSolution class has been called!" << std::endl;
}

//constructor
CVRPSolution::CVRPSolution(std::vector<int> sol, std::set<int> custs, std::map<int, int> custTodemand, double cost, int satNode, int maxCapacity, int demandSatisfied, int numRoutes) {
	customers = custs;
	customerTodemand = custTodemand;
	solCost = cost;
	satelliteNode = satNode;
	maxRouteCapacity = maxCapacity;
	totalDemandSatisfied = demandSatisfied;
	numberOfRoutes = numRoutes;
	if (!sol.empty())
		solution = sol;
	else {
		std::cout << "The solution vector is empty! CVRP Solution object could not be formed." << std::endl;
	}
}

//copy constructor
CVRPSolution::CVRPSolution(const CVRPSolution& febSol) {
	customers = febSol.customers;
	customerTodemand = febSol.customerTodemand;
	solCost = febSol.solCost;
	satelliteNode = febSol.satelliteNode;
	totalDemandSatisfied = febSol.totalDemandSatisfied;
	numberOfRoutes = febSol.numberOfRoutes;
	if (!febSol.solution.empty())
		solution = febSol.solution;
	else {
		std::cout << "The solution vector is empty! CVRP Solution object could not be formed." << std::endl;
	}
}

//prints the feasible solution object
void CVRPSolution::showSolution() const {
	if (!solution.empty()) {
		std::cout << "The solution is : ";
		for (const auto& iter : solution) {
			std::cout << iter << " ";
		}
		std::cout << ";" << std::endl;
		std::cout << "The cost of the solution is : " << solCost << std::endl;
		std::cout << "The satellite node id is : " << satelliteNode << std::endl;
		std::cout << "The amount of the demand satisfied is : " << totalDemandSatisfied << std::endl;
		std::cout << "The number of vechile routes generated is : " << numberOfRoutes << std::endl;
		std::cout << "The number of customers satisfied is : " << customers.size() << std::endl;
	}
	else
		std::cout << "The solution is empty!" << std::endl;
}

//prints customers
void CVRPSolution::showCustomers() {
	std::cout << "The customers ID are : " << std::endl;
	for (auto it : customers) {
		std::cout << it << " ";
	}
	std::cout << ";" << std::endl;
}

//prints cutomer to demand map
void CVRPSolution::showCustomerTodemandMap() {
	std::cout << "The customer to demand map is : " << std::endl;
	for (auto& it : customerTodemand) {
		std::cout << "cutomer : " << it.first << "--->" << " demand : " << it.second << std::endl;
	}
}

//prints solution cost
void CVRPSolution::showSolCost() {
	std::cout << "The cost of the solution is : " << solCost << std::endl;
}

//prints separator integer value
void CVRPSolution::showSatelliteNode() {
	std::cout << "The satellite node id is : " << satelliteNode << std::endl;
}

//prints maximum route capacity in terms of demand satisfied
void CVRPSolution::showMaxRouteCapacity() {
	std::cout << "The maximum route capacity is : " << maxRouteCapacity << std::endl;
}

//prints the amount of demand has been satisfied by the solution
void CVRPSolution::showTotalDemandSatisfied() {
	std::cout << "The amount of demand satisfied is : " << totalDemandSatisfied << std::endl;
}

//prints the number of routes in the solution
void CVRPSolution::showNumberOfRoutes() {
	std::cout << "The number of routes in the solution is : " << numberOfRoutes << std::endl;
}

//returns solution vector
std::vector<int> CVRPSolution::getSolution() {
	return solution;
}

//returns separator integer value
int CVRPSolution::getSatelliteNode() {
	return satelliteNode;
}

//returns cost of the solution
double CVRPSolution::getCost() {
	return solCost;
}

//returns max route capacity
int CVRPSolution::getMaxRouteCapacity() {
	return maxRouteCapacity;
}

//returns customer set
std::set<int> CVRPSolution::getCustomers() {
	return customers;
}

//returns customer to demand map
std::map<int, int> CVRPSolution::getCustomerTodemandMap() {
	return customerTodemand;
}

//returns the amount of demand satisfied
int CVRPSolution::getTotalDemandSatisfied() {
	return totalDemandSatisfied;
}

//returns the number of routes present in the solution
int CVRPSolution::getNumberOfRoutes() {
	return numberOfRoutes;
}

//default constructor
TwoEchelonSolution::TwoEchelonSolution() {
	std::cout << "The default constructor of the TwoEchelnSolution has been called!" << std::endl;
}

//copy constructor
TwoEchelonSolution::TwoEchelonSolution(const TwoEchelonSolution& locSol) {
	solutionFitness = locSol.solutionFitness;
	numberOfActiveSatellites = locSol.numberOfActiveSatellites;
	customerToDemandMap = locSol.customerToDemandMap;
	distanceMatrix = locSol.distanceMatrix;
	customerNodes = locSol.customerNodes;
	satelliteNodes = locSol.satelliteNodes;
	firstEchelonSolution = locSol.firstEchelonSolution;
	secondEchelonSolutions = locSol.secondEchelonSolutions;
}

//constructor
TwoEchelonSolution::TwoEchelonSolution(double solutionFitness, int numActiveSatellites, std::map<int, int> cusToDemandMap, std::map<int, int> satToDemandMap, std::vector<std::vector<double>> disMatrix, std::set<int> cusNodes, std::set<int> satNodes, CVRPSolution firstEchelonSol, std::list<CVRPSolution> secondEchelonSols) {
	solutionFitness = solutionFitness;
	numberOfActiveSatellites = numActiveSatellites;
	customerToDemandMap = cusToDemandMap;
	satelliteToDemandMap = satToDemandMap;
	distanceMatrix = disMatrix;
	customerNodes = cusNodes;
	satelliteNodes = satNodes;
	firstEchelonSolution = firstEchelonSol;
	secondEchelonSolutions = secondEchelonSols;
}

//returns fitness value of the solution
int TwoEchelonSolution::getSolutionFitness() {
	return solutionFitness;
}

//returns the number of total satellites
int TwoEchelonSolution::getNumberOfActiveSatellites() {
	return numberOfActiveSatellites;
}

//returns customer to demand map
std::map<int, int> TwoEchelonSolution::getCustomerToDemandMap() {
	return customerToDemandMap;
}

//returns satellite to demand map
std::map<int, int> TwoEchelonSolution::getSatelliteToDemandMap() {
	return satelliteToDemandMap;
}

//returns distance  matrix
std::vector<std::vector<double>> TwoEchelonSolution::getDistanceMatrix() {
	return distanceMatrix;
}

//returns customer nodes
std::set<int> TwoEchelonSolution::getCustomerNodes() {
	return customerNodes;
}

std::set<int> TwoEchelonSolution::getSatelliteNodes() {
	return satelliteNodes;
}

//returns first echelon solution
CVRPSolution TwoEchelonSolution::getFirstEchelonSolution() {
	return firstEchelonSolution;
}

//returns second echelon solutions
std::list<CVRPSolution> TwoEchelonSolution::getSecondEchelonSolutions() {
	return secondEchelonSolutions;
}

//prints two echelon solution
void TwoEchelonSolution::showTwoEchelonSolution() {
	//needs immplementation
}

//default constructor
Initialsolution::Initialsolution() {
	std::cout << "The default constructor of the initial solution has been called!" << std::endl;
}

//copy constructor
Initialsolution::Initialsolution(const Initialsolution& initSol) {
	std::cout << "The copy constructor of the initial solution has been called!" << std::endl;
	febSol = initSol.febSol;
	chrom  = initSol.chrom;
	cvrpSol = initSol.cvrpSol;
	problemParameters = initSol.problemParameters;
	twoEchelonSolution = initSol.twoEchelonSolution;
}

//constructor
Initialsolution::Initialsolution(FeasibleSolution febSol, Chromosome chrom, CVRPSolution cvrpSol, ProblemParameters problemParameters, TwoEchelonSolution twoEchelonSolution) {
	febSol = febSol;
	chrom = chrom;
	cvrpSol = cvrpSol;
	problemParameters = problemParameters;
	twoEchelonSolution = twoEchelonSolution;
}

//returns feasible solution
FeasibleSolution Initialsolution::getFeasibleSolution() {
	return febSol;
}

//returns chromosome
Chromosome Initialsolution::getChromosome() {
	return chrom;
}

//returns cvrpsoluton
CVRPSolution Initialsolution::getCVRPSolution() {
	return cvrpSol;
}

//returns problem parameters
ProblemParameters Initialsolution::getProblemParameters() {
	return problemParameters;
}

//returns two echelon solution
TwoEchelonSolution Initialsolution::getTwoEchelonSolution() {
	return twoEchelonSolution;
}

