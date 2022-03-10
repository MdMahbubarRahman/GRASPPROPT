#include "Initialsolution.h"

//Implement Initial solution algorihtm here 

//default constructor
CustomerToSatelliteDistance::CustomerToSatelliteDistance() {
	customerID = 0;
	satelliteID = 0;
	distance = 0;
}

//copy constructor
CustomerToSatelliteDistance::CustomerToSatelliteDistance(const CustomerToSatelliteDistance & cusToSatDis) {
	customerID = cusToSatDis.customerID;
	satelliteID = cusToSatDis.satelliteID;
	distance = cusToSatDis.distance;
}

//constructor
CustomerToSatelliteDistance::CustomerToSatelliteDistance(int cusID, int satID, double disCost) {
	customerID = cusID;
	satelliteID = satID;
	distance = disCost;
}

//returns customer id
int CustomerToSatelliteDistance::getCustomerID() {
	return customerID;
}

//returns satellite id
int CustomerToSatelliteDistance::getSatelliteID() {
	return satelliteID;
}

//returns customer to sat distance
double CustomerToSatelliteDistance::getDistance() {
	return distance;
}

//prints customer to satellite distance
void CustomerToSatelliteDistance::showCustomerToSatelliteDistance() {
	std::cout << "The customer ID is : " << customerID << std::endl;
	std::cout << "The satellite ID is : " << satelliteID << std::endl;
	std::cout << "The distance is : " << distance << std::endl;
}

//distant operator()/comparator
bool Distant::operator()(CustomerToSatelliteDistance & a, CustomerToSatelliteDistance & b) {
	return (a.getDistance() > b.getDistance());
}

//default constructor
ProblemParameters::ProblemParameters() {
	std::cout << "The default constructor of the problem parameter class has been called!" << std::endl;
}

//copy constructor
ProblemParameters::ProblemParameters(const ProblemParameters& probParam) {
	customerToDemandMap = probParam.customerToDemandMap;
	roadDistanceMatrix = probParam.roadDistanceMatrix;
	aerialDistanceMatrix = probParam.aerialDistanceMatrix;
	customerNodes = probParam.customerNodes;
	satelliteNodes = probParam.satelliteNodes;
	customersMustServeByFirstEchelon = probParam.customersMustServeByFirstEchelon;
	firstEchelonVehicleCapacityLimit = probParam.firstEchelonVehicleCapacityLimit;
	secondEchelonVehicleCapacityLimit = probParam.secondEchelonVehicleCapacityLimit;
	maxNumberOfVehicleInFirstEchelon = probParam.maxNumberOfVehicleInFirstEchelon;
	maxNumberOfVehicleInSecondEchelon = probParam.maxNumberOfVehicleInSecondEchelon;
}

//constructor
ProblemParameters::ProblemParameters(std::map<int, int> cusToDemandMap, std::vector<std::vector<double>> roadDistMatrix, std::vector<std::vector<double>> aerialDistMatrix, std::set<int> cusNodes, std::set<int> satNodes, std::set<int> cusMustServeByFirstEchelon, int firstEchelonVehicleCapLimit, int secondEchelonVehicleCapLimit, int maxNumOfVehicleInFirstEchelon, int maxNumOfVehicleInSecondEchelon) {
	customerToDemandMap = cusToDemandMap;
	roadDistanceMatrix = roadDistMatrix;
	aerialDistanceMatrix = aerialDistMatrix;
	customerNodes = cusNodes;
	satelliteNodes = satNodes;
	customersMustServeByFirstEchelon = cusMustServeByFirstEchelon;
	firstEchelonVehicleCapacityLimit = firstEchelonVehicleCapLimit;
	secondEchelonVehicleCapacityLimit = secondEchelonVehicleCapLimit;
	maxNumberOfVehicleInFirstEchelon = maxNumOfVehicleInFirstEchelon;
	maxNumberOfVehicleInSecondEchelon = maxNumOfVehicleInSecondEchelon;
}

//prints customer to demand map
void ProblemParameters::showCustomerToDemandMap() {
	std::cout << "The customer to demand map is given below : " << std::endl;
	for (auto & it: customerToDemandMap) {
		std::cout << "customer : " << it.first << "-->  demand : " << it.second << std::endl;
	}
}

//prints distance matrix
void ProblemParameters::showRoadDistanceMatrix() {
	int i = 0;
	std::cout << "\nneeds implementation." << std::endl;
	/*
	for (auto & it: distanceMatrix) {
		for (auto iter: it) {
			std::cout<<"distance from node : "<<i
		}
		i++;
	}
	*/
}

//prints distance matrix
void ProblemParameters::showAerialDistanceMatrix() {
	int i = 0;
	std::cout << "\nneeds implementation." << std::endl;
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
std::vector<std::vector<double>> ProblemParameters::getRoadDistanceMatrix() {
	return roadDistanceMatrix;
}

//returns distance matrix
std::vector<std::vector<double>> ProblemParameters::getAerialDistanceMatrix() {
	return aerialDistanceMatrix;
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
CVRPSolution::CVRPSolution(const CVRPSolution & febSol) {
	customers = febSol.customers;
	customerTodemand = febSol.customerTodemand;
	solCost = febSol.solCost;
	satelliteNode = febSol.satelliteNode;
	totalDemandSatisfied = febSol.totalDemandSatisfied;
	numberOfRoutes = febSol.numberOfRoutes;
	maxRouteCapacity = febSol.maxRouteCapacity;
	if (!febSol.solution.empty())
		solution = febSol.solution;
	else {
		std::cout << "The solution vector is empty! CVRP Solution object could not be formed." << std::endl;
	}
}

//constructor
CVRPSolution::CVRPSolution(Chromosome & chromSol, std::set<int> cuss, std::map<int, int> cusTodemand) {
	solution = chromSol.getChromosomeRepresentation();
	maxRouteCapacity = chromSol.getMaxRouteCapacity();
	customers = cuss;
	customerTodemand = cusTodemand;
	solCost = chromSol.getFitness();
	satelliteNode = chromSol.getSepInt();
	totalDemandSatisfied = 0;
	for (auto it : cuss) {
		totalDemandSatisfied += cusTodemand[it];
	}
	int routeID = 1;
	bool flag = false;
	for (auto it = solution.begin(); it != solution.end(); ++it) {
		if (*it != satelliteNode) {
			flag = true;
		}
		else {
			if (flag) {
				routeID += 1;
				flag = false;
			}
		}
	}
	flag == true ? numberOfRoutes = routeID : numberOfRoutes = (routeID - 1);
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
		std::cout << "The maximum capacity of a route is : " << maxRouteCapacity << std::endl;
		//no need to show customer to demand
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
	satelliteToDemandMap = locSol.satelliteToDemandMap;
	roadDistanceMatrix = locSol.roadDistanceMatrix;
	aerialDistanceMatrix = locSol.aerialDistanceMatrix;
	customerNodes = locSol.customerNodes;
	satelliteNodes = locSol.satelliteNodes;
	customersDedicatedToDepot = locSol.customersDedicatedToDepot;
	firstEchelonSolution = locSol.firstEchelonSolution;
	secondEchelonSolutions = locSol.secondEchelonSolutions;
}

//constructor
TwoEchelonSolution::TwoEchelonSolution(double solutionFitness, int numActiveSatellites, std::map<int, int> cusToDemandMap, std::map<int, int> satToDemandMap, std::vector<std::vector<double>> roadDisMatrix, std::vector<std::vector<double>> aerialDisMatrix, std::set<int> cusNodes, std::set<int> satNodes, std::set<int> dedicatedCus, CVRPSolution firstEchelonSol, std::list<CVRPSolution> secondEchelonSols) {
	solutionFitness = solutionFitness;
	numberOfActiveSatellites = numActiveSatellites;
	customerToDemandMap = cusToDemandMap;
	satelliteToDemandMap = satToDemandMap;
	roadDistanceMatrix = roadDisMatrix;
	aerialDistanceMatrix = aerialDisMatrix;
	customerNodes = cusNodes;
	satelliteNodes = satNodes;
	customersDedicatedToDepot = dedicatedCus;
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
std::vector<std::vector<double>> TwoEchelonSolution::getRoadDistanceMatrix() {
	return roadDistanceMatrix;
}

//returns distance  matrix
std::vector<std::vector<double>> TwoEchelonSolution::getAerialDistanceMatrix() {
	return aerialDistanceMatrix;
}

//returns customer nodes
std::set<int> TwoEchelonSolution::getCustomerNodes() {
	return customerNodes;
}

std::set<int> TwoEchelonSolution::getSatelliteNodes() {
	return satelliteNodes;
}

std::set<int> TwoEchelonSolution::getCustomersDedicatedToDepot() {
	return customersDedicatedToDepot;
}

//returns first echelon solution
CVRPSolution TwoEchelonSolution::getFirstEchelonSolution() {
	return firstEchelonSolution;
}

//returns second echelon solutions
std::list<CVRPSolution> TwoEchelonSolution::getSecondEchelonSolutions() {
	return secondEchelonSolutions;
}

//inserts first echelon solution
void TwoEchelonSolution::insertFirstEchelonSolution(CVRPSolution cvrp) {
	firstEchelonSolution = cvrp;
}

//inserts second echelon solution
void TwoEchelonSolution::insertSecondEchelonSolution(CVRPSolution cvrp) {
	secondEchelonSolutions.push_back(cvrp);
}

//prints two echelon solution
void TwoEchelonSolution::showTwoEchelonSolution() {
	std::cout << "\nThe cost of the two echelon solution is : " << solutionFitness << std::endl;
	std::cout << "\nThe first echelon solution is : " << std::endl;
	firstEchelonSolution.showSolution();
	std::cout << "\nThe second echelon solutions are : " << std::endl;
	int i = 1;
	for (auto &it: secondEchelonSolutions) {
		std::cout << "\nSecond echelon cvrp sol no : " << i << std::endl;
		it.showSolution();
		i++;
	}
}

//clears second echelon solutions
void TwoEchelonSolution::clearSecondEchelonSolutionList() {
	secondEchelonSolutions.clear();
}

//populates two echelon solution
void TwoEchelonSolution::populateTwoEchelonSolution(std::map<int, int> cusToDemandMap, std::map<int, int> satToDemandMap, std::vector<std::vector<double>> roadDisMatrix, std::vector<std::vector<double>> aerialDisMatrix, std::set<int> cusNodes, std::set<int> satNodes, std::set<int> dedicatedCusToDepo) {
	numberOfActiveSatellites = secondEchelonSolutions.size();
	customerToDemandMap = cusToDemandMap;
	satelliteToDemandMap = satToDemandMap;
	roadDistanceMatrix = roadDisMatrix;
	aerialDistanceMatrix = aerialDisMatrix;
	customerNodes = cusNodes;
	satelliteNodes = satNodes;
	customersDedicatedToDepot = dedicatedCusToDepo;
	solutionFitness = firstEchelonSolution.getCost();
	for (auto & it : secondEchelonSolutions) {
		solutionFitness += it.getCost();
	}
}

//default constructor
Initialsolution::Initialsolution() {
	std::cout << "The default constructor of the initial solution has been called!" << std::endl;
	isProblemFeasible = true;
}

//copy constructor
Initialsolution::Initialsolution(const Initialsolution& initSol) {
	std::cout << "The copy constructor of the initial solution has been called!" << std::endl;
	isProblemFeasible = initSol.isProblemFeasible;
	chrom  = initSol.chrom;
	cvrpSol = initSol.cvrpSol;
	problemParameters = initSol.problemParameters;
	twoEchelonSolution = initSol.twoEchelonSolution;
}

//constructor
Initialsolution::Initialsolution(Chromosome chrm, CVRPSolution cvrpSl, ProblemParameters problemParam, TwoEchelonSolution twoEchelonSol) {
	isProblemFeasible = true;
	chrom = chrm;
	cvrpSol = cvrpSl;
	problemParameters = problemParam;
	twoEchelonSolution = twoEchelonSol;
}

//constructor
Initialsolution::Initialsolution(ProblemParameters probParameters) {
	problemParameters = probParameters;
	isProblemFeasible = true;
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

//returns satellite to customer map
std::multimap<int, int> Initialsolution::getSatelliteToCustomerMap() {
	return satelliteToCustomerMap;
}

//rerturns satellite to demand map
std::map<int, int> Initialsolution::getSatelliteToDemandMap() {
	return satelliteToDemandMap;
}

//returns free customers
std::set<int> Initialsolution::getFreeCustomers() {
	return freeCustomers;
}

//returns cutomer to satellite distance priority queue
std::priority_queue<CustomerToSatelliteDistance, std::vector<CustomerToSatelliteDistance>, Distant> Initialsolution::getCustomerToSatDistList() {
	return customerToSatDistList;
}

//returns feasibility status
bool Initialsolution::getFeasibilityStatus() {
	return isProblemFeasible;
}

//map customers to its nearest satelliltes
void Initialsolution::mapCustomersToSatellites() {
	std::set<int> customers = problemParameters.getCustomerNodes();
	std::set<int> satellites = problemParameters.getSatelliteNodes();
	std::vector<std::vector<double>> roadDistance = problemParameters.getRoadDistanceMatrix();
	std::vector<std::vector<double>> aerialDistance = problemParameters.getAerialDistanceMatrix();
	std::set<int> firstEchelonCustomers = problemParameters.getCustomersMustServeByFirstEchelon();
	//depot node is set to 0 by default
	//assign customers to its nearest satellites
	for (auto it : customers) {
		//std::cout << it << std::endl;
		if (firstEchelonCustomers.find(it) != firstEchelonCustomers.end()) {
			satelliteToCustomerMap.insert(std::make_pair(0, it));
			std::cout << "First echelon priority found!" << std::endl;
		}
		else {
			double cost = INFINITY;
			double variableCost = 0;
			int sat = 0;
			for (auto iter : satellites) {
				if (iter == 0) {//depot sat == 0 one of the assumptions
					variableCost = roadDistance[it][iter];
				}
				else {
					variableCost = aerialDistance[it][iter];
				}
				//assign satellite to the current customer
				if (variableCost < cost) {
					cost = variableCost;
					sat = iter;
				}
			}
			satelliteToCustomerMap.insert(std::make_pair(sat, it));
		}	
	}
	std::cout << "\nThe customer to satellite map is done" << std::endl;
	std::cout << "The map is : " << std::endl;
	for (auto & it: satelliteToCustomerMap) {
		std::cout << "satellite : " << it.first << " customer : " << it.second << std::endl;
	}
}

//update customers clusters so that they satisfy vehicle capacity and max vehicle constraints
void Initialsolution::makeFeasibleCustomersCluster() {
	int maxNoFirstEchelonVehicles = problemParameters.getMaxNumberOfVehicleInFirstEchelon();
	int maxNoSecondEchelonVehicles = problemParameters.getMaxNumberOfVehicleInSecondEchelon();
	int firstVehicleCapLimit = problemParameters.getFirstEchelonVehicleCapacityLimit();
	int secondVehicleCapLimit = problemParameters.getSecondEchelonVehicleCapacityLimit();
	std::map<int, int> customerToDemandMap = problemParameters.getCustomerToDemandMap();
	std::set<int> satellites = problemParameters.getSatelliteNodes();
	std::vector<std::vector<double>> roadDistance = problemParameters.getRoadDistanceMatrix();
	std::vector<std::vector<double>> aerialDistance = problemParameters.getAerialDistanceMatrix();
	//feasibility check of the clusters and update if necessary
	for (auto it : satellites) {
		if (it != 0) {
			auto itlow = satelliteToCustomerMap.lower_bound(it);
			auto ithigh = satelliteToCustomerMap.upper_bound(it);
			int capacity = 0;
			int size = 0;
			for (auto & iter = itlow; iter != ithigh; ++iter) {
				//std::cout << "satellite : " << (*iter).first << " customer : " << (*iter).second << std::endl;
				capacity += customerToDemandMap[(*iter).second];
				//std::cout << "capacity : " << capacity << std::endl;
				CustomerToSatelliteDistance csd((*iter).second, (*iter).first, aerialDistance[(*iter).second][(*iter).first]);
				//csd.showCustomerToSatelliteDistance();
				customerToSatDistList.push(csd);
				size++;
			}
			if (capacity <= (secondVehicleCapLimit*maxNoSecondEchelonVehicles)) {
				std::cout << "The customer assignment is feasible." << std::endl;
				satelliteToDemandMap.insert(std::pair<int, int>(it, capacity));
				while (!customerToSatDistList.empty()) {
					customerToSatDistList.pop();
				}
			}
			else {
				capacity = 0;
				for (int i = 0; i < size; ++i) {
					CustomerToSatelliteDistance csd = customerToSatDistList.top();
					std::cout << "In the else looop" << std::endl;
					if ((capacity + customerToDemandMap[csd.getCustomerID()]) <= (secondVehicleCapLimit * maxNoSecondEchelonVehicles)) {
						capacity += customerToDemandMap[csd.getCustomerID()];
						customerToSatDistList.pop();
					}
					else {
						satelliteToDemandMap.insert(std::pair<int, int>(it, capacity));
						break;
					}	
				}
			}
			while (!customerToSatDistList.empty()) {
				CustomerToSatelliteDistance csd = customerToSatDistList.top();
				auto it = satelliteToCustomerMap.lower_bound(csd.getSatelliteID());
				for (auto & iter = it; iter != satelliteToCustomerMap.upper_bound(csd.getSatelliteID()); ++iter) {
					if ((*iter).second == csd.getCustomerID()) {
						satelliteToCustomerMap.erase(iter);
						break;
					}
				}
				customerToSatDistList.pop();	
				freeCustomers.insert(csd.getCustomerID());
			}
		}
		else {
			//calculate customers demand load to the depot
			auto itlow = satelliteToCustomerMap.lower_bound(it);
			auto ithigh = satelliteToCustomerMap.upper_bound(it);
			int capacity = 0;
			for (auto & iter = itlow; iter != ithigh; ++iter) {
				capacity += customerToDemandMap[(*iter).second];
			}
			satelliteToDemandMap.insert(std::pair<int, int>(it, capacity));
		}
	}
	//show free customers
	std::cout << "\nFree customers are : " << std::endl;
	for (auto it : freeCustomers) {
		std::cout << it << std::endl;
	}
	//show total demand
	int demand = 0;
	for (auto & it : satelliteToDemandMap) {
		demand += it.second;
	}
	std::cout << "\nInitial allocated total demand is : " << demand << std::endl;
	//assign free customers to satellites
	std::map<int, int> satelliteToAvailableCapacityMap;
	for (auto it: satellites) {
		if (it != 0) {
			satelliteToAvailableCapacityMap.insert(std::pair<int, int>(it, (maxNoSecondEchelonVehicles * secondVehicleCapLimit) - satelliteToDemandMap[it]));
			std::cout << "satellitle : " << it << " satisfied demand : " << satelliteToDemandMap[it] << std::endl;
		}
		else {
			satelliteToAvailableCapacityMap.insert(std::pair<int, int>(it, (maxNoFirstEchelonVehicles * firstVehicleCapLimit) - satelliteToDemandMap[it]));
			std::cout << "depot : " << it << " satisfied demand : " << satelliteToDemandMap[it] << std::endl;
		}
	}
	//show available capacity in each satellite
	std::cout << "\nInitial satellite to available capacity map : " << std::endl;
	for (auto it : satellites) {
		std::cout << "satellilte : " << it << " available capacity : " << satelliteToAvailableCapacityMap[it] << std::endl;
	}
	while (!freeCustomers.empty()){
		auto cus = freeCustomers.begin();
		double cost = INFINITY;
		double difCost = 0;
		int sat = INFINITY;
		for (auto it: satellites) {
			if (it == 0) {
				difCost = roadDistance[*cus][it];
			}
			else {
				difCost = aerialDistance[*cus][it];
			}			
			if (difCost < cost && (customerToDemandMap[*cus] <= satelliteToAvailableCapacityMap[it])) {
				cost = difCost;
				sat = it;
			}
		}
		if (satellites.find(sat) != satellites.end()) {
			int currentLoad = satelliteToDemandMap[sat] + customerToDemandMap[*cus];
			satelliteToDemandMap.erase(sat);
			satelliteToDemandMap.insert(std::pair<int, int>(sat, currentLoad));
			int currentCap = satelliteToAvailableCapacityMap[sat] - customerToDemandMap[*cus];
			satelliteToAvailableCapacityMap.erase(sat);
			satelliteToAvailableCapacityMap.insert(std::pair<int, int>(sat, currentCap));
			satelliteToCustomerMap.insert(std::make_pair(sat, *cus));
			freeCustomers.erase(cus);
		}
		else {
			std::cout << "\nThe Two Echelon CVRP problem is infeasible! The satellites available capacity is very limited!" << std::endl;
			isProblemFeasible = false;
			break;
		}
	}
	//show satellite to customer cluster
	std::cout << "\nShow feasible satellite to customer map : " << std::endl;
	for (auto & it: satelliteToCustomerMap) {
		std::cout << "satellite : " << it.first << " customer : " << it.second<<" demand : " << customerToDemandMap[it.second] << std::endl;
	}
	std::cout << "\nIntermediate satellite to available capacity map : " << std::endl;
	for (auto it : satellites) {
		std::cout << "satellilte : " << it << " available capacity : " << satelliteToAvailableCapacityMap[it] << std::endl;
	}
	std::wcout << "\nIntermediate satellite to demand map : " << std::endl;
	for (auto it: satellites) {
		std::cout << "satellite : " << it << " satisfied demand : " << satelliteToDemandMap[it] << std::endl;
	}
	//update depot assignment with satellites
	int totalDemand = 0;
	for (auto it: satellites) {
		if (it != 0) {
			totalDemand += satelliteToDemandMap[it];
			satelliteToCustomerMap.insert(std::make_pair(0, it));
			customerToDemandMap.insert(std::pair<int, int>(it, satelliteToDemandMap[it]));
		}
	}
	totalDemand += satelliteToDemandMap[0];
	satelliteToDemandMap.erase(0);
	satelliteToDemandMap.insert(std::pair<int, int>(0, totalDemand));
	//show final satellite to demand distribution
	std::wcout << "\nFinal satellite to demand map : " << std::endl;
	for (auto it : satellites) {
		std::cout << "satellite : " << it << " satisfied demand : " << satelliteToDemandMap[it] << std::endl;
	}
	//show final satellite to customer distribution
	std::cout << "\nFinal satellite to customer map : " << std::endl;
	for (auto& it : satelliteToCustomerMap) {
		std::cout << "satellite : " << it.first << " customer : " << it.second << " demand : " << customerToDemandMap[it.second] << std::endl;
	}
}

//generates cvrp solution
void Initialsolution::generateCVRPSolutions() {
	std::vector<std::vector<double>> roadDistance = problemParameters.getRoadDistanceMatrix();
	std::vector<std::vector<double>> aerialDistance = problemParameters.getAerialDistanceMatrix();
	int capLimitFirst = problemParameters.getFirstEchelonVehicleCapacityLimit();
	int capLimitSecond = problemParameters.getSecondEchelonVehicleCapacityLimit();
	std::map<int, int> cusToDemand = problemParameters.getCustomerToDemandMap();
	//update cusToDemand map
	for (auto & it : satelliteToDemandMap) {
		if (it.first != 0) {
			cusToDemand.insert(it);
		}
	}
	//check feasibility
	if (isProblemFeasible == false) {
		std::cout << "\nThe problem is infeasible!" << std::endl;
		std::cout << "The algorithm terminated." << std::endl;
	}
	else {
		//solve cvrp problem
		twoEchelonSolution.clearSecondEchelonSolutionList();
		for (auto iter : problemParameters.getSatelliteNodes()) {
			int satID = iter;
			std::vector<int> customerSet;
			std::vector<std::vector<double>> distance;
			for (auto it = satelliteToCustomerMap.lower_bound(iter); it != satelliteToCustomerMap.upper_bound(iter); ++it) {
				customerSet.push_back((*it).second);
			}
			int capLimit = 0;
			satID == 0 ? capLimit = capLimitFirst : capLimit = capLimitSecond;
			satID == 0 ? distance = roadDistance : distance = aerialDistance;
			std::cout << "\nsatellite : " << satID << " capacity limit for a route : " << capLimit << std::endl;
			std::cout << "\nThe customers assigned to this satellite are : " << std::endl;
			for (auto it: customerSet) {
				std::cout << "customer : " << it <<" demand : "<<cusToDemand[it]<< std::endl;
			}
			std::cout << "\nGenetic Algorithm started" << std::endl;
			Geneticalgorithm ga(satID, capLimit, cusToDemand, distance, customerSet);
			ga.runGeneticAlgorithm();
			std::cout << "\nGenetic Algoirthm terminated" << std::endl;
			chrom = ga.getGASolution();
			std::set<int> customers;
			for (auto it: customerSet) {
				customers.insert(it);
			}
			cvrpSol = CVRPSolution(chrom, customers, cusToDemand);
			satID == 0 ? twoEchelonSolution.insertFirstEchelonSolution(cvrpSol) : twoEchelonSolution.insertSecondEchelonSolution(cvrpSol);
		}
	}
}

//generates two echelon solution
void Initialsolution::generateTwoEchelonSolution() {
	twoEchelonSolution.populateTwoEchelonSolution(problemParameters.getCustomerToDemandMap(), satelliteToDemandMap, problemParameters.getRoadDistanceMatrix(), problemParameters.getAerialDistanceMatrix(), problemParameters.getCustomerNodes(), problemParameters.getSatelliteNodes(), problemParameters.getCustomersMustServeByFirstEchelon());
}

//run initial solution algoirthm to get initial two echelon solution
void Initialsolution::runInitialSolution() {
	mapCustomersToSatellites();
	makeFeasibleCustomersCluster();
	generateCVRPSolutions();
	generateTwoEchelonSolution();
	twoEchelonSolution.showTwoEchelonSolution();
}


