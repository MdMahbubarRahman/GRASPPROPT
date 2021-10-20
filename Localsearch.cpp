#include "Localsearch.h"

//Implement Local Search algorithm here

//default constructor
CVRPSolution::CVRPSolution() {
	std::cout << "The default constructor of the CVRPSolution class has been called!" << std::endl;
}

//constructor
CVRPSolution::CVRPSolution(std::vector<int> sol, std::set<int> custs, std::map<int, int> custTodemand, double cost, int satNode, int maxCapacity, int demandSatisfied, int numRoutes){
	customers		    	= custs;
	customerTodemand        = custTodemand;
	solCost                 = cost;
	satelliteNode           = satNode;
	maxRouteCapacity        = maxCapacity;
	totalDemandSatisfied    = demandSatisfied;
	numberOfRoutes          = numRoutes;
	if (!sol.empty())
		solution = sol;
	else {
		std::cout << "The solution vector is empty! CVRP Solution object could not be formed." << std::endl;
	}
}

//copy constructor
CVRPSolution::CVRPSolution(const CVRPSolution& febSol) {
	customers				= febSol.customers;
	customerTodemand        = febSol.customerTodemand;
	solCost                 = febSol.solCost;
	satelliteNode           = febSol.satelliteNode;
	totalDemandSatisfied    = febSol.totalDemandSatisfied;
	numberOfRoutes          = febSol.numberOfRoutes;
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
	for (auto & it: customerTodemand) {
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
CustomerDepotDifferentialCost::CustomerDepotDifferentialCost() {
	customerID = 0;
	currentDepot = 0;
	potentialDepot = 0;
	differentialCost = 0;
}

//copy constructor
CustomerDepotDifferentialCost::CustomerDepotDifferentialCost(const CustomerDepotDifferentialCost& cusDifCost) {
	customerID = cusDifCost.customerID;
	currentDepot = cusDifCost.currentDepot;
	potentialDepot = cusDifCost.potentialDepot;
	differentialCost = cusDifCost.differentialCost;
}

//constructor
CustomerDepotDifferentialCost::CustomerDepotDifferentialCost(int cusID, int curDepot, int potDepot, double diffCost) {
	customerID = cusID;
	currentDepot = curDepot;
	potentialDepot = potDepot;
	differentialCost = diffCost;
}

//returns customer id
int CustomerDepotDifferentialCost::getCustomerID() {
	return customerID;
}

//returns current depot/satellite id
int CustomerDepotDifferentialCost::getCurrentDepot() {
	return currentDepot;
}

//returns potential depot/satellite id
int CustomerDepotDifferentialCost::getPotentialDepot() {
	return potentialDepot;
}

//returns differential cost
double CustomerDepotDifferentialCost::getDifferentialCost() {
	return differentialCost;
}

//shows customer to depots differential cost assignment
void CustomerDepotDifferentialCost::showCustomerDepotDifferentialCost() {
	std::cout << "The customer ID is : " << customerID << std::endl;
	std::cout << "The current depot/satellite ID is : " << currentDepot << std::endl;
	std::cout << "The potential depot/satellite ID is : " << potentialDepot << std::endl;
	std::cout << "The differential cost is : " << differentialCost << std::endl;
}

//comparator
bool Compare::operator()(CustomerDepotDifferentialCost & a, CustomerDepotDifferentialCost & b) {
	return (a.getDifferentialCost() > b.getDifferentialCost());
}


//default constructor
TwoEchelonSolution::TwoEchelonSolution() {
	std::cout << "The default constructor of the TwoEchelnSolution has been called!" << std::endl;
}

//copy constructor
TwoEchelonSolution::TwoEchelonSolution(const TwoEchelonSolution & locSol) {
	numberOfActiveSatellites = locSol.numberOfActiveSatellites;
	customerToDemandMap      = locSol.customerToDemandMap;
	distanceMatrix           = locSol.distanceMatrix;
	customerNodes			 = locSol.customerNodes;
	satelliteNodes			 = locSol.satelliteNodes;
	firstEchelonSolution	 = locSol.firstEchelonSolution;
	secondEchelonSolutions	 = locSol.secondEchelonSolutions;
}

//constructor
TwoEchelonSolution::TwoEchelonSolution(int numActiveSatellites, std::map<int, int> cusToDemandMap, std::map<int, int> satToDemandMap, std::vector<std::vector<double>> disMatrix, std::set<int> cusNodes, std::set<int> satNodes, CVRPSolution firstEchelonSol, std::list<CVRPSolution> secondEchelonSols){
	numberOfActiveSatellites = numActiveSatellites;
	customerToDemandMap = cusToDemandMap;
	satelliteToDemandMap = satToDemandMap;
	distanceMatrix = disMatrix;
	customerNodes = cusNodes;
	satelliteNodes = satNodes;
	firstEchelonSolution = firstEchelonSol;
	secondEchelonSolutions = secondEchelonSols;
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

//default constructor
Localsearch::Localsearch() {
	isNodeAdditionFeasible = false;
}

//copy constructor
Localsearch::Localsearch(const Localsearch& locsrch) {
	currentSolution = locsrch.currentSolution;
	bestSolution = locsrch.bestSolution;
	orderedCustomerList = locsrch.orderedCustomerList;
	firstEchelonSolution = locsrch.firstEchelonSolution;
	currentSatelliteSolution = locsrch.currentSatelliteSolution;
	potentialSatelliteSolution = locsrch.potentialSatelliteSolution;
	isNodeAdditionFeasible = locsrch.isNodeAdditionFeasible;
}

//constructor
Localsearch::Localsearch(TwoEchelonSolution currentSolution, TwoEchelonSolution bestSolution) {
	isNodeAdditionFeasible = false;
}

//returns current solution
TwoEchelonSolution Localsearch::getCurrentSolution(){
	return currentSolution;
}

//returns best solution
TwoEchelonSolution Localsearch::getBestSolution(){
	return bestSolution;
}

//show current solution
void Localsearch::showCurrentSolution() {
	//needs implementation
}

//show best solution
void Localsearch::showBestSolution() {
	//needs implementation
}

//run local search algorithm
void Localsearch::runLocalSearch() {
	//needs implementation
}

//orders customers based on reassignment costs
void Localsearch::customersReassignmentOrder() {
	//create customer to depot/satellite map
	std::map<int, int> customerTodepotDict;
	std::vector<int> firstEchelonSol = currentSolution.getFirstEchelonSolution().getSolution();
	int depotID = currentSolution.getFirstEchelonSolution().getSatelliteNode();
	for (auto it: firstEchelonSol) {
		if (currentSolution.getCustomerNodes().find(it) != currentSolution.getCustomerNodes().end()) {
			customerTodepotDict.insert(std::pair<int, int>(it, depotID));
		}
	}
	std::list<CVRPSolution> secondEchelonSolutions = currentSolution.getSecondEchelonSolutions();
	for (auto & it : secondEchelonSolutions) {
		int satelliteID = it.getSatelliteNode();
		for (auto iter : it.getSolution()) {
			if (iter != satelliteID) {
				customerTodepotDict.insert(std::pair<int, int>(iter, satelliteID));
			}
		}
	}
	//populate the customer order list/priority queue
	for (auto &it : customerTodepotDict) {
		int custID = it.first;
		int currentDepot = it.second;
		int potentialDepo = 0;
		double cost = 0;
		double currCost = currentSolution.getDistanceMatrix()[custID][currentDepot];
		double diffCost = currCost;
		for (auto sat : currentSolution.getSatelliteNodes()) {
			if (sat != currentDepot) {
				cost = currentSolution.getDistanceMatrix()[custID][sat]- currCost;
				if (cost < diffCost) {
					diffCost = cost;
					potentialDepo = sat;
				}
			}
		}
		CustomerDepotDifferentialCost cddc(custID, currentDepot, potentialDepo, diffCost);
		orderedCustomerList.push(cddc);
	}
}

//check whether the potential satellite can take an additional node
void Localsearch::checkNodeAdditionFeasibility() {
	int NodeDemand = 2;//need to extract this from demand vector
	std::vector<int> solution = potentialSatelliteSolution.getSolution();
	int capLimit = potentialSatelliteSolution.getMaxRouteCapacity();
	int satelliteID = potentialSatelliteSolution.getSatelliteNode();
	std::map<int, int> demand = currentSolution.getCustomerToDemandMap();
	//check feasibility
	int capacity = 0;
	for (auto it : solution) {
		if (it != satelliteID) {
			capacity += demand[it];
		}
		else {
			if (capacity != 0) {
				if (capacity+NodeDemand <= capLimit) {
					isNodeAdditionFeasible = true;
					break;
				}
			}
			capacity = 0;
		}
	}
}

//adds node to the potential satellite
void Localsearch::addNodeToThePotentialSatellite() {
	int node = 0;//node to be inserted
	int nodeDemand = 2;//needs update
	std::list<std::vector<int>> saturatedRoutes;
	std::list<std::vector<int>> unSaturatedRoutes;
	std::list<std::vector<int>> updatedUnSaturatedRoutes;
	std::vector<int> solution = potentialSatelliteSolution.getSolution();
	int satelliteID = potentialSatelliteSolution.getSatelliteNode();
	double cost = potentialSatelliteSolution.getCost();
	int capLimit = potentialSatelliteSolution.getMaxRouteCapacity();
	std::vector<std::vector<double>> distance = currentSolution.getDistanceMatrix();
	std::map<int, int> customerToDemand = potentialSatelliteSolution.getCustomerTodemandMap();
	//find the potential routes for insertion
	std::vector<int> route;
	int capacity = 0;
	for (auto it: solution) {
		if (route.size() == 0) {
			route.push_back(satelliteID);
		}
		if (it != satelliteID) {
			capacity += customerToDemand[it];
		}
		if (it == satelliteID && capacity != 0) {
			if (capacity+nodeDemand <= capLimit) {
				route.push_back(satelliteID);
				unSaturatedRoutes.push_back(route);
				route.clear();
				capacity = 0;
			}
			else {
				route.push_back(satelliteID);
				saturatedRoutes.push_back(route);
				route.clear();
				capacity = 0;
			}
		}
	}
	//insert node to the best position
	double costDiff = 100000;//big number
	int insertAt = 0;
	int routeNo = 0;
	int counter = 0;
	for (auto & iter: unSaturatedRoutes) {
		//calculate route cost
		double currentRouteCost = 0;
		for (int i = 0; i < (iter.size()-1); ++i) {
			currentRouteCost += distance[iter.at(i)][iter.at(i+1)];
		}
		//calculate route cost after insertion
		double newCost = 0;
		double newRouteCost = 100000;
		int insert = 0;
		for (int i = 0; i < (iter.size() - 1); ++i) {
			newCost = currentRouteCost - distance[iter.at(i)][iter.at(i + 1)] + distance[iter.at(i)][node] + distance[node][iter.at(i+1)];
			if (newCost < newRouteCost) {
				newRouteCost = newCost;
				insert = i;
			}
		}
		double diff = newRouteCost - currentRouteCost;
		if (diff < costDiff) {
			costDiff = diff;
			routeNo = counter;
			insertAt = insert;
		}
		counter++;
	}
	//update the best route
	counter = 0;
	std::vector<int> updatedRoute;
	for (auto & iter : unSaturatedRoutes) {
		if (counter == routeNo) {
			int i = 0;
			for (auto it: iter) {
				updatedRoute.push_back(it);
				if (i == insertAt) {
					updatedRoute.push_back(node);
				}
				i++;
			}
		}
		else {
			updatedUnSaturatedRoutes.push_back(iter);
		}
		counter++;
	}
	updatedUnSaturatedRoutes.push_back(updatedRoute);
	//updated the potential satellite solution
	std::vector<int> chrom_rep;
	for (auto &it: saturatedRoutes) {
		for (auto iter: it) {
			chrom_rep.push_back(iter);
		}
	}
	for (auto& it : updatedUnSaturatedRoutes) {
		for (auto iter : it) {
			chrom_rep.push_back(iter);
		}
	}
	cost = cost + costDiff;
	//Chromosome newChrom(chrom_rep, cost, sepInt, true, capLimit);
	//tabu search on newChrom can be applied
	//potentialSatelliteSolution = newChrom;
}

//deletes node from the current satellite
void Localsearch::deleteNodeFromTheCurrentSatellite() {
	std::vector<int> sol = currentSatelliteSolution.getSolution();
	int satelliteID = currentSatelliteSolution.getSatelliteNode();
	int capLimit = currentSatelliteSolution.getMaxRouteCapacity();
	double cost = currentSatelliteSolution.getCost();
	std::vector<std::vector<double>> distance = currentSolution.getDistanceMatrix();
	int node = 0;//node to be deleted
	std::vector<int> newSol;
	int i = 0;
	for (auto it: sol) {
		if (it != node) {
			newSol.push_back(it);
		}
		else {
			cost = cost - distance[sol.at(i - 1)][node] - distance[node][sol.at(i + 1)] + distance[sol.at(i - 1)][sol.at(i + 1)];
		}
		i++;
	}
	//CVRPSolution newChrom(newSol, cost, sepInt, true, capLimit);
	//tabu search on newChrom can be applied
	//currentSatelliteSolution = newChrom;
}

//updates first echelon solution 
void Localsearch::updateFirstEchelonSolution() {
	//NOTE: total demand field could be added in the chromosome or feasible solution class.
	//update demand for each satellite stations
	//check feasibility of the first echelon solution with the updated demand
	//if infeasible then update to feasible solution
	//run tabu search to optimize the current solution
	//
	//
	//
}


//compute objectve value of the new solution after local search
void Localsearch::computeObjectiveValueOfTheNewSolution() {
	//Local solution class should have a objective value field.
}