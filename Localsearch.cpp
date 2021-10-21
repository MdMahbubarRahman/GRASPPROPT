#include "Localsearch.h"

//Implement Local Search algorithm here


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