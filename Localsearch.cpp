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
	currentSolutionImproved = false;
	firstEchelonVehicleCapacityLimit = 0;
	secondEchelonVehicleCapacityLimit = 0;
	maxNumberOfVehicleInFirstEchelon = 0;
	maxNumberOfVehicleInSecondEchelon = 0;
}

//copy constructor
Localsearch::Localsearch(const Localsearch& locsrch) {
	currentSolution = locsrch.currentSolution;
	bestSolution = locsrch.bestSolution;
	orderedCustomerList = locsrch.orderedCustomerList;
	firstEchelonSolution = locsrch.firstEchelonSolution;
	currentSatelliteSolution = locsrch.currentSatelliteSolution;
	potentialSatelliteSolution = locsrch.potentialSatelliteSolution;
	currentSolutionImproved = locsrch.currentSolutionImproved;
	firstEchelonVehicleCapacityLimit = locsrch.firstEchelonVehicleCapacityLimit;
	secondEchelonVehicleCapacityLimit = locsrch.secondEchelonVehicleCapacityLimit;
	maxNumberOfVehicleInFirstEchelon = locsrch.maxNumberOfVehicleInFirstEchelon;
	maxNumberOfVehicleInSecondEchelon = locsrch.maxNumberOfVehicleInSecondEchelon;
}

//constructor
Localsearch::Localsearch(TwoEchelonSolution currSol, ProblemParameters probParams) {
	currentSolutionImproved = false;
	currentSolution = currSol;
	firstEchelonSolution = currSol.getFirstEchelonSolution();
	firstEchelonVehicleCapacityLimit = probParams.getFirstEchelonVehicleCapacityLimit();
	secondEchelonVehicleCapacityLimit = probParams.getSecondEchelonVehicleCapacityLimit();
	maxNumberOfVehicleInFirstEchelon = probParams.getMaxNumberOfVehicleInFirstEchelon();
	maxNumberOfVehicleInSecondEchelon = probParams.getMaxNumberOfVehicleInSecondEchelon();
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
	currentSolution.showTwoEchelonSolution();
}

//show best solution
void Localsearch::showBestSolution() {
	bestSolution.showTwoEchelonSolution();
}

//orders customers based on reassignment costs
void Localsearch::createCustomersPriorityQueue() {
	//create customer to depot/satellite map/Mohammad Ali
	std::map<int, int> customerToDepotDict;
	std::vector<int> firstEchelonSol = currentSolution.getFirstEchelonSolution().getSolution();
	int depotID = currentSolution.getFirstEchelonSolution().getSatelliteNode();
	std::set<int> dedicatedToDepotCustomers = currentSolution.getCustomersDedicatedToDepot();
	for (auto it : currentSolution.getFirstEchelonSolution().getCustomers()) {
		auto iter = dedicatedToDepotCustomers.find(it);
		if (iter == dedicatedToDepotCustomers.end()) {
			customerToDepotDict.insert(std::pair<int, int>(it, depotID));
		}
	}
	std::list<CVRPSolution> secondEchelonSolutions = currentSolution.getSecondEchelonSolutions();
	for (auto & it : secondEchelonSolutions) {
		int satelliteID = it.getSatelliteNode();
		for (auto iter : it.getCustomers()) {
			customerToDepotDict.insert(std::pair<int, int>(iter, satelliteID));
		}
	}
	//show the customer to depot map
	/*
	for (auto &it: customerToDepotDict) {
		std::cout << "customer : " << it.first << " depot/sat : " << it.second << std::endl;
	}
	*/
	//populate the customer order list/priority queue
	for (auto &it : customerToDepotDict) {
		int custID = it.first;
		std::set<int> satelliteNodes = currentSolution.getSatelliteNodes();
		auto iterator = satelliteNodes.find(custID);
		if (iterator == satelliteNodes.end()) {
			int currentDepot = it.second;
			int potentialDepo = 0;
			double cost = 0;
			double currCost = 0;
			currentDepot == 0 ? currCost = currentSolution.getRoadDistanceMatrix()[custID][currentDepot] : currCost = currentSolution.getAerialDistanceMatrix()[custID][currentDepot];
			double diffCost = INFINITY;
			for (auto sat : currentSolution.getSatelliteNodes()) {
				if (sat != currentDepot) {
					if (sat != custID) {
						if (sat == 0) {
							cost = currentSolution.getRoadDistanceMatrix()[custID][sat] - currCost;
						}             //RC_ij = c_ik-c_ij
						else {
							cost = currentSolution.getAerialDistanceMatrix()[custID][sat] - currCost;
						}
						if (cost < diffCost) {
							diffCost = cost;
							potentialDepo = sat;
						}
					}
				}
			}
			CustomerDepotDifferentialCost cddc(custID, currentDepot, potentialDepo, diffCost);
			orderedCustomerList.push(cddc);
		}
	}
}

//run local search algorithm
void Localsearch::runLocalSearch() {
	int counter = 0;
	while (!orderedCustomerList.empty()){
		counter++;
		CustomerDepotDifferentialCost cddc = orderedCustomerList.top();
		orderedCustomerList.pop();
		int cusID = cddc.getCustomerID();
		int currentDepot = cddc.getCurrentDepot();
		int potentialDepot = cddc.getPotentialDepot();
		double diffCost = cddc.getDifferentialCost();
		std::list<CVRPSolution> secondEchelonSolutions = currentSolution.getSecondEchelonSolutions();
		firstEchelonSolution = currentSolution.getFirstEchelonSolution();
		//populate with current sols
		if (currentDepot != currentSolution.getFirstEchelonSolution().getSatelliteNode()) {
			for (auto& it : secondEchelonSolutions) {
				if (it.getSatelliteNode() == currentDepot) {
					currentSatelliteSolution = it;
				}
			}
			if(potentialDepot != currentSolution.getFirstEchelonSolution().getSatelliteNode()){
				for(auto &it: secondEchelonSolutions){
					if (it.getSatelliteNode() == potentialDepot)
						potentialSatelliteSolution = it;
				}
			}
			else {
				potentialSatelliteSolution = firstEchelonSolution;
			}
		}
		else {
			currentSatelliteSolution = currentSolution.getFirstEchelonSolution();
			for (auto& it : secondEchelonSolutions) {
				if (it.getSatelliteNode() == potentialDepot) {
					potentialSatelliteSolution = it;
				}
			}
		}
		//update potential satellite  
		CVRPSolution updatedPotentialSatSol;   
		if (potentialDepot != firstEchelonSolution.getSatelliteNode()) {
			int demand = currentSatelliteSolution.getCustomerTodemandMap()[cusID];
			int satisfiedDemand = potentialSatelliteSolution.getTotalDemandSatisfied();
			//Check whether the potential satellite can adapt the new customer?
			if ((demand + satisfiedDemand) <= (secondEchelonVehicleCapacityLimit * maxNumberOfVehicleInSecondEchelon)) {
				std::map<int, int> cusToDemandMapUp = potentialSatelliteSolution.getCustomerTodemandMap();
				cusToDemandMapUp.insert(std::pair<int, int>(cusID, demand));
				std::vector<int> cusCluster;
				for (auto it : potentialSatelliteSolution.getCustomers()) {
					cusCluster.push_back(it);
				}
				cusCluster.push_back(cusID);
				Geneticalgorithm gap(potentialDepot, potentialSatelliteSolution.getMaxRouteCapacity(), cusToDemandMapUp, currentSolution.getAerialDistanceMatrix(), cusCluster);
				gap.runGeneticAlgorithm();
				Chromosome chromp = gap.getGASolution();
				std::set<int> customersp;
				for (auto it : cusCluster) {
					customersp.insert(it);
				}
				updatedPotentialSatSol = CVRPSolution(chromp, customersp, cusToDemandMapUp);
			}
			else {
				continue;
			}
		}
		//update current satellite
		CVRPSolution updatedCurrentSatSol;
		if (currentDepot != firstEchelonSolution.getSatelliteNode()) {
			std::map<int, int> cusToDemandMapCur = currentSatelliteSolution.getCustomerTodemandMap();
			std::map<int, int>::iterator iter;
			for (auto it = cusToDemandMapCur.begin(); it != cusToDemandMapCur.end(); ++it) {
				if ((*it).first == cusID) {
					iter = it;
					break;
				}
			}
			cusToDemandMapCur.erase(iter);
			std::vector<int> cusClusterCur;
			for (auto it : currentSatelliteSolution.getCustomers()) {
				if (it != cusID) {
					cusClusterCur.push_back(it);
				}
			}
			Geneticalgorithm gac(currentDepot, currentSatelliteSolution.getMaxRouteCapacity(), cusToDemandMapCur, currentSolution.getAerialDistanceMatrix(), cusClusterCur);
			gac.runGeneticAlgorithm();
			Chromosome chromc = gac.getGASolution();
			std::set<int> customersc;
			for (auto it : cusClusterCur) {
				customersc.insert(it);
			}
			updatedCurrentSatSol = CVRPSolution(chromc, customersc, cusToDemandMapCur);
		}
		//update first echelon solution
		int depot = firstEchelonSolution.getSatelliteNode();
		CVRPSolution updatedFirstEchelonSolution;
		if (depot == currentDepot) {
			int potentialSatDemand = updatedPotentialSatSol.getTotalDemandSatisfied();
			std::map<int, int> cusToDemandMap = firstEchelonSolution.getCustomerTodemandMap();
			std::map<int, int>::iterator iterp, iterd;
			for (auto it = cusToDemandMap.begin(); it != cusToDemandMap.end(); ++it) {
				if ((*it).first == potentialDepot) {
					iterp = it;
				}
				if ((*it).first == cusID) {
					iterd = it;
				}
			}
			cusToDemandMap.erase(iterp);
			cusToDemandMap.erase(iterd);
			if (potentialSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(potentialDepot, potentialSatDemand));
			}
			std::vector<int> customersVec;
			for (auto it : firstEchelonSolution.getCustomers()) {
				if (it != cusID) {
					customersVec.push_back(it);
				}
			}
			Geneticalgorithm gad(depot, firstEchelonSolution.getMaxRouteCapacity(), cusToDemandMap, currentSolution.getRoadDistanceMatrix(), customersVec);
			gad.runGeneticAlgorithm();
			Chromosome chromd = gad.getGASolution();
			std::set<int> customersd;
			for (auto it : customersVec) {
				customersd.insert(it);
			}
			updatedFirstEchelonSolution = CVRPSolution(chromd, customersd, cusToDemandMap);
			updatedCurrentSatSol = CVRPSolution(chromd, customersd, cusToDemandMap);
		}
		else if (depot == potentialDepot) {
			int currentSatDemand = updatedCurrentSatSol.getTotalDemandSatisfied();
			std::map<int, int> cusToDemandMap = firstEchelonSolution.getCustomerTodemandMap();
			std::map<int, int>::iterator iterc;
			for (auto it = cusToDemandMap.begin(); it != cusToDemandMap.end(); ++it) {
				if ((*it).first == currentDepot) {
					iterc = it;
				}
			}
			cusToDemandMap.erase(iterc);
			if (currentSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(currentDepot, currentSatDemand));
			}
			std::vector<int> customersVec;
			for (auto it : firstEchelonSolution.getCustomers()) {
				customersVec.push_back(it);
			}
			customersVec.push_back(cusID);
			int cusDemand = currentSolution.getCustomerToDemandMap()[cusID];
			cusToDemandMap.insert(std::pair<int, int>(cusID, cusDemand));
			Geneticalgorithm gad(depot, firstEchelonSolution.getMaxRouteCapacity(), cusToDemandMap, currentSolution.getRoadDistanceMatrix(), customersVec);
			gad.runGeneticAlgorithm();
			Chromosome chromd = gad.getGASolution();
			std::set<int> customersd;
			for (auto it : customersVec) {
				customersd.insert(it);
			}
			updatedFirstEchelonSolution = CVRPSolution(chromd, customersd, cusToDemandMap);
			updatedPotentialSatSol = CVRPSolution(chromd, customersd, cusToDemandMap);
		}
		else {
			int currentSatDemand = updatedCurrentSatSol.getTotalDemandSatisfied();
			int potentialSatDemand = updatedPotentialSatSol.getTotalDemandSatisfied();
			std::map<int, int> cusToDemandMap = firstEchelonSolution.getCustomerTodemandMap();
			std::map<int, int>::iterator iterc, iterp;
			for (auto it = cusToDemandMap.begin(); it != cusToDemandMap.end(); ++it) {
				if ((*it).first == currentDepot) {
					iterc = it;
				}
				if ((*it).first == potentialDepot) {
					iterp = it;
				}
			}
			cusToDemandMap.erase(iterc);
			cusToDemandMap.erase(iterp);
			std::vector<int> customersVec;
			for (auto it : firstEchelonSolution.getCustomers()) {
				if (it == currentDepot || it == potentialDepot) {
					continue;
				}
				else {
					customersVec.push_back(it);
				}
			}
			if (currentSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(currentDepot, currentSatDemand));
				customersVec.push_back(currentDepot);
			}
			if (potentialSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(potentialDepot, potentialSatDemand));
				customersVec.push_back(potentialDepot);
			}
			Geneticalgorithm gad(depot, firstEchelonSolution.getMaxRouteCapacity(), cusToDemandMap, currentSolution.getRoadDistanceMatrix(), customersVec);
			gad.runGeneticAlgorithm();
			Chromosome chromd = gad.getGASolution();
			std::set<int> customersd;
			for (auto it : customersVec) {
				customersd.insert(it);
			}
			updatedFirstEchelonSolution = CVRPSolution(chromd, customersd, cusToDemandMap);
		}
		//update twoechelon solution
		TwoEchelonSolution newTwoEchelonSol;
		newTwoEchelonSol.insertFirstEchelonSolution(updatedFirstEchelonSolution);
		std::map<int, int> satToDemandMap;
		satToDemandMap.insert(std::pair<int, int>(updatedFirstEchelonSolution.getSatelliteNode(), updatedFirstEchelonSolution.getTotalDemandSatisfied()));
		for (auto & it : secondEchelonSolutions) {
			if (it.getSatelliteNode() == currentDepot){
				newTwoEchelonSol.insertSecondEchelonSolution(updatedCurrentSatSol);
				satToDemandMap.insert(std::pair<int, int>(it.getSatelliteNode(), updatedCurrentSatSol.getTotalDemandSatisfied()));
			}
			else if (it.getSatelliteNode() == potentialDepot){
				newTwoEchelonSol.insertSecondEchelonSolution(updatedPotentialSatSol);
				satToDemandMap.insert(std::pair<int, int>(it.getSatelliteNode(), updatedPotentialSatSol.getTotalDemandSatisfied()));
			}
			else {
				newTwoEchelonSol.insertSecondEchelonSolution(it);
				satToDemandMap.insert(std::pair<int, int>(it.getSatelliteNode(), it.getTotalDemandSatisfied()));
			}
		}
		newTwoEchelonSol.populateTwoEchelonSolution(currentSolution.getCustomerToDemandMap(), satToDemandMap, currentSolution.getRoadDistanceMatrix(), currentSolution.getAerialDistanceMatrix(), currentSolution.getCustomerNodes(), currentSolution.getSatelliteNodes(), currentSolution.getCustomersDedicatedToDepot());
		//check with terminating conditions and/or update local solution
		if (newTwoEchelonSol.getSolutionFitness() < currentSolution.getSolutionFitness()) {
			currentSolution = newTwoEchelonSol;
			std::cout << "\nShow the customer shift parameters" << std::endl;
			cddc.showCustomerDepotDifferentialCost();
			currentSolutionImproved = true;
			std::cout << "\nNumber of neighbour investigated : " << counter << std::endl;
			break;
		}
		else {
			double newCost = newTwoEchelonSol.getSolutionFitness();
			double currentCost = currentSolution.getSolutionFitness();
			double costPercentChange = 0;
			costPercentChange = (newCost - currentCost) / currentCost;
			if (costPercentChange > 0.05) {
				std::cout << "\nNew Solution cost : " << newCost << std::endl;
				std::cout << "\nCurrent Solution cost : " << currentCost << std::endl;
				std::cout << "\nNumber of neighbour investigated : " << counter << std::endl;
				break;
			}
		}
	}
	std::cout << "\nNumber of neighbour investigated at all : " << counter << std::endl;
}

//returns whether the current solution improved
bool Localsearch::isCurrentSolutionImproved() {
	return currentSolutionImproved;
}


