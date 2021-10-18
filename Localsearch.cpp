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

//returns current depot/sattellite id
int CustomerDepotDifferentialCost::getCurrentDepot() {
	return currentDepot;
}

//returns potential depot/sattellite id
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
	std::cout << "The current depot/sattellite ID is : " << currentDepot << std::endl;
	std::cout << "The potential depot/sattellite ID is : " << potentialDepot << std::endl;
	std::cout << "The differential cost is : " << differentialCost << std::endl;
}


//comparator
bool Comparator::operator()(CustomerDepotDifferentialCost & a, CustomerDepotDifferentialCost & b) {
	return (a.getDifferentialCost() > b.getDifferentialCost());
}


//default constructor
LocalSolution::LocalSolution() {
	numberOfSattellites = 0;
	depotNode = 0;
}

//copy constructor
LocalSolution::LocalSolution(const LocalSolution& locSol) {
	numberOfSattellites		 = locSol.numberOfSattellites;
	depotNode				 = locSol.depotNode;
	demands					 = locSol.demands;
	distances				 = locSol.distances;//cost
	customerNodes			 = locSol.customerNodes;
	sattelliteNodes			 = locSol.sattelliteNodes;
	firstEchelonRoutes		 = locSol.firstEchelonRoutes;
	secondEchelonRoutes		 = locSol.secondEchelonRoutes;
	firstEchelonSolution	 = locSol.firstEchelonSolution;
	secondEchelonSolutions	 = locSol.secondEchelonSolutions;
}

//constructor
LocalSolution::LocalSolution(Chromosome firstEchelonSol, std::list<Chromosome> secondEchelonSols, std::vector<int> demands, std::vector<std::vector<double>> distances, std::vector<int> customerNodes, std::vector<int> sattelliteNodes) {
	numberOfSattellites = secondEchelonSols.size();
	depotNode = firstEchelonSol.getSepInt();
	demands = demands;
	distances = distances;
	customerNodes = customerNodes;
	sattelliteNodes = sattelliteNodes;
	std::vector<int> firstRouts = firstEchelonSol.getChromosomeRepresentation();
	for (auto it: firstRouts) {
		firstEchelonRoutes.push_back(it);
	}
	std::list<Chromosome> secondRouts = secondEchelonSols;
	for (auto & iter: secondRouts) {
		std::list<int> route;
		for (auto it: iter.getChromosomeRepresentation()) {
			route.push_back(it);
		}
		secondEchelonRoutes.push_back(route);
	}
	firstEchelonSolution = firstEchelonSol;
	secondEchelonSolutions = secondEchelonSols;
}

//returns the number of total sattellites
int LocalSolution::getNumberOfSattellites() {
	return numberOfSattellites;
}

//returns the depot node number
int LocalSolution::getDepotNode() {
	return depotNode;
}

//returns demand vector
std::vector<int> LocalSolution::getDemands() {
	return demands;
}

//returns distance matrix
std::vector<std::vector<double>> LocalSolution::getDistances() {
	return distances;
}

//returns customer nodes
std::set<int> LocalSolution::getCustomerNodes() {
	return customerNodes;
}

std::set<int> LocalSolution::getSattelliteNodes() {
	return sattelliteNodes;
}

//returns first echelon routes
std::list<int> LocalSolution::getFirstEchelonRoutes() {
	return firstEchelonRoutes;
}

//retuns second echelon routes
std::list<std::list<int>> LocalSolution::getSecondEchelonRoutes() {
	return secondEchelonRoutes;
}

//returns first echelon solution
Chromosome LocalSolution::getFirstEchelonSolution() {
	return firstEchelonSolution;
}

//returns second echelon solutions
std::list<Chromosome> LocalSolution::getSecondEchelonSolutions() {
	return secondEchelonSolutions;
}

//updates the local solution
void LocalSolution::updateLocalSolution(Chromosome firstEchelonSolution, std::list<Chromosome> secondEchelonSolutions, std::list<int> firstEchelonRoutes, std::list<std::list<int>> secondEchelonRoutes) {
	//needs implementation
}


//default constructor
Localsearch::Localsearch() {

}

//copy constructor
Localsearch::Localsearch(const Localsearch& locsrch) {

}

//constructor
Localsearch::Localsearch(LocalSolution currentSolution, LocalSolution bestSolution) {

}

//returns current solution
LocalSolution Localsearch::getCurrentSolution(){
	return currentSolution;
}

//returns best solution
LocalSolution Localsearch::getBestSolution(){
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
void Localsearch::orderCustomersBasedOnReassignmentCosts() {
	//create a customer depot map
	std::map<int, int> customerTodepotDict;
	std::vector<int> firstEchelonSol = currentSolution.getFirstEchelonSolution().getChromosomeRepresentation();
	int firstSepInt = currentSolution.getFirstEchelonSolution().getSepInt();
	for (auto it: firstEchelonSol) {
		if (currentSolution.getCustomerNodes().find(it) != currentSolution.getCustomerNodes().end()) {
			customerTodepotDict.insert(std::pair<int, int>(it, 0));
		}
	}
	std::list<Chromosome> secondEchelonSolutions = currentSolution.getSecondEchelonSolutions();
	for (auto & it : secondEchelonSolutions) {
		int sepInt = it.getSepInt();//sattellite id
		for (auto iter : it.getChromosomeRepresentation()) {
			if (iter != sepInt) {
				customerTodepotDict.insert(std::pair<int, int>(iter, sepInt));
			}
		}
	}
	//create a struct for customer id original depot potential depot differential cost
	



	//create a comparator for the struct/class considering differential cost

	//populate the priority queue

	//needs immplementation
	/*
	cutomer id
	differential cost 
	first sattellite/depot id
	second sattellite/depot id
	heap or priority queue can be used for this purpose
	customer id sattellite id map or dictionary
	need to define a comparator may be
	*/
}






