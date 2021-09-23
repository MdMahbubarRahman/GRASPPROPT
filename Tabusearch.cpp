#include "Tabusearch.h"


//default constructor
TabuList::TabuList() {}

//constructor
TabuList::TabuList(int routeNum, int customerID) {
	tabuList.insert(std::make_pair(routeNum, customerID));
}

//shows content of the tabulist
void TabuList::showTabuList() const {
	if (tabuList.size() == 0) {
		std::cout << "The tabu list is empty!" << std::endl;
	}
	else {
		std::cout << "Contents of the tabu list are :" << std::endl;
		for (auto& elm : tabuList) {
			std::cout << elm.first << "  " << elm.second << std::endl;
		}
	}
}

//adds tabu move to the tabu list
void TabuList::addToTabuList(int routeNum, int customerID) {
	tabuList.insert(std::make_pair(routeNum, customerID));
}

//deletes contents of the tabu list
void TabuList::clearTabuList() {
	tabuList.clear();
}

//checks whether a move is in the list
bool TabuList::checkInTabuList(int routeNum, int customerID) {
	bool flag = false;
	auto itlow = tabuList.lower_bound(routeNum);  
	auto itup = tabuList.upper_bound(routeNum);   
	if (itlow != tabuList.end())
		for (auto &it = itlow; it != itup; it++) {
			if ((*it).second == customerID)
				return flag = true;
		}
	else
		return flag;	
}

//default constructor
AspirationCriteria::AspirationCriteria() {
	currentBestCost = 100000000.0;//Big number
}

//updates aspiration criteria
void AspirationCriteria::updateCurrentBestCost(double cost) {
	if (cost < currentBestCost)
		currentBestCost = cost;
}

//shows current aspiration criteria
double AspirationCriteria::showCurrentBestCost() {
	return currentBestCost;
}

//default constructor
FeasibleSolution::FeasibleSolution() {
	cost = 1000000.0;
}

//constructor
FeasibleSolution::FeasibleSolution(std::vector<int> sol, double cost) {
	solution = sol;
	cost = cost;
}

//FeasibleSolution::FeasibleSolution(const FeasibleSolution & sol) {
//	TODO: implementation needs to be done.
//}

//prints the solution
void FeasibleSolution::showSolution() const {
	if (solution.size() != 0) {
		std::cout << "The solution is : " << std::endl;
		for (auto& iter : solution) {
			std::cout << iter << std::endl;
		}
	}
	else
		std::cout << "The solution is empty!" << std::endl;
}

//prints cost of the solution
void FeasibleSolution::showCost() {
	std::cout << cost << std::endl; 
}

//returns cost of the solution
double FeasibleSolution::getCost() {
	return cost;
}

//clears the feasible solution vector
void FeasibleSolution::clearSolution() {

}

//default constructor
Neighbourhood::Neighbourhood() {

}

// copy constructor
Neighbourhood::Neighbourhood(const Neighbourhood& neighbour) {


}

//inserts a new neighbour/solution
void Neighbourhood::insertNeighbour(FeasibleSolution neighbour) {

}

//shows neighbouring solutions
void Neighbourhood::showNeighbours() {

}

//inserts neighbour solutions to the kneighbour solution list
void Neighbourhood::insertToKNeighbour(std::list<FeasibleSolution> neighbourSolution) {


}

//shows all k neighbour solutions
void Neighbourhood::showKNeighbours() {


}

//returns the best solution from the neighbourhood
FeasibleSolution Neighbourhood::getBestFromNeighbour() {

}

//returns the best solution from the K neighbourhood
FeasibleSolution Neighbourhood::getBestFromKNeighbour() {

}

//clears the list of neighbour solutions
void Neighbourhood::clearNeighbourSolution() {


}

//clears the list of K neighbour solutions
void Neighbourhood::clearKNeighbourSolution() {


}

//default constructor
Tabusearch::Tabusearch() {
	k = 2;
}

//constructor
Tabusearch::Tabusearch(FeasibleSolution feasblesol) {
	k = 2;
}

//copy constructor
Tabusearch::Tabusearch(const Tabusearch& tabusrch) {
	k = 2;
}

//updates incumbent solution
void Tabusearch::updateIncumbentSolution() {

}

//generates Neighbour solutions
//k-chain moves uses Add and Drop as well as 1-Swap methods
void Tabusearch::generateKNeighbourSolutions() {

}

//runs tabu search algoirthm
void Tabusearch::tabuSearchRun(FeasibleSolution febSol) {

}

//performs TSP heuristics on incumbent solution
void Tabusearch::performTspHeuristics() {

}

//clears tabu search class data
void Tabusearch::clearTabuSearchData() {

}
