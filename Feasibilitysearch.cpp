#include "Feasibilitysearch.h"

//Implement Feasibility search algorithm here.

//default constructor
Route::Route() {
	routeCost = 0;
	capacityServed = 0;
	sepInt = 0;
}

//constructor
Route::Route(std::list<int> rtSol, double rtCost, int capServed, int spInt) {
	routeSol = rtSol;
	routeCost = rtCost;
	capacityServed = capServed;
	sepInt = spInt;
}

//copy constructor
Route::Route(const Route& route) {
	routeSol = route.routeSol;
	routeCost = route.routeCost;
	capacityServed = route.capacityServed;
	sepInt = route.sepInt;
}

//returns route solution
std::list<int> Route::getRouteSol() {
	return routeSol;
}

//returns route cost
double Route::getRouteCost() {
	return routeCost;
}

//returns capacity
int Route::getCapacityServed() {
	return capacityServed;
}

//returns sep int
int Route::getSepInt() {
	return sepInt;
}

//prints the route properties
void Route::showRoute() {
	std::cout << "The satellite ID : " << sepInt << std::endl;
	std::cout << "The route is : " << std::endl;
	for (auto it: routeSol) {
		std::cout << it << " ";
	}
	std::cout << ";" << std::endl;
	std::cout << "The cost of the route is : " << routeCost << std::endl;
	std::cout << "The amount of capacity served is : " << capacityServed << std::endl;
}

void Route::updateRoute(int cus, int cap, double additionalCost) {
	routeSol.push_back(cus);
	capacityServed += cap;
	routeCost += additionalCost;
}

//default constructor
Feasibilitysearch::Feasibilitysearch() {
	isFeasible = false;
}

//constructor
Feasibilitysearch::Feasibilitysearch(ProblemParameters probParams, TwoEchelonSolution currSol) {
	problemParams = probParams;
	currentSolution = currSol;
	isFeasible = false;
}

//copy constructor
Feasibilitysearch::Feasibilitysearch(const Feasibilitysearch& febSearch) {
	problemParams = febSearch.problemParams;
	currentSolution = febSearch.currentSolution;
	feasibilitySearchSolution = febSearch.feasibilitySearchSolution;
	isFeasible = febSearch.isFeasible;
}

//returns feasibility status
bool Feasibilitysearch::getFeasibilityStatus() {
	return isFeasible;
}

//returns current solution
TwoEchelonSolution Feasibilitysearch::getCurrentSolution() {
	return currentSolution;
}

//returns feasibility search solution
TwoEchelonSolution Feasibilitysearch::getFeasibilitySearchSolution() {
	return feasibilitySearchSolution;
}

//prints feasibility status
void Feasibilitysearch::showFeasibilityStatus() {
	std::cout << "Feasibility status : " << isFeasible << std::endl;
}

//prints current solution
void Feasibilitysearch::showCurrentSolution() {
	currentSolution.showTwoEchelonSolution();
}

//runs feasibility search algorithm
void Feasibilitysearch::runFeasibilitySearch() {
	bool firstEchelonSolFeasible = false;
	TwoEchelonSolution febSrchSol;
	std::list<CVRPSolution> secondEchelonSols = currentSolution.getSecondEchelonSolutions();
	std::list<CVRPSolution> infeasibleSecondEchelonSols;
	//check initial feasibility
	currentSolution.getFirstEchelonSolution().getNumberOfRoutes() <= problemParams.getMaxNumberOfVehicleInFirstEchelon() ? firstEchelonSolFeasible = true : firstEchelonSolFeasible = false;
	for (auto & it : secondEchelonSols) {
		if (it.getNumberOfRoutes() > problemParams.getMaxNumberOfVehicleInSecondEchelon()) {
			infeasibleSecondEchelonSols.push_back(it);
		}
	}
	(!infeasibleSecondEchelonSols.empty() || firstEchelonSolFeasible == false) ? isFeasible = false : isFeasible = true;
	//feasibility search
	if (isFeasible == false) {
		std::set<int> freeCustomersSet;
		CVRPSolution firstEchelonCVRPSol;
		CVRPSolution secondEchelonCVRPSol;
		std::list<Route> firstEchelonRoutesList;
		std::list<std::list<Route>> secondEchelonRoutesLists;
		//populate first echelon routes list
		if (firstEchelonSolFeasible == false || firstEchelonSolFeasible == true) {
			int sepInt = currentSolution.getFirstEchelonSolution().getSatelliteNode();
			std::vector<int> firstEchelonSol = currentSolution.getFirstEchelonSolution().getSolution();
			std::list<int> routeSol;
			double cost = 0;
			int capacity = 0;
			int prev = sepInt;
			for (auto it : firstEchelonSol) {
				if (it != sepInt && routeSol.size() == 0) {
					cost = currentSolution.getDistanceMatrix()[sepInt][it];
					capacity = currentSolution.getCustomerToDemandMap()[it];
					routeSol.push_back(it);
					prev = it;
				}
				if (it != sepInt && routeSol.size() > 0) {
					cost += currentSolution.getDistanceMatrix()[prev][it];
					capacity += currentSolution.getCustomerToDemandMap()[it];
					routeSol.push_back(it);
					prev = it;
				}
				if (it == sepInt && routeSol.size() > 0) {
					cost += currentSolution.getDistanceMatrix()[prev][it];
					Route route(routeSol, cost, capacity, sepInt);
					firstEchelonRoutesList.push_back(route);
					capacity = 0;
					cost = 0;
					routeSol.clear();
					prev = sepInt;
				}
			}
			if (routeSol.size() > 0) {
				cost += currentSolution.getDistanceMatrix()[prev][sepInt];
				Route route(routeSol, cost, capacity, sepInt);
				firstEchelonRoutesList.push_back(route);
			}
		}
		//populate second echelon routes list
		for (auto & iter : infeasibleSecondEchelonSols) {
			std::list<Route> satelliteRouteList;
			int sepInt = iter.getSatelliteNode();
			std::vector<int> solution = iter.getSolution();
			std::list<int> routeSol;
			double cost = 0;
			int capacity = 0;
			int prev = sepInt;
			for (auto it : solution) {
				if (it != sepInt && routeSol.size() == 0) {
					cost = currentSolution.getDistanceMatrix()[sepInt][it];
					capacity = currentSolution.getCustomerToDemandMap()[it];
					routeSol.push_back(it);
					prev = it;
				}
				if (it != sepInt && routeSol.size() > 0) {
					cost += currentSolution.getDistanceMatrix()[prev][it];
					capacity += currentSolution.getCustomerToDemandMap()[it];
					routeSol.push_back(it);
					prev = it;
				}
				if (it == sepInt && routeSol.size() > 0) {
					cost += currentSolution.getDistanceMatrix()[prev][it];
					Route route(routeSol, cost, capacity, sepInt);
					satelliteRouteList.push_back(route);
					capacity = 0;
					cost = 0;
					routeSol.clear();
					prev = sepInt;
				}
			}
			if (routeSol.size() > 0) {
				cost += currentSolution.getDistanceMatrix()[prev][sepInt];
				Route route(routeSol, cost, capacity, sepInt);
				satelliteRouteList.push_back(route);
			}	
			secondEchelonRoutesLists.push_back(satelliteRouteList);
		}
		//make infeasible satellite solutions feasible
		for (auto & iter : secondEchelonRoutesLists) {
			std::set<int> customers;
			while (iter.size() > problemParams.getMaxNumberOfVehicleInSecondEchelon()) {
				std::list<Route>::iterator itr;
				int cap = 10000;
				for (auto it = iter.begin(); it != iter.end(); ++it) {
					if ((*it).getCapacityServed() < cap) {
						cap = (*it).getCapacityServed();
						itr = it;
					}
				}
				for (auto cus : (*itr).getRouteSol()) {
					customers.insert(cus);
				}
				iter.erase(itr);
			}
			//insert free customers in the routes if possible
			std::set<int> deleteableItem;
			for (auto cus : customers) {
				int cap = currentSolution.getCustomerToDemandMap()[cus];
				for (auto & it : iter) {
					if (cap <= (problemParams.getSecondEchelonVehicleCapacityLimit() - it.getCapacityServed())) {
						int lastElmt = it.getRouteSol().back();
						int sepInt = it.getSepInt();
						double addCost = currentSolution.getDistanceMatrix()[lastElmt][cus] + currentSolution.getDistanceMatrix()[cus][sepInt] - currentSolution.getDistanceMatrix()[lastElmt][sepInt];
						it.updateRoute(cus, cap, addCost);
						deleteableItem.insert(cus);
						break;
					}
				}
			}
			//delete the inserted customers
			for (auto it: deleteableItem) {
				customers.erase(it);
			}
			//update free customers list 
			for (auto it: customers) {
				freeCustomersSet.insert(it);
			}
		}
		//make first echelon solution feasible
		while (firstEchelonRoutesList.size() > problemParams.getMaxNumberOfVehicleInFirstEchelon()) {
			std::list<Route>::iterator itr;
			int cap = 10000;
			for (auto it = firstEchelonRoutesList.begin(); it != firstEchelonRoutesList.end(); ++it) {
				if ((*it).getCapacityServed() < cap) {
					cap = (*it).getCapacityServed();
					itr = it;
				}
			}
			for (auto cus : (*itr).getRouteSol()) {
				freeCustomersSet.insert(cus);
			}
			firstEchelonRoutesList.erase(itr);
		}
		//insert free customers in the routes if possible
		std::set<int> deleteableItem;
		for (auto cus : freeCustomersSet) {
			int cap = currentSolution.getCustomerToDemandMap()[cus];
			for (auto & it : firstEchelonRoutesList) {
				if (cap <= (problemParams.getFirstEchelonVehicleCapacityLimit() - it.getCapacityServed())) {
					int lastElmt = it.getRouteSol().back();
					int sepInt = it.getSepInt();
					double addCost = currentSolution.getDistanceMatrix()[lastElmt][cus] + currentSolution.getDistanceMatrix()[cus][sepInt] - currentSolution.getDistanceMatrix()[lastElmt][sepInt];
					it.updateRoute(cus, cap, addCost);
					deleteableItem.insert(cus);
					break;
				}
			}
		}
		//delete the inserted customers
		for (auto it : deleteableItem) {
			freeCustomersSet.erase(it);
		}
		//update final feasibility status
		freeCustomersSet.size() == 0 ? isFeasible = true : isFeasible = false;
		//update second echelon cvrp sols and first echelon cvrp sol
	}
	else {
		feasibilitySearchSolution = currentSolution;
	}
}




