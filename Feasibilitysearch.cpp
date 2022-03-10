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
		else {
			//std::cout << "CVRP solution is being inserted in the new two echelon solution!" << std::endl;
			febSrchSol.insertSecondEchelonSolution(it);
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
		std::cout << "List of routes for first echelon and second echelon cvrp solutions have been created." << std::endl;
		//populate first echelon routes list
		//if (firstEchelonSolFeasible == false || firstEchelonSolFeasible == true) {
		std::cout << "Show the cvrp solution for the first echelon : " << std::endl;
		currentSolution.getFirstEchelonSolution().showSolution();
		int sepInt = currentSolution.getFirstEchelonSolution().getSatelliteNode();
		std::vector<int> firstEchelonSol = currentSolution.getFirstEchelonSolution().getSolution();
		std::list<int> routeSol;
		double cost = 0;
		int capacity = 0;
		int prev = sepInt;
		std::set<int> satNodes = currentSolution.getSatelliteNodes();
		for (auto it : firstEchelonSol) {
			if (it != sepInt) {
				cost += currentSolution.getRoadDistanceMatrix()[prev][it];
				auto itrr = satNodes.find(it);
				if (itrr != satNodes.end()) {
					capacity += currentSolution.getSatelliteToDemandMap()[it];
				}
				else {
					capacity += currentSolution.getCustomerToDemandMap()[it];
				}
				routeSol.push_back(it);
				prev = it;
			}
			if (it == sepInt) {
				if (routeSol.size() != 0) {
					cost += currentSolution.getRoadDistanceMatrix()[prev][it];
					Route route(routeSol, cost, capacity, sepInt);
					firstEchelonRoutesList.push_back(route);
					capacity = 0;
					cost = 0;
					routeSol.clear();
					prev = sepInt;
				}	
			}
		}
		if (routeSol.size() > 0) {
			cost += currentSolution.getRoadDistanceMatrix()[prev][sepInt];
			Route route(routeSol, cost, capacity, sepInt);
			firstEchelonRoutesList.push_back(route);
		}
		//show the routes for the first echelon
		/*
		for (auto &it : firstEchelonRoutesList) {
			std::cout << "show the route : " << std::endl;
			it.showRoute();
		}
		*/
		//}
		//show infeasible second echelon cvrp sols
		/*
		for (auto & iter: infeasibleSecondEchelonSols) {
			std::cout << "show second echelon cvrp infeasible sol : " << std::endl;
			iter.showSolution();
		}
		*/
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
				if (it != sepInt) {
					cost += currentSolution.getAerialDistanceMatrix()[prev][it];
					capacity += currentSolution.getCustomerToDemandMap()[it];
					routeSol.push_back(it);
					prev = it;
				}
				if (it == sepInt) {
					if (routeSol.size() != 0) {
						cost += currentSolution.getAerialDistanceMatrix()[prev][it];
						Route route(routeSol, cost, capacity, sepInt);
						satelliteRouteList.push_back(route);
						capacity = 0;
						cost = 0;
						routeSol.clear();
						prev = sepInt;
					}
				}
			}
			if (routeSol.size() > 0) {
				cost += currentSolution.getAerialDistanceMatrix()[prev][sepInt];
				Route route(routeSol, cost, capacity, sepInt);
				satelliteRouteList.push_back(route);
			}	
			/*
			std::cout << "show second echelon route : " << std::endl;
			for (auto & it: satelliteRouteList) {
				std::cout << "show route : " << std::endl;
				it.showRoute();
			}
			*/
			secondEchelonRoutesLists.push_back(satelliteRouteList);
		}
		//make infeasible satellite solutions feasible
		for (auto & iter : secondEchelonRoutesLists) {
			std::set<int> customers;
			std::cout << "\nInside the second echelon routes lists" << std::endl;
			while (iter.size() > problemParams.getMaxNumberOfVehicleInSecondEchelon()) {
				std::list<Route>::iterator itr;
				int cap = 100000;
				//std::cout << "\nInside while loop" << std::endl;
				//remove minimum capacity serving route from the solution
				for (auto it = iter.begin(); it != iter.end(); ++it) {
					int routeCap = (*it).getCapacityServed();
					if (routeCap < cap) { 
						cap = routeCap;
						itr = it;
					}
				}
				for (auto cus : (*itr).getRouteSol()) {
					customers.insert(cus);
				}
				iter.erase(itr);
				//std::cout << "\nAfter Erase operation" << std::endl;
			}
			std::cout << "\nInsert free customers in the routes if possible" << std::endl;
			std::set<int> deleteableItem;
			for (auto cus : customers) {
				int cap = currentSolution.getCustomerToDemandMap()[cus];
				for (auto & it : iter) {
					if (cap <= (problemParams.getSecondEchelonVehicleCapacityLimit() - it.getCapacityServed())) {
						int lastElmt = it.getRouteSol().back();
						int sepInt = it.getSepInt();
						double addCost = currentSolution.getAerialDistanceMatrix()[lastElmt][cus] + currentSolution.getAerialDistanceMatrix()[cus][sepInt] - currentSolution.getAerialDistanceMatrix()[lastElmt][sepInt];
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
			
			for (auto & it : iter) {
				std::cout << "show the updated routes : " << std::endl;
				it.showRoute();
			}
			
		}
		//make first echelon solution feasible
		while (firstEchelonRoutesList.size() > problemParams.getMaxNumberOfVehicleInFirstEchelon()) {
			std::list<Route>::iterator itr;
			int cap = 100000;
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
					double addCost = currentSolution.getRoadDistanceMatrix()[lastElmt][cus] + currentSolution.getRoadDistanceMatrix()[cus][sepInt] - currentSolution.getRoadDistanceMatrix()[lastElmt][sepInt];
					it.updateRoute(cus, cap, addCost);
					deleteableItem.insert(cus);
					break;
				}
			}
		}
		/*
		std::cout << "show the updated first echelon routesList : " << std::endl;
		for (auto & it: firstEchelonRoutesList) {
			std::cout << "show updated route : " << std::endl;
			it.showRoute();
		}
		*/
		//delete the inserted customers
		for (auto it : deleteableItem) {
			freeCustomersSet.erase(it);
		}
		//update final feasibility status
		freeCustomersSet.size() == 0 ? isFeasible = true : isFeasible = false;
		std::cout << "\nThe feasibility status is : " << isFeasible << std::endl;
		//update second echelon cvrp sols and first echelon cvrp sol
		if (isFeasible == true) {
			//update second echelon cvrps
			for (auto& it : secondEchelonRoutesLists) {
				std::vector<int> solution;
				std::map<int, int> cusToDemand;
				std::set<int> customer;
				int capacity = 0;
				double cost = 0;
				int satID = 0;
				for (auto& iter : it) {
					satID = iter.getSepInt();
					for (auto itr : iter.getRouteSol()) {
						solution.push_back(itr);
						cusToDemand.insert(std::pair<int, int>(itr, currentSolution.getCustomerToDemandMap()[itr]));
						customer.insert(itr);
					}
					solution.push_back(satID);
					capacity += iter.getCapacityServed();
					cost += iter.getRouteCost();
				}
				CVRPSolution cvrpSol(solution, customer, cusToDemand, cost, satID, problemParams.getSecondEchelonVehicleCapacityLimit(), capacity, it.size());
				febSrchSol.insertSecondEchelonSolution(cvrpSol);
			}
			//update first echelon cvrp
			std::vector<int> solution;
			std::map<int, int> cusToDemandMap;
			std::set<int> firstEchelonCustomersSet;
			int firstCapacity = 0;
			double cost = 0;
			int satID = 0;
			for (auto& iter : firstEchelonRoutesList) {
				satID = iter.getSepInt();
				for (auto itr : iter.getRouteSol()) {
					solution.push_back(itr);
					firstEchelonCustomersSet.insert(itr);
				}
				solution.push_back(satID);
				cost += iter.getRouteCost();
			}
			std::list<CVRPSolution> secondSolList = febSrchSol.getSecondEchelonSolutions();
			std::set<int> satNodes = currentSolution.getSatelliteNodes();
			for (auto cus : firstEchelonCustomersSet) {
				auto itr = satNodes.find(cus);
				if (itr != satNodes.end()) {
					for (auto & sat : secondSolList) {
						if (sat.getSatelliteNode() == cus) {
							cusToDemandMap.insert(std::pair<int, int>(cus, sat.getTotalDemandSatisfied()));
							break;
						}
					}
				}
				else {
					cusToDemandMap.insert(std::pair<int, int>(cus, currentSolution.getCustomerToDemandMap()[cus]));
				}
			}
			for (auto & it: cusToDemandMap) {
				firstCapacity += it.second;
			}
			CVRPSolution firstCvrpSol(solution, firstEchelonCustomersSet, cusToDemandMap, cost, satID, problemParams.getFirstEchelonVehicleCapacityLimit(), firstCapacity, firstEchelonRoutesList.size());
			febSrchSol.insertFirstEchelonSolution(firstCvrpSol);
			//update new twoechelon feasible solution
			std::map<int, int> satToDemandMap;
			satToDemandMap.insert(std::pair<int, int>(firstCvrpSol.getSatelliteNode(), firstCvrpSol.getTotalDemandSatisfied()));
			for (auto & it : secondSolList) {
				satToDemandMap.insert(std::pair<int, int>(it.getSatelliteNode(), it.getTotalDemandSatisfied()));
			}
			febSrchSol.populateTwoEchelonSolution(currentSolution.getCustomerToDemandMap(), satToDemandMap, currentSolution.getRoadDistanceMatrix(), currentSolution.getAerialDistanceMatrix(), currentSolution.getCustomerNodes(), currentSolution.getSatelliteNodes(), currentSolution.getCustomersDedicatedToDepot());
			feasibilitySearchSolution = febSrchSol;
		}
	}
	else {
		feasibilitySearchSolution = currentSolution;
		std::cout << "The initial solution is feasible. No need for feasibility search!" << std::endl;
	}
}




