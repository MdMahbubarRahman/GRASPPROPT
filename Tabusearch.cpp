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
		for (const auto& elm : tabuList) {
			std::cout << elm.first << "  " << elm.second << std::endl;
		}
	}
}

//adds tabu move to the tabu list
void TabuList::addToTabuList(int routeNum, int customerID) {
	tabuList.insert(std::make_pair(routeNum, customerID));
}

//return contents of the tabu list
std::multimap<int, int> TabuList::getTabuList(){
	return tabuList;
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
	solution.push_back(0);
}

//constructor
FeasibleSolution::FeasibleSolution(std::vector<int> sol, double cost, int sepIntVal) {
	solCost = cost;
	separatorIntVal = sepIntVal;
	if (!sol.empty())
		for (int i = 0; i < sol.size(); ++i) {
			solution.push_back(sol.at(i));
		}
	else
		std::cout << "The solution vector is empty! Feasible Solution object could not be formed." << std::endl;
}

//copy constructor
FeasibleSolution::FeasibleSolution(const FeasibleSolution & febSol) {
	solCost = febSol.solCost;
	separatorIntVal = febSol.separatorIntVal;
	if (!febSol.solution.empty())
		for (int i = 0; i < febSol.solution.size(); ++i) {
			solution.push_back(febSol.solution.at(i));
		}
	else
		std::cout << "The solution vector is empty! Copy constructor could not copy the solution." << std::endl;
}

//prints the feasible solution object
void FeasibleSolution::showSolution() const {
	if (solution.size() != 0) {
		std::cout << "The solution is : " << std::endl;
		for (const auto& iter : solution) {
			std::cout << iter << ", "; 
		}
		std::cout << std::endl;
		std::cout << "The cost of the solution is : "<< solCost <<std::endl;
		std::cout << "The separtor interger value is : " << separatorIntVal << std::endl;
	}
	else
		std::cout << "The solution is empty!" << std::endl;
}

//returns solution vector
std::vector<int> FeasibleSolution::getSolution() {
	return solution;
}

//returns separator integer value
int FeasibleSolution::getSeparatorIntVal() {
	return separatorIntVal;
}

//returns cost of the solution
double FeasibleSolution::getCost() {
	return solCost;
}

// copy constructor
Neighbourhood::Neighbourhood(const Neighbourhood& neighbour) {
	if (!neighbour.neighbourSolution.empty()) {
		for (auto it : neighbour.neighbourSolution)
			neighbourSolution.push_back(it);
	}
	else
		std::cout << "The neighbourhood solution is empty! The copy constructor fails." << std::endl;

	if (!neighbour.kNeighbourSolution.empty()) {
		for (auto it : neighbour.kNeighbourSolution)
			kNeighbourSolution.push_back(it);
	}
	else
		std::cout << "The K neighbourhood solution is empty! The copy constructor fails." << std::endl;
}

//inserts a new neighbour/solution
void Neighbourhood::insertToNeighbour(FeasibleSolution neighbour) {
	neighbourSolution.push_back(neighbour);
}

//inserts neighbour solutions to the kneighbour solution list
void Neighbourhood::insertToKNeighbour() {
	if (neighbourSolution.size() != 0) {
		for (auto it : neighbourSolution) {
			kNeighbourSolution.push_back((it));
		}
		neighbourSolution.clear();
	}
	else
		std::cout << "The neighbour solution list is empty!" << std::endl;
}

//shows neighbouring solutions
void Neighbourhood::showNeighbours() {
	if (neighbourSolution.size() != 0) {
		for (const auto &it : neighbourSolution) {
			it.showSolution();
		}
	}
}

//shows all k neighbour solutions
void Neighbourhood::showKNeighbours() {
	if (kNeighbourSolution.size() != 0) {
		for (const auto &it : kNeighbourSolution) {
			it.showSolution();
		}
	}
}

//returns the best solution from the neighbourhood
FeasibleSolution Neighbourhood::getBestFromNeighbour() {
	double cost = 10000000.0;
	std::list<FeasibleSolution>::iterator iter;
	if (neighbourSolution.size() != 0) {
		for (auto it = neighbourSolution.begin(); it != neighbourSolution.end(); ++it) {
			if ((*it).getCost() < cost) {
				cost = (*it).getCost();
				iter = it;
			}	
		}
		(*iter).showSolution();
	}
	return (*iter);
}

//returns the best solution from the K neighbourhood
//linear search method implemented, could be better by heap implementation
FeasibleSolution Neighbourhood::getBestFromKNeighbour() {
	double cost = 10000000.0;
	std::list<FeasibleSolution>::iterator iter;
	if (kNeighbourSolution.size() != 0) {
		for (auto it = kNeighbourSolution.begin(); it != kNeighbourSolution.end(); ++it) {
			if ((*it).getCost() < cost) {
				cost = (*it).getCost();
				iter = it;
			}
		}
		(*iter).showSolution();
	}
	return (*iter);
}

//default constructor
Tabusearch::Tabusearch() {}

//copy constructor
Tabusearch::Tabusearch(const Tabusearch& tabusrch) {
	k = 2;
	initialSolution = tabusrch.initialSolution;
	incumbentSolution = tabusrch.incumbentSolution;
	iterationBestSolution = tabusrch.iterationBestSolution;
}

//updates incumbent solution
void Tabusearch::updateIncumbentSolution() {
	if (iterationBestSolution.getCost() < incumbentSolution.getCost())
		incumbentSolution = iterationBestSolution;
}

//generates routes and route to customers map
void Tabusearch::generateRouteCustomerMap(FeasibleSolution febSol) {
	bool flag = false;
	if (!routCustomerMap.empty())
		routCustomerMap.clear();
	if (!listOfRoutes.empty())
		listOfRoutes.clear();
	int routeID = 1;
	int separatorIntVal = febSol.getSeparatorIntVal();
	int size = febSol.getSolution().size();
	std::vector<int> routeVector;
	for (int i = 0; i < size; ++i) {
		int val = febSol.getSolution().at(i);
		if (val != separatorIntVal) {
			routCustomerMap.insert(std::make_pair(routeID, val));
			routeVector.push_back(val);
			flag = true;
		}
		else {
			if (flag) {
				routeID += 1;
				flag = false;
				listOfRoutes.push_back(routeVector);
				routeVector.clear();
			}
		}	
	}
	numberOfRoutes = routeID;
}

//generates Neighbour solutions
//k-chain moves uses Add and Drop as well as 1-Swap methods
void Tabusearch::generateKNeighbourSolutions() {
	int numberOfCustomers = 10;	
	int maxCapacity = 6;
	Neighbourhood neighbourHood;
	//distance distribution
	//this cost matrix needs to be declared in the global scope of the ts
	double distanceMatrix[11][11];
	for (int i = 0; i<11; ++i)
		for (int j = 0; j < 11; ++j) {
			distanceMatrix[i][j] = 10;
		}
	//customer demand distribution
	//this demand array needs to be declared in the global scope of the ts
	int demandVector[11];
	for (int i = 0; i < 11; ++i) {
		if (i == 0)
			demandVector[i] = 0;
		else
			demandVector[i] = 2; //uniform demand
	}
	for (int i = 0; i < k; ++i) {//k for k chain
		//generate route to customer map
		generateRouteCustomerMap(iterationBestSolution);
		//randomly choose routes for Add and Drop heuristic
		std::random_device rd; 
		std::mt19937 gen(rd()); 
		std::uniform_int_distribution<> distr(1, numberOfRoutes); 
		int dropFromRoute = distr(gen);
		int addToRoute = distr(gen);
		if (numberOfRoutes >= 2) {
			while (dropFromRoute == addToRoute) {
				addToRoute = distr(gen);
			}
		}
		else {
			std::cout << "There is only one route in this solution! Hence add and drop heuristic is invalid." << std::endl;
			break;
		}
		//generate neighbour solutions
		for (auto it = routCustomerMap.lower_bound(dropFromRoute); it != routCustomerMap.upper_bound(dropFromRoute); ++it) {
			if (tabuList.getTabuList().size() != 0) {
				//check whether the move is tabu
				if (!tabuList.checkInTabuList(addToRoute, (*it).second)) {
					std::cout << "The move is not tabu" << std::endl;
					//check capacity constraint 
					int capacity = 0;
					for (auto iter = routCustomerMap.lower_bound(addToRoute); iter != routCustomerMap.upper_bound(addToRoute); ++iter) {
						capacity += demandVector[(*iter).second];
					}
					capacity += demandVector[(*it).second];
					if (capacity <= maxCapacity) {
						std::cout << "The capacity constraint is satisfied!" << std::endl;
						//generate neighbour solution by Add and Drop heuristic
						int iterator = 1;
						double currentAddRouteCost = 0.0;
						double newAddRouteCost = 0.0;
						double currentDropRouteCost = 0.0;
						double newDropRouteCost = 0.0;
						std::list<std::vector<int>> newRoutes = listOfRoutes;
						double newCost = iterationBestSolution.getCost();
						int sepInt = iterationBestSolution.getSeparatorIntVal();
						//generate feasible solution
						for (auto iter = newRoutes.begin(); iter != newRoutes.end(); ++iter) {
							if (iterator == addToRoute) {
								//calculate current Add Route cost
								for (int j = 0; j < (*iter).size(); ++j) {
									if (j == 0) {
										currentAddRouteCost += distanceMatrix[0][(*iter).at(j)];
									}
									else {
										currentAddRouteCost += distanceMatrix[(*iter).at(j-1)][(*iter).at(j)];
									}
								}
								currentAddRouteCost += distanceMatrix[(*iter).at((*iter).size()-1)][0];
								//find the minimum cost entry position 
								int l = 0;
								double cost = 0.0;
								for (int j = 0; j < (*iter).size(); ++j) {
									if (j == 0) {
										cost = currentAddRouteCost - distanceMatrix[0][(*iter).at(j)] + distanceMatrix[0][(*it).second] + distanceMatrix[(*it).second][(*iter).at(j)];
										newAddRouteCost = cost;
										l = j;
									}
									else{
										cost = currentAddRouteCost - distanceMatrix[(*iter).at(j-1)][(*iter).at(j)] + distanceMatrix[(*iter).at(j-1)][(*it).second] + distanceMatrix[(*it).second][(*iter).at(j)];
										if (cost < newAddRouteCost) {
											newAddRouteCost = cost;
											l = j;
										}
									}								
								}
								cost = currentAddRouteCost - distanceMatrix[(*iter).at((*iter).size())][0] + distanceMatrix[(*iter).at((*iter).size())][(*it).second] + distanceMatrix[(*it).second][0];
								if (cost < newAddRouteCost) {
									newAddRouteCost = cost;
									l = (*iter).size();
								}
								//update Add route  
								auto& elm = *iter;
								std::list<int> oldRoute;
								for (int m = 0; m < (*iter).size(); ++m) {
									oldRoute.push_back((*iter).at(m));
								}
								if (l == (*iter).size()) {
									oldRoute.push_back((*it).second);
								}
								else {
									int n = 0;
									for (auto oit = oldRoute.begin(); oit != oldRoute.end(); ++oit) {
										if (n == l) {
											oldRoute.insert(oit, (*it).second);
											break;
										}
										n++;
									}
								}
								elm.clear();
								for (auto uit = oldRoute.begin(); uit != oldRoute.end(); ++uit)
									elm.push_back(*uit);			
							}
							if (iterator == dropFromRoute) {
								int dropNode = (*it).second;
								//calculate current Drop Route cost
								for (int j = 0; j < (*iter).size(); ++j) {
									if (j == 0) {
										currentDropRouteCost += distanceMatrix[0][(*iter).at(j)];
									}
									else {
										currentDropRouteCost += distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)];
									}
								}
								currentDropRouteCost += distanceMatrix[(*iter).at((*iter).size() - 1)][0];
								//calculate updated Drop Route cost
								int l = 0;
								for (int j = 0; j < (*iter).size(); ++j) {
									if ((*iter).size()==1) {
										newDropRouteCost = 0.0;
										break;
									}
									else {
										if (j == 0 && (*iter).at(j) == dropNode) {
											newDropRouteCost = currentDropRouteCost - distanceMatrix[0][(*iter).at(j)] - distanceMatrix[(*iter).at(j)][(*iter).at(j + 1)] + distanceMatrix[0][(*iter).at(j + 1)];
											l = j;
											break;
										}
										else if (j == ((*iter).size() - 1) && (*iter).at(j) == dropNode) {
											newDropRouteCost = currentDropRouteCost - distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)] - distanceMatrix[(*iter).at(j)][0] + distanceMatrix[(*iter).at(j - 1)][0];
											l = j;
											break;
										}
										else if ((*iter).at(j) == dropNode) {
											newDropRouteCost = currentDropRouteCost - distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)] - distanceMatrix[(*iter).at(j)][(*iter).at(j + 1)] + distanceMatrix[(*iter).at(j - 1)][(*iter).at(j + 1)];
											l = j;
											break;
										}
										else
											continue;
									}
								}
								//update Drop route  
								auto& elm = *iter;
								std::list<int> oldRoute;
								for (int m = 0; m < (*iter).size(); ++m) {
									oldRoute.push_back((*iter).at(m));
								}
								int n = 0;
								for (auto oit = oldRoute.begin(); oit != oldRoute.end(); ++oit) {
									if (n == l) {
										oldRoute.erase(oit);
										break;
									}
									n++;
								}
								elm.clear();
								if (!oldRoute.empty()) {
									for (auto uit = oldRoute.begin(); uit != oldRoute.end(); ++uit)
										elm.push_back(*uit);
								}
								else {
									elm.push_back(iterationBestSolution.getSeparatorIntVal());//temporary! 
								}	
							}
							iterator += 1;
						}
						//generate solution vector
						std::vector<int> solution;
						for (auto iter = newRoutes.begin(); iter != newRoutes.end(); ++iter) {
							for (int p = 0; p < (*iter).size(); ++p) {
								solution.push_back((*iter).at(p));
							}
							solution.push_back(sepInt);
						}
						//update solution cost
						newCost = newCost - currentAddRouteCost - currentAddRouteCost + newAddRouteCost + newDropRouteCost;
						//new neighbour
						FeasibleSolution neighbour(solution, newCost, sepInt);
						neighbourHood.insertToNeighbour(neighbour);		
					}
					else {
						std::cout << "The capacity constraint is violated!" << std::endl;
						continue;
					}
				}
				else {
					std::cout << "The move is tabu!" << std::endl;
					continue;
				}
			}
			else {
				std::cout << "The tabu list is empty!" << std::endl;
				//check capacity constraint
				int capacity = 0;
				for (auto iter = routCustomerMap.lower_bound(addToRoute); iter != routCustomerMap.upper_bound(addToRoute); ++iter) {
					capacity += demandVector[(*iter).second];
				}
				capacity += demandVector[(*it).second];
				if (capacity <= maxCapacity) {
					std::cout << "The capacity constraint is satisfied!" << std::endl;
					//generate neighbour solution by Add and Drop heuristic
					int iterator = 1;
					double currentAddRouteCost = 0.0;
					double newAddRouteCost = 0.0;
					double currentDropRouteCost = 0.0;
					double newDropRouteCost = 0.0;
					std::list<std::vector<int>> newRoutes = listOfRoutes;
					double newCost = iterationBestSolution.getCost();
					int sepInt = iterationBestSolution.getSeparatorIntVal();
					//generate feasible solution
					for (auto iter = newRoutes.begin(); iter != newRoutes.end(); ++iter) {
						if (iterator == addToRoute) {
							//calculate current Add Route cost
							for (int j = 0; j < (*iter).size(); ++j) {
								if (j == 0) {
									currentAddRouteCost += distanceMatrix[0][(*iter).at(j)];
								}
								else {
									currentAddRouteCost += distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)];
								}
							}
							currentAddRouteCost += distanceMatrix[(*iter).at((*iter).size() - 1)][0];
							//find the minimum cost entry position 
							int l = 0;
							double cost = 0.0;
							for (int j = 0; j < (*iter).size(); ++j) {
								if (j == 0) {
									cost = currentAddRouteCost - distanceMatrix[0][(*iter).at(j)] + distanceMatrix[0][(*it).second] + distanceMatrix[(*it).second][(*iter).at(j)];
									newAddRouteCost = cost;
									l = j;
								}
								else {
									cost = currentAddRouteCost - distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)] + distanceMatrix[(*iter).at(j - 1)][(*it).second] + distanceMatrix[(*it).second][(*iter).at(j)];
									if (cost < newAddRouteCost) {
										newAddRouteCost = cost;
										l = j;
									}
								}
							}
							cost = currentAddRouteCost - distanceMatrix[(*iter).at((*iter).size())][0] + distanceMatrix[(*iter).at((*iter).size())][(*it).second] + distanceMatrix[(*it).second][0];
							if (cost < newAddRouteCost) {
								newAddRouteCost = cost;
								l = (*iter).size();
							}
							//update Add route  
							auto& elm = *iter;
							std::list<int> oldRoute;
							for (int m = 0; m < (*iter).size(); ++m) {
								oldRoute.push_back((*iter).at(m));
							}
							if (l == (*iter).size()) {
								oldRoute.push_back((*it).second);
							}
							else {
								int n = 0;
								for (auto oit = oldRoute.begin(); oit != oldRoute.end(); ++oit) {
									if (n == l) {
										oldRoute.insert(oit, (*it).second);
										break;
									}
									n++;
								}
							}
							elm.clear();
							for (auto uit = oldRoute.begin(); uit != oldRoute.end(); ++uit)
								elm.push_back(*uit);
						}
						if (iterator == dropFromRoute) {
							int dropNode = (*it).second;
							//calculate current Drop Route cost
							for (int j = 0; j < (*iter).size(); ++j) {
								if (j == 0) {
									currentDropRouteCost += distanceMatrix[0][(*iter).at(j)];
								}
								else {
									currentDropRouteCost += distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)];
								}
							}
							currentDropRouteCost += distanceMatrix[(*iter).at((*iter).size() - 1)][0];
							//calculate updated Drop Route cost
							int l = 0;
							for (int j = 0; j < (*iter).size(); ++j) {
								if ((*iter).size() == 1) {
									newDropRouteCost = 0.0;
									break;
								}
								else {
									if (j == 0 && (*iter).at(j) == dropNode) {
										newDropRouteCost = currentDropRouteCost - distanceMatrix[0][(*iter).at(j)] - distanceMatrix[(*iter).at(j)][(*iter).at(j + 1)] + distanceMatrix[0][(*iter).at(j + 1)];
										l = j;
										break;
									}
									else if (j == ((*iter).size() - 1) && (*iter).at(j) == dropNode) {
										newDropRouteCost = currentDropRouteCost - distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)] - distanceMatrix[(*iter).at(j)][0] + distanceMatrix[(*iter).at(j - 1)][0];
										l = j;
										break;
									}
									else if ((*iter).at(j) == dropNode) {
										newDropRouteCost = currentDropRouteCost - distanceMatrix[(*iter).at(j - 1)][(*iter).at(j)] - distanceMatrix[(*iter).at(j)][(*iter).at(j + 1)] + distanceMatrix[(*iter).at(j - 1)][(*iter).at(j + 1)];
										l = j;
										break;
									}
									else
										continue;
								}
							}
							//update Drop route  
							auto& elm = *iter;
							std::list<int> oldRoute;
							for (int m = 0; m < (*iter).size(); ++m) {
								oldRoute.push_back((*iter).at(m));
							}
							int n = 0;
							for (auto oit = oldRoute.begin(); oit != oldRoute.end(); ++oit) {
								if (n == l) {
									oldRoute.erase(oit);
									break;
								}
								n++;
							}
							elm.clear();
							if (!oldRoute.empty()) {
								for (auto uit = oldRoute.begin(); uit != oldRoute.end(); ++uit)
									elm.push_back(*uit);
							}
							else {
								elm.push_back(iterationBestSolution.getSeparatorIntVal());//temporary! 
							}
						}
						iterator += 1;
					}
					//generate solution vector
					std::vector<int> solution;
					for (auto iter = newRoutes.begin(); iter != newRoutes.end(); ++iter) {
						for (int p = 0; p < (*iter).size(); ++p) {
							solution.push_back((*iter).at(p));
						}
						solution.push_back(sepInt);
					}
					//update solution cost
					newCost = newCost - currentAddRouteCost - currentAddRouteCost + newAddRouteCost + newDropRouteCost;
					//new neighbour
					FeasibleSolution neighbour(solution, newCost, sepInt);
					neighbourHood.insertToNeighbour(neighbour);
				}
			}
		}
		neighbourHood.insertToKNeighbour();
	}
	iterationBestSolution = neighbourHood.getBestFromKNeighbour();
}

//runs tabu search algoirthm
void Tabusearch::tabuSearchRun(FeasibleSolution febSol) {
	int MaxIter = 10;
	TabuList tabuList;
	AspirationCriteria aspCria;
	initialSolution = febSol;
	for (int i = 0; i < MaxIter; ++i) {
		generateKNeighbourSolutions();
		updateIncumbentSolution();
	}
}

//performs TSP heuristics on incumbent solution
void Tabusearch::performTspHeuristics() {

}

