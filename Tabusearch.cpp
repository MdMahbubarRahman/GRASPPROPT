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
	if (!tabuList.empty()) {
		auto itlow = tabuList.lower_bound(routeNum);
		auto itup = tabuList.upper_bound(routeNum);
		if (itlow != tabuList.end()) {
			for (auto& it = itlow; it != itup; it++) {
				if ((*it).second == customerID) {
					flag = true;
					break;
				}

			}
		}
	}
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
FeasibleSolution::FeasibleSolution(std::list<int> sol, double cost, int sepIntVal, int srcRoute, int addNode){
	solCost			 = cost;
	separatorIntVal  = sepIntVal;
	sourceRoute		 = srcRoute;
	addedNode		 = addNode;
	if (!sol.empty())
		for (auto it = sol.begin(); it != sol.end(); ++it) {
			solution.push_back(*it);
		}
	else {
		std::cout << "The solution vector is empty! Feasible Solution object could not be formed." << std::endl;
	}	
}

//copy constructor
FeasibleSolution::FeasibleSolution(const FeasibleSolution & febSol) {
	solCost			= febSol.solCost;
	separatorIntVal = febSol.separatorIntVal;
	sourceRoute	    = febSol.sourceRoute;
	addedNode	    = febSol.addedNode;
	if (!febSol.solution.empty())
		for (auto it = febSol.solution.begin(); it != febSol.solution.end(); ++it) {
			solution.push_back(*it);
		}
	else {
		std::cout << "The solution vector is empty! Feasible Solution object could not be formed." << std::endl;
	}
}

//prints the feasible solution object
void FeasibleSolution::showSolution() const {
	if (!solution.empty()) {
		std::cout << "The solution is : " << std::endl;
		for (const auto& iter : solution) {
			std::cout << iter << ", "; 
		}
		std::cout << std::endl;
		std::cout << "The cost of the solution is : "<< solCost <<std::endl;
		std::cout << "The separtor interger value is : " << separatorIntVal << std::endl;
		std::cout << "The added Node is : " << addedNode << std::endl;
		std::cout << "The node is added from Route : " << sourceRoute << std::endl;
	}
	else
		std::cout << "The solution is empty!" << std::endl;
}

//returns solution vector
std::list<int> FeasibleSolution::getSolution() {
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

//return the outsourced route
int FeasibleSolution::getSourceRoute() {
	return sourceRoute;
}
int FeasibleSolution::getAddedNode() {
	return addedNode;
}

// copy constructor
Neighbourhood::Neighbourhood(const Neighbourhood& neighbour) {
	Neighbourhood neighHood;
	if (!neighbour.neighbourSolution.empty()) {
		for (auto it : neighbour.neighbourSolution)
			neighHood.neighbourSolution.push_back(it);
	}
	else
		std::cout << "The neighbourhood solution is empty! The copy constructor fails." << std::endl;

	if (!neighbour.kNeighbourSolution.empty()) {
		for (auto it : neighbour.kNeighbourSolution)
			neighHood.kNeighbourSolution.push_back(it);
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
	if (!neighbourSolution.empty()) {
		for (auto it : neighbourSolution) {
			kNeighbourSolution.push_back((it));
		}
		neighbourSolution.clear();
	}
	else {
		std::cout << "The neighbour solution list is empty!" << std::endl;
	}	
}

//shows neighbouring solutions
void Neighbourhood::showNeighbours() {
	if (!neighbourSolution.empty()) {
		for (const auto &it : neighbourSolution) {
			it.showSolution();
		}
	}
}

//shows all k neighbour solutions
void Neighbourhood::showKNeighbours() {
	if (!kNeighbourSolution.empty()) {
		for (const auto &it : kNeighbourSolution) {
			it.showSolution();
		}
	}
}

//returns the best solution from the neighbourhood
FeasibleSolution Neighbourhood::getBestFromNeighbour() {
	double cost = 10000000.0;
	std::list<FeasibleSolution>::iterator iter;
	if (!neighbourSolution.empty()) {
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
	if (!kNeighbourSolution.empty()) {
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
Tabusearch::Tabusearch(FeasibleSolution febSol, std::vector<int> demandVec, std::vector<std::vector<double>> disMat) {
	k = 2;
	numberOfRoutes = 0;
	maxRouteCapacity = 10;

	for (int i = 0; i < demandVec.size(); ++i) {
		demandVector.push_back(demandVec.at(i));
	}
	for (double i = 0; i < disMat.size(); ++i) {
		distanceMatrix.push_back(disMat.at(i));
	}
}

//copy constructor
//this one has potential problem
Tabusearch::Tabusearch(const Tabusearch& tabusrch) {
	k = tabusrch.k;
	numberOfRoutes = tabusrch.numberOfRoutes;
	initialSolution = tabusrch.initialSolution;
	incumbentSolution = tabusrch.incumbentSolution;
	iterationBestSolution = tabusrch.iterationBestSolution;
	distanceMatrix = tabusrch.distanceMatrix;
	demandVector = tabusrch.demandVector;
	tabuList = tabusrch.tabuList;
	routCustomerMap = tabusrch.routCustomerMap;
	listOfRoutes = tabusrch.listOfRoutes;
}

//updates incumbent solution
void Tabusearch::updateIncumbentSolution() {
	if (iterationBestSolution.getCost() < incumbentSolution.getCost())
		incumbentSolution = iterationBestSolution;
}

//generates routes and route to customers map
void Tabusearch::generateRouteCustomerMap(FeasibleSolution febSol) {
	if (!routCustomerMap.empty())
		routCustomerMap.clear();
	if (!listOfRoutes.empty())
		listOfRoutes.clear();
	int routeID = 1;
	int separatorIntVal = febSol.getSeparatorIntVal();
	int size = febSol.getSolution().size();
	std::list<int> routeVector;
	bool flag = false;
	std::list<int> solution = febSol.getSolution();
	for (auto it = solution.begin(); it != solution.end(); ++it) {
		int val = *it;
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

//generate neighbour solution by Add and Drop heuristic
FeasibleSolution Tabusearch::generateNeighbourByAddDrop(std::list<std::list<int>> newRoutes, double newCost, int sepInt, int addToRoute, int dropFromRoute, int dropNode) {
	double currentAddRouteCost = 0.0;
	double newAddRouteCost = 0.0;
	double currentDropRouteCost = 0.0;
	double newDropRouteCost = 0.0;
	int iterator = 1;
	for (auto &iter : newRoutes) {
		if (iterator == addToRoute) {
			//calculate current Add Route cost
			int val = 0;
			for (auto it = iter.begin(); it != iter.end(); ++it) {
				currentAddRouteCost += distanceMatrix[val][(*it)];
				val = *it;
			}
			currentAddRouteCost += distanceMatrix[val][0];
			//find the minimum cost entry position 
			int l = 0;
			int j = 0;
			double cost = 0.0;
			int vall = 0;
			for (auto it = iter.begin(); it != iter.end(); ++it) {
				if (j == 0) {
					cost = currentAddRouteCost - distanceMatrix[0][(*it)] + distanceMatrix[0][dropNode] + distanceMatrix[dropNode][(*it)];
					newAddRouteCost = cost;
					l = j;
					vall = *it;
					j++;
				}
				else {
					cost = currentAddRouteCost - distanceMatrix[vall][(*it)] + distanceMatrix[vall][dropNode] + distanceMatrix[dropNode][(*it)];
					if (cost < newAddRouteCost) {
						newAddRouteCost = cost;
						l = j;
					}
					j++;
					vall = *it;
				}
			}
			cost = currentAddRouteCost - distanceMatrix[vall][0] + distanceMatrix[vall][dropNode] + distanceMatrix[dropNode][0];
			if (cost < newAddRouteCost) {
				newAddRouteCost = cost;
				l = iter.size();
			}
			//update Add route  
			int n = 0;
			if(l==iter.size()){
				iter.push_back(dropNode);			
			}
			else {
				for (auto it = iter.begin(); it != iter.end(); ++it) {
					if (n == l) {
						iter.insert(it, dropNode);
						break;
					}
					n++;
				}
			}
		}
		if (iterator == dropFromRoute) {
			//calculate current Drop Route cost
			int val = 0;
			for (auto it = iter.begin(); it != iter.end(); ++it) {
				currentDropRouteCost += distanceMatrix[val][(*it)];
				val = *it;
			}
			currentDropRouteCost += distanceMatrix[val][0];
			//calculate updated Drop Route cost
			int l = 0;
			int j = 0;
			int vall = 0;
			int preVal = 0;
			int postVal = 0;
			bool flag = false;
			for (auto it = iter.begin(); it != iter.end(); ++it) {
				if (iter.size() == 1) {
					newDropRouteCost = 0.0;
					iter.erase(it);
					iter.push_back(sepInt+1);
					break;
				}
				else {
					if (j == 0 && flag == false) {
						vall = *it;
					}
					if (j >= 1 && flag == false) {
						preVal = vall;
						vall = *it;
					}
					if ((*it) == dropNode && flag == false) {
						flag = true;
						l = j;
					}
					if ((*it) != dropNode && flag == true) {
						postVal = *it;
						break;
					}
					j++;
				}					
			}
			if (iter.size() != 1) {
				newDropRouteCost = currentDropRouteCost - distanceMatrix[preVal][vall] - distanceMatrix[vall][postVal] + distanceMatrix[preVal][postVal];
			}
			//update Drop route  
			if (flag != false) {
				int n = 0;
				for (auto it = iter.begin(); it != iter.end(); ++it) {
					if (n == l) {
						iter.erase(it);
						break;
					}
					n++;
				}
			}
		}
		iterator += 1;
	}
	//generate solution vector
	std::list<int> solution;
	for (auto iter = newRoutes.begin(); iter != newRoutes.end(); ++iter) {
		for (auto it = (*iter).begin();  it != (*iter).end(); ++it) {
			solution.push_back(*it);
		}
		solution.push_back(sepInt);
	}
	newCost = newCost - currentAddRouteCost - currentDropRouteCost + newAddRouteCost + newDropRouteCost;
	FeasibleSolution neighbour(solution, newCost, sepInt, dropFromRoute, dropNode);
	return neighbour;
}

//associating cost with each route would speed up the computation
FeasibleSolution Tabusearch::generateNeighbourByOneSwap(std::list<std::list<int>> newRoutes, double newCost, int sepInt, int firstRoute, int secondRoute) {
	double currentFirstRouteCost = 0.0;
	double currentSecondRouteCost = 0.0;
	//new route vectors
	std::vector<int> firstRt;
	std::vector<int> secondRt;
	int iterator = 1;
	for (auto &iter : newRoutes) {
		if (iterator == firstRoute) {
			int val = 0;
			for (auto it = iter.begin(); it != iter.end(); ++it) {
				currentFirstRouteCost += distanceMatrix[val][(*it)];
				val = *it;
				firstRt.push_back(val);
			}
			currentFirstRouteCost += distanceMatrix[val][0];			
		}
		if (iterator == secondRoute) {
			int val = 0;
			for (auto it = iter.begin(); it != iter.end(); ++it) {
				currentSecondRouteCost += distanceMatrix[val][(*it)];
				val = *it;
				secondRt.push_back(val);
			}
			currentSecondRouteCost += distanceMatrix[val][0];			
		}
		iterator += 1;
	}
	//find out minimum cost swap pair
	double firstCost = 0;
	double secondCost = 0;
	double cost = 0;
	struct SwapInfo {
		int firstNodeLocation = 0;
		int secondNodeLocation = 0;
		int firstNode = 0;
		int secondNode = 0;
		double finalCost = 0;
		double newFirstRouteCost = 0.0;
		double newSecondRouteCost = 0.0;
		void updateSwapInfo(int fstNodeL, int sndNodeL, int fstNode, int sndNode, double finlCost, double newFstRtCost, double newSndRtCost) {
			firstNodeLocation = fstNodeL;
			secondNodeLocation = sndNodeL;
			firstNode = fstNode;
			secondNode = sndNode;
			finalCost = finlCost;
			newFirstRouteCost = newFstRtCost;
			newSecondRouteCost = newSndRtCost;
		}
	};
	SwapInfo swapInfo;
	swapInfo.finalCost = currentFirstRouteCost + currentSecondRouteCost;
	int firstPreVal = 0;
	int firstPostVal = 0;
	int secondPreVal = 0;
	int secondPostVal = 0;
	//n^2 computational complexity
	for (int i = 0; i < firstRt.size(); ++i) {
		for (int j = 0; j < secondRt.size(); ++j) {
			if (firstRt.size() ==1) {
				firstPreVal = 0;
				firstPostVal = 0;
			}
			if (secondRt.size() == 1) {
				secondPreVal = 0;
				secondPostVal = 0;
			}
			if (firstRt.size() > 1 && i==0) {
				firstPreVal = 0;
				firstPostVal = firstRt.at(i+1);
			}
			if (firstRt.size() > 1 && i == (firstRt.size()-1)) {
				firstPreVal = firstRt.at(i - 1);
				firstPostVal = 0;
			}
			if (firstRt.size() > 1 && i < (firstRt.size() - 1) && i > 0) {
				firstPreVal = firstRt.at(i - 1);
				firstPostVal = firstRt.at(i + 1);
			}
			if (secondRt.size() > 1 && j == 0) {
				secondPreVal = 0;
				secondPostVal = secondRt.at(j + 1);
			}
			if (secondRt.size() > 1 && j == (secondRt.size() - 1)) {
				secondPreVal = secondRt.at(j - 1);
				secondPostVal = 0;
			}
			if (secondRt.size() > 1 && j < (secondRt.size() - 1) && j > 0) {
				secondPreVal = secondRt.at(j - 1);
				secondPostVal = secondRt.at(j + 1);
			}
			firstCost = currentFirstRouteCost - double(distanceMatrix[firstPreVal][firstRt.at(i)]) - double(distanceMatrix[firstRt.at(i)][firstPostVal]) + double(distanceMatrix[firstPreVal][secondRt.at(j)]) + double(distanceMatrix[secondRt.at(j)][firstPostVal]);
			secondCost = currentSecondRouteCost - double(distanceMatrix[secondPreVal][secondRt.at(j)]) - double(distanceMatrix[secondRt.at(j)][secondPostVal]) + double(distanceMatrix[secondPreVal][firstRt.at(i)]) + double(distanceMatrix[firstRt.at(i)][secondPostVal]);
			cost = firstCost + secondCost;
			if (cost < swapInfo.finalCost) {
				swapInfo.updateSwapInfo(i, j, firstRt.at(i), secondRt.at(j), cost, firstCost, secondCost);
			}
		}
	}
	//update routes and generate new solution
	std::list<int> solution;
	iterator = 1;
	for (auto &iter : newRoutes) {
		if (iterator == firstRoute) {
			int i = 0;
			for (auto& it : iter) {
				if (i == swapInfo.firstNodeLocation && swapInfo.firstNode == it) {
					it = swapInfo.secondNode;
				}
				solution.push_back(it);
				i++;
			}
		}
		else if (iterator == secondRoute) {
			int j = 0;
			for (auto& it : iter) {
				if (j == swapInfo.secondNodeLocation && swapInfo.secondNode == it) {
					it = swapInfo.firstNode;
				}
				solution.push_back(it);
				j++;
			}
		}
		else {
			for (auto& it : iter) {
				solution.push_back(it);
			}
		} 
		solution.push_back(sepInt);
		iterator += 1;
	}
	newCost = newCost - currentFirstRouteCost - currentSecondRouteCost + swapInfo.newFirstRouteCost + swapInfo.newSecondRouteCost;
	FeasibleSolution neighbour(solution, newCost, sepInt, firstRoute, swapInfo.firstNode);
	return neighbour;
}




//if possible use reference to avoid creating local copy!
//generates Neighbour solutions
void Tabusearch::generateKChainNeighbourSolutions() {
	Neighbourhood neighbourHood;
	for (int i = 0; i < k; ++i) {//k for k chain
		//generate route to customer map
		generateRouteCustomerMap(iterationBestSolution);
		//randomly choose routes for Add and Drop heuristic
		std::random_device rd; 
		std::mt19937 gen(rd()); 
		std::uniform_int_distribution<> distr(1, numberOfRoutes); 
		int dropFromRoute = distr(gen);//need to exclude null routes
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
					if (capacity <= maxRouteCapacity) {
						std::cout << "The capacity constraint is satisfied!" << std::endl;
						//generate neighbour solution by Add and Drop heuristic
						std::list<std::list<int>> newRoutes = listOfRoutes;
						double newCost = iterationBestSolution.getCost();
						int sepInt = iterationBestSolution.getSeparatorIntVal();
						int dropNode = (*it).second;
						FeasibleSolution neighbour = generateNeighbourByAddDrop(newRoutes, newCost, sepInt, addToRoute, dropFromRoute, dropNode);
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
				if (capacity <= maxRouteCapacity) {
					std::cout << "The capacity constraint is satisfied!" << std::endl;
					//generate neighbour solution by Add and Drop heuristic
					std::list<std::list<int>> newRoutes = listOfRoutes;
					double newCost = iterationBestSolution.getCost();
					int sepInt = iterationBestSolution.getSeparatorIntVal();
					int dropNode = (*it).second;
					FeasibleSolution neighbour = generateNeighbourByAddDrop(newRoutes, newCost, sepInt, addToRoute, dropFromRoute, dropNode);
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
		generateKChainNeighbourSolutions();
		updateIncumbentSolution();
	}
}

//performs TSP heuristics on incumbent solution
void Tabusearch::performTspHeuristics() {

}

