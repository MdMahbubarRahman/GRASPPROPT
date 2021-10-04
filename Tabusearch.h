#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <random>
#include <array>
#include <chrono>

using namespace std::chrono;
#ifndef TABUSEARCH_H
#define TABUSEARCH_H

/*
TODO:Based on "EVE-OPT: a hybrid algorithm for the capacitated vehicle routing problem" 
paper implement following functions and algorithm:
	1. Add_and_drop function
	2. 1-Swap function
	3. N neighbour solution generating function 
	4. Decent heuristic algorithm
	5. k-chain-moves algorithm
	6. Tabusearch algorithm

More insight to develop above functions and algoirthms can be found from following papers:
	 1. Cyclic transfer algorithm for multivehicle routing and scheduling problems by
	 PM Thompson, HN Psaraftis.
	 2. TS2PACK: A two-level tabu search for the three-dimensional bin packing problem by
	 Teodor Gabriel Crainic a, Guido Perboli b,*, Roberto Tadei b
	 3. The Theory of Cyclic Transfers by
	 Paul M. Thompson and James B. Orlin
*/



/*
Bblue-print of the Tabu search algorithm

INPUT: A Feasible solution

Start: randomly choose a target vehicle in order to generate its neighbourhood

Generate neighbourhood: using Add and Drop and 1-swap neighbourhood method
K-chain-method: improve the accuracy of the neighbourhood using k-chain method
using k-chain move neighbourhood will be generated where both add and drop and 1-swap method will
be used.

Set of Solutions: check for best solution from the set of solutuions

TSP Heuristics: improve the best solution using tsp heuristics

Stoping criteria: aspiration criteria and stopping criteria

OUTPUT: A Feasible solution

Tabu List: it is a global list which will
track which customer is removed from which vehicle/route.

Aspiration criteria: this criteria will override tabu action plans
Stopping criteria: number of iteration elapsed without improvement.
*/

//class tabulist 
class TabuList {
private:
	std::multimap<int, int> tabuList;
public:
	TabuList(); 
	TabuList(const TabuList & tbList);
	TabuList(int routeNum, int customerID);
	void showTabuList() const;
	void addToTabuList(int routeNum, int customerID);
	std::multimap<int, int> getTabuList();
	bool checkInTabuList(int routeNumber, int customerID);  
};


//class aspiration criteria
class AspirationCriteria {
private:
	double currentBestCost;
public:
	AspirationCriteria();
	AspirationCriteria(const AspirationCriteria & aspCriteria);
	void updateCurrentBestCost(double cost);
	void showCurrentBestCost();
	double getCurrentBestCost();
};


//class feasible solution 
class FeasibleSolution {
private:
	std::list<int> solution;
	double solCost = 0;
	int separatorIntVal = 0;
	int sourceRoute = 0;
	int addedNode = 0;
public:
	FeasibleSolution();//default constructor
	FeasibleSolution(std::list<int> sol, double cost, int sepIntVal, int sourceRoute, int addedNode);//constructor
	FeasibleSolution(const FeasibleSolution & sol); //copy constructor
	void showSolution() const;
	std::list<int> getSolution();
	int getSeparatorIntVal();
	double getCost();
	int getSourceRoute();
	int getAddedNode();
};


//class neighbourhood
class Neighbourhood {
private:
	std::list<FeasibleSolution> neighbourSolution;
	std::list<FeasibleSolution> kNeighbourSolution;
public:
	Neighbourhood() {}//default constructor
	Neighbourhood(const Neighbourhood & neighbour);// copy constructor
	void insertToNeighbour(FeasibleSolution febSolution);
	void insertToKNeighbour();
	void showNeighbours();
	void showKNeighbours();
	std::list<FeasibleSolution> getNeighbourSolutions();
	std::list<FeasibleSolution> getKNeighbourSolutions();
	FeasibleSolution getBestFromNeighbour();
	FeasibleSolution getBestFromKNeighbour();
};


class Tabusearch{
private:
	int kChain;// k for k-chain
	int swapChain;
	int dropFromRoute;
	int addToRoute;
	int numberOfRoutes;
	int maxRouteCapacity;
	std::vector<int> demandVector;
	std::vector<std::vector<double>> distanceMatrix;
	TabuList tabuList;
	FeasibleSolution initialSolution;//input solution to the tabu search algoirthm
	FeasibleSolution incumbentSolution;//current best solution
	FeasibleSolution iterationBestSolution;//output solution from the tabu search algoirthm
	std::multimap<int, int> routCustomerMap;
	std::list<std::list<int>> listOfRoutes;
	Neighbourhood neighbourHood;
	AspirationCriteria aspCriteria;
public:	
	Tabusearch();//default constructor
	Tabusearch(FeasibleSolution febSol, std::vector<int> demandVec, std::vector<std::vector<double>> disMat, int kChain, int swapChain, int maxRouteCapacity);//constructor
	Tabusearch(const Tabusearch& tabusrch);//copy constructor
	void updateIncumbentSolution();
	void generateRouteCustomerMap(FeasibleSolution febSol);
	void selectRandomAddAndDropRoutes();
	FeasibleSolution generateNeighbourByAddDrop(std::list<std::list<int>> listOfRoutes, double cost, int separatorInt, int addToRoute, int dropFromRoute, int dropNode);
	FeasibleSolution generateNeighbourByOneSwap(std::list<std::list<int>> listOfRoutes, double cost, int separatorInt, int firstRoute, int secondRoute);
	void generateKChainNeighbourSolutions();//k-chain moves uses Add and Drop method
	void generateOneSwapSolutions();
	void runTabuSearch();//This is the function that will perform the tabu search
	void performTspHeuristics();//perform this to improve the incumbent solution.
	FeasibleSolution getInitialSolution();
	FeasibleSolution getIncumbentSolution();
	FeasibleSolution getIterationBestSolution();
	std::multimap<int, int> getRoutCustomerMap();
	std::list<std::list<int>> getListOfRoutes();
	TabuList getTabuList();
	int getNumberOfRoutes();
	Neighbourhood getNeighbourHood();
	AspirationCriteria getAspirationCriteria();
	void showTabuSolution();
};


#endif // !TABUSEARCH_H




