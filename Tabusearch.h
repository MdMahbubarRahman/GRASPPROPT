#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <random>
#include <array>


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
	void updateCurrentBestCost(double cost);
	double showCurrentBestCost();
};


//class feasible solution 
class FeasibleSolution {
private:
	std::vector<int> solution;
	double solCost = 0;
	int separatorIntVal = 0;
public:
	FeasibleSolution();//default constructor
	FeasibleSolution(std::vector<int> sol, double cost, int sepIntVal);//constructor
	FeasibleSolution(const FeasibleSolution & sol); //copy constructor
	void showSolution() const;
	std::vector<int> getSolution();
	int getSeparatorIntVal();
	double getCost();
};


//class neighbourhood
class Neighbourhood {
private:
	std::list<FeasibleSolution> neighbourSolution;
	std::list<FeasibleSolution> kNeighbourSolution;
public:
	Neighbourhood() {}//default constructor
	Neighbourhood(const Neighbourhood & neighbour);// copy constructor
	//Neighbourhood operator=(const Neighbourhood & neighHood) const;//operator= overloading
	void insertToNeighbour(FeasibleSolution febSolution);
	void insertToKNeighbour();//call of this function makes the neighbourSolution list empty
	void showNeighbours();
	void showKNeighbours();
	FeasibleSolution getBestFromNeighbour();
	FeasibleSolution getBestFromKNeighbour();
};


class Tabusearch{
private:
	int k=2;// k for k-chain
	FeasibleSolution initialSolution;//input solution to the tabu search algoirthm
	FeasibleSolution incumbentSolution;//current best solution
	FeasibleSolution iterationBestSolution;//output solution from the tabu search algoirthm
	std::multimap<int, int> routCustomerMap;
	int numberOfRoutes = 0;
	TabuList tabuList;
	std::list<std::vector<int>> listOfRoutes;
	double distanceMatrix[100][100];
	int demandVector[100];
public:	
	Tabusearch();//default constructor
	Tabusearch(const Tabusearch& tabusrch);//copy constructor
	void updateIncumbentSolution();
	void generateRouteCustomerMap(FeasibleSolution febSol);
	FeasibleSolution generateNeighbour(std::list<std::vector<int>> listOfRoutes, double cost, int separatorInt, int addToRoute, int dropFromRoute, int dropNode);
	void generateKNeighbourSolutions();//k-chain moves uses Add and Drop method
	void tabuSearchRun(FeasibleSolution febSol);//This is the function that will perform the tabu search
	void performTspHeuristics();//perform this to improve the incumbent solution.
};


#endif // !TABUSEARCH_H




