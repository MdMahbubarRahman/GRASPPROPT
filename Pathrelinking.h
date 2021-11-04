#include <iostream>
#include <queue>
#include <list>
#include <vector>

#include "Initialsolution.h"

#ifndef PATHRELINKING_H
#define PATHRELINKING_H 


/*
TODO: Develop Path Relinking algoirthm based on following papers:
	1. GRASP with Path Relinking for the Two-Echelon Vehicle Routing Problem.
	2. The two-echelon vehicle routing problem by S Mancini 
*/

//customer status class
class CustomerStatus {
private:
	int cusID;
	int satIDInLocalSol;
	int satIDInBestSol;
	double distanceDiff;
public:
	CustomerStatus();
	CustomerStatus(const CustomerStatus & cusStat);
	CustomerStatus(int cusID, int satIDInLocalSol, int satIDInBestSol, double distanceDiff);
	int getCusID();
	int getLocalSatID();
	int getBestSatID();
	double getDistanceDiff();
	void showCustomerStatus();
};

//comparator
class SortTool {
public:
	bool operator()(CustomerStatus &a , CustomerStatus &b);
};

//path relinking class
class Pathrelinking {
private:
	TwoEchelonSolution localSolution;
	TwoEchelonSolution bestSolution;
	TwoEchelonSolution intermediateSolution;
	TwoEchelonSolution currentBestSolution;
	std::priority_queue<CustomerStatus, std::vector<CustomerStatus>, SortTool> customersToBeReassigned;
public:
	Pathrelinking();
	Pathrelinking(const Pathrelinking & pathRelkn);
	Pathrelinking(TwoEchelonSolution localSol, TwoEchelonSolution bestSol);
	TwoEchelonSolution getLocalSolution();
	TwoEchelonSolution getInitialBestSolution();
	TwoEchelonSolution getIntermediateSolution();
	TwoEchelonSolution getCurrentBestSolution();
	void showLocalSolution();
	void showBestSolution();
	void showIntermediateSolution();
	void runPathRelinking();
	void showCurrentBestSolution();
};


#endif // !PATHRELINKING_H



