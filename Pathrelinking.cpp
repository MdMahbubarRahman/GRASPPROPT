#include "Pathrelinking.h"

//Implement Path Relinking algoirthm here

//default constructor
CustomerStatus::CustomerStatus() {
	cusID = 0;
	satIDInBestSol = 0;
	satIDInLocalSol = 0;
	distanceDiff = 0;
}

//copy constructor
CustomerStatus::CustomerStatus(const CustomerStatus& cusStat) {
	cusID = cusStat.cusID;
	satIDInBestSol = cusStat.satIDInBestSol;
	satIDInLocalSol = cusStat.satIDInLocalSol;
	distanceDiff = cusStat.distanceDiff;
}

//constructor
CustomerStatus::CustomerStatus(int cID, int satIDLocalSol, int satIDBestSol, double distDiff) {
	cusID = cID;
	satIDInBestSol = satIDBestSol;
	satIDInLocalSol = satIDLocalSol;
	distanceDiff = distDiff;
}

//returns customer id
int CustomerStatus::getCusID() {
	return cusID;
}

//returns satellite id for the local solution
int CustomerStatus::getLocalSatID() {
	return satIDInLocalSol;
}

//returns satellite id for the best solution
int CustomerStatus::getBestSatID() {
	return satIDInBestSol;
}

//returns relative distance difference
double CustomerStatus::getDistanceDiff() {
	return distanceDiff;
}

//prints customer's status
void CustomerStatus::showCustomerStatus() {
	std::cout << "\nShow customer ID : " << cusID << std::endl;
	std::cout << "Satellite ID in the local solution : " << satIDInLocalSol << std::endl;
	std::cout << "Satellite ID in the best solution : " << satIDInBestSol << std::endl;
	std::cout << "Relative distance difference : " << distanceDiff << std::endl;
}

//comparator
bool SortTool::operator()(CustomerStatus& a, CustomerStatus& b) {
	return (a.getDistanceDiff() > b.getDistanceDiff());
}

//default constructor
Pathrelinking::Pathrelinking() {
	std::cout << "The default constructor of the path relinking class has been called!" << std::endl;
}

//copy constructor
Pathrelinking::Pathrelinking(const Pathrelinking& pathRelkn) {
	localSolution = pathRelkn.localSolution;
	bestSolution = pathRelkn.bestSolution;
	intermediateSolution = pathRelkn.intermediateSolution;
	customersToBeReassigned = pathRelkn.customersToBeReassigned;
}

//constructor
Pathrelinking::Pathrelinking(TwoEchelonSolution localSol, TwoEchelonSolution bestSol) {
	localSolution = localSol;
	bestSolution = bestSol;
}

//returns local solution
TwoEchelonSolution Pathrelinking::getLocalSolution() {
	return localSolution;
}

//returns best solution
TwoEchelonSolution Pathrelinking::getInitialBestSolution() {
	return bestSolution;
}

//returns intermediate solution
TwoEchelonSolution Pathrelinking::getIntermediateSolution() {
	return intermediateSolution;
}

//returns best solution
TwoEchelonSolution Pathrelinking::getCurrentBestSolution() {
	runPathRelinking();
	return bestSolution;
}

//prints local solution
void Pathrelinking::showLocalSolution() {
	localSolution.showTwoEchelonSolution();
}

//prints best solution
void Pathrelinking::showBestSolution() {
	bestSolution.showTwoEchelonSolution();
}

//prints intermediate solution
void Pathrelinking::showIntermediateSolution() {
	intermediateSolution.showTwoEchelonSolution();
}

//runs the path relinking algorithm
void Pathrelinking::runPathRelinking() {
	//needs implementation
}



