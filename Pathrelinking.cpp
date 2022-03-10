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
	currentBestSolution = pathRelkn.currentBestSolution;
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
	return currentBestSolution;
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
	intermediateSolution = localSolution;
	currentBestSolution = localSolution;
	std::map<int, int> localCusToSatMap;
	std::map<int, int> bestCusToSatMap;
	std::set<int> localFirstEchelonCustomers = localSolution.getFirstEchelonSolution().getCustomers();
	std::set<int> bestFirstEchelonCustomers = bestSolution.getFirstEchelonSolution().getCustomers();
	std::list<CVRPSolution> localSecondEchelonSolutions = localSolution.getSecondEchelonSolutions();
	std::list<CVRPSolution> bestSecondEchelonSolutions = bestSolution.getSecondEchelonSolutions();
	int localDepot = localSolution.getFirstEchelonSolution().getSatelliteNode();
	int bestDepot = bestSolution.getFirstEchelonSolution().getSatelliteNode();
	//populate customer to satellite map for both solutions
	for (auto it: localFirstEchelonCustomers) {
		localCusToSatMap.insert(std::pair<int, int>(it, localDepot));
	}
	for (auto it : bestFirstEchelonCustomers) {
		bestCusToSatMap.insert(std::pair<int, int>(it, bestDepot));
	}
	for (auto & it: localSecondEchelonSolutions) {
		int satID = it.getSatelliteNode();
		std::set<int> customers = it.getCustomers();
		for (auto cus: customers) {
			localCusToSatMap.insert(std::pair<int, int>(cus, satID));
		}
	}
	for (auto & it : bestSecondEchelonSolutions) {
		int satID = it.getSatelliteNode();
		std::set<int> customers = it.getCustomers();
		for (auto cus : customers) {
			bestCusToSatMap.insert(std::pair<int, int>(cus, satID));
		}
	}
	//perform customers reassignment
	std::set<int> allCustomers = localSolution.getCustomerNodes();
	for (auto it : allCustomers) {
		if (localCusToSatMap[it] != bestCusToSatMap[it]) {
			double bestDist = 0;
			double localDist = 0;
			//cus to sat distance for local solution
			if (localCusToSatMap[it] == localDepot) {
				localDist = localSolution.getRoadDistanceMatrix()[it][localDepot];
			}
			else {
				localDist = localSolution.getAerialDistanceMatrix()[it][localCusToSatMap[it]];
			}
			//cus to sat distance for best solution
			if (bestCusToSatMap[it] == bestDepot) {
				bestDist = bestSolution.getRoadDistanceMatrix()[it][bestDepot];
			}
			else {
				bestDist = bestSolution.getAerialDistanceMatrix()[it][bestCusToSatMap[it]];
			}
			double dist = abs(bestDist - localDist);
			CustomerStatus cs(it, localCusToSatMap[it], bestCusToSatMap[it], dist);
			customersToBeReassigned.push(cs);
		}
	}
	std::cout << "The size of the priority queue is : " << customersToBeReassigned.size() << std::endl;
	//perform path relinking and optimization
	int maxPathIter = 0;
	while (!customersToBeReassigned.empty()) {
		CustomerStatus cusStatus = customersToBeReassigned.top();
		customersToBeReassigned.pop();
		int cusID = cusStatus.getCusID();
		int dropFromSat = cusStatus.getLocalSatID();
		int addToSat = cusStatus.getBestSatID();
		double diffCost = cusStatus.getDistanceDiff();
		CVRPSolution customerToDropFromCVRPSol;
		CVRPSolution customerToAddToCVRPSol;
		CVRPSolution firstEchelonCVRPSol = intermediateSolution.getFirstEchelonSolution();
		std::list<CVRPSolution> secondEchelonCVRPSols = intermediateSolution.getSecondEchelonSolutions();
		if (dropFromSat != firstEchelonCVRPSol.getSatelliteNode()) {
			for (auto & it : secondEchelonCVRPSols) {
				if (it.getSatelliteNode() == dropFromSat) {
					customerToDropFromCVRPSol = it;
				}
			}
			if (addToSat != firstEchelonCVRPSol.getSatelliteNode()) {
				for (auto & it : secondEchelonCVRPSols) {
					if (it.getSatelliteNode() == addToSat)
						customerToAddToCVRPSol = it;
				}
			}
			else {
				customerToAddToCVRPSol = firstEchelonCVRPSol;
			}
		}
		else {
			customerToDropFromCVRPSol = firstEchelonCVRPSol;
			for (auto & it : secondEchelonCVRPSols) {
				if (it.getSatelliteNode() == addToSat) {
					customerToAddToCVRPSol = it;
				}
			}
		}		
		//update addTo satellite
		CVRPSolution updatedAddToCVRPSol;
		if (addToSat != firstEchelonCVRPSol.getSatelliteNode()) {
			int demand = intermediateSolution.getCustomerToDemandMap()[cusID];
			std::map<int, int> cusToDemandMapUp = customerToAddToCVRPSol.getCustomerTodemandMap();
			cusToDemandMapUp.insert(std::pair<int, int>(cusID, demand));
			std::vector<int> cusCluster;
			for (auto it : customerToAddToCVRPSol.getCustomers()) {
				cusCluster.push_back(it);
			}
			cusCluster.push_back(cusID);
			//need to check whether added capacity makes the satellite infeasible?
			Geneticalgorithm gap(addToSat, customerToAddToCVRPSol.getMaxRouteCapacity(), cusToDemandMapUp, intermediateSolution.getAerialDistanceMatrix(), cusCluster);
			gap.runGeneticAlgorithm();
			Chromosome chromp = gap.getGASolution();
			std::set<int> customersp;
			for (auto it : cusCluster) {
				customersp.insert(it);
			}
			updatedAddToCVRPSol = CVRPSolution(chromp, customersp, cusToDemandMapUp);
		}
		//update dropFrom satellite
		CVRPSolution updatedDropFromCVRPSol;
		if (dropFromSat != firstEchelonCVRPSol.getSatelliteNode()) {
			std::map<int, int> cusToDemandMapCur = customerToDropFromCVRPSol.getCustomerTodemandMap();
			std::map<int, int>::iterator iter;
			for (auto it = cusToDemandMapCur.begin(); it != cusToDemandMapCur.end(); ++it) {
				if ((*it).first == cusID) {
					iter = it;
					break;
				}
			}
			cusToDemandMapCur.erase(iter);
			std::vector<int> cusClusterCur;
			for (auto it : customerToDropFromCVRPSol.getCustomers()) {
				if (it != cusID) {
					cusClusterCur.push_back(it);
				}
			}
			Geneticalgorithm gac(dropFromSat, customerToDropFromCVRPSol.getMaxRouteCapacity(), cusToDemandMapCur, intermediateSolution.getAerialDistanceMatrix(), cusClusterCur);
			gac.runGeneticAlgorithm();
			Chromosome chromc = gac.getGASolution();
			std::set<int> customersc;
			for (auto it : cusClusterCur) {
				customersc.insert(it);
			}
			updatedDropFromCVRPSol = CVRPSolution(chromc, customersc, cusToDemandMapCur);
		}
		//update first echelon solution
		int depot = firstEchelonCVRPSol.getSatelliteNode();
		CVRPSolution updatedFirstEchelonCVRPSol;
		if (depot == dropFromSat) {
			int addToSatDemand = updatedAddToCVRPSol.getTotalDemandSatisfied();
			std::map<int, int> cusToDemandMap = firstEchelonCVRPSol.getCustomerTodemandMap();
			std::map<int, int>::iterator iterp, iterd;
			for (auto it = cusToDemandMap.begin(); it != cusToDemandMap.end(); ++it) {
				if ((*it).first == addToSat) {
					iterp = it;
				}
				if ((*it).first == cusID) {
					iterd = it;
				}
			}
			cusToDemandMap.erase(iterp);
			cusToDemandMap.erase(iterd);
			if (addToSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(addToSat, addToSatDemand));
			}
			std::vector<int> customersVec;
			for (auto it : firstEchelonCVRPSol.getCustomers()) {
				if (it != cusID) {
					customersVec.push_back(it);
				}
			}
			Geneticalgorithm gad(depot, firstEchelonCVRPSol.getMaxRouteCapacity(), cusToDemandMap, intermediateSolution.getRoadDistanceMatrix(), customersVec);
			gad.runGeneticAlgorithm();
			Chromosome chromd = gad.getGASolution();
			std::set<int> customersd;
			for (auto it : customersVec) {
				customersd.insert(it);
			}
			updatedFirstEchelonCVRPSol = CVRPSolution(chromd, customersd, cusToDemandMap);
			updatedDropFromCVRPSol = CVRPSolution(chromd, customersd, cusToDemandMap);
		}
		else if (depot == addToSat) {
			int dropFromSatDemand = updatedDropFromCVRPSol.getTotalDemandSatisfied();
			std::map<int, int> cusToDemandMap = firstEchelonCVRPSol.getCustomerTodemandMap();
			std::map<int, int>::iterator iterc;
			for (auto it = cusToDemandMap.begin(); it != cusToDemandMap.end(); ++it) {
				if ((*it).first == dropFromSat) {
					iterc = it;
				}
			}
			cusToDemandMap.erase(iterc);
			if (dropFromSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(dropFromSat, dropFromSatDemand));
			}
			std::vector<int> customersVec;
			for (auto it : firstEchelonCVRPSol.getCustomers()) {
				customersVec.push_back(it);
			}
			customersVec.push_back(cusID);
			int cusDemand = intermediateSolution.getCustomerToDemandMap()[cusID];
			cusToDemandMap.insert(std::pair<int, int>(cusID, cusDemand));
			Geneticalgorithm gad(depot, firstEchelonCVRPSol.getMaxRouteCapacity(), cusToDemandMap, intermediateSolution.getRoadDistanceMatrix(), customersVec);
			gad.runGeneticAlgorithm();
			Chromosome chromd = gad.getGASolution();
			std::set<int> customersd;
			for (auto it : customersVec) {
				customersd.insert(it);
			}
			updatedFirstEchelonCVRPSol = CVRPSolution(chromd, customersd, cusToDemandMap);
			updatedAddToCVRPSol = CVRPSolution(chromd, customersd, cusToDemandMap);
		}
		else {
			int dropFromSatDemand = updatedDropFromCVRPSol.getTotalDemandSatisfied();
			int addToSatDemand = updatedAddToCVRPSol.getTotalDemandSatisfied();
			std::map<int, int> cusToDemandMap = firstEchelonCVRPSol.getCustomerTodemandMap();
			std::map<int, int>::iterator iterc, iterp;
			for (auto it = cusToDemandMap.begin(); it != cusToDemandMap.end(); ++it) {
				if ((*it).first == dropFromSat) {
					iterc = it;
				}
				if ((*it).first == addToSat) {
					iterp = it;
				}
			}
			cusToDemandMap.erase(iterc);
			cusToDemandMap.erase(iterp);
			std::vector<int> customersVec;
			for (auto it : firstEchelonCVRPSol.getCustomers()) {
				if (it == dropFromSat || it == addToSat) {
					continue;
				}
				else {
					customersVec.push_back(it);
				}
			}
			if (dropFromSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(dropFromSat, dropFromSatDemand));
				customersVec.push_back(dropFromSat);
			}
			if (addToSatDemand != 0) {
				cusToDemandMap.insert(std::pair<int, int>(addToSat, addToSatDemand));
				customersVec.push_back(addToSat);
			}
			Geneticalgorithm gad(depot, firstEchelonCVRPSol.getMaxRouteCapacity(), cusToDemandMap, intermediateSolution.getRoadDistanceMatrix(), customersVec);
			gad.runGeneticAlgorithm();
			Chromosome chromd = gad.getGASolution();
			std::set<int> customersd;
			for (auto it : customersVec) {
				customersd.insert(it);
			}
			updatedFirstEchelonCVRPSol = CVRPSolution(chromd, customersd, cusToDemandMap);
		}
		//update twoechelon solution
		TwoEchelonSolution newTwoEchelonSol;
		newTwoEchelonSol.insertFirstEchelonSolution(updatedFirstEchelonCVRPSol);
		std::map<int, int> satToDemandMap;
		satToDemandMap.insert(std::pair<int, int>(updatedFirstEchelonCVRPSol.getSatelliteNode(), updatedFirstEchelonCVRPSol.getTotalDemandSatisfied()));
		for (auto & it : secondEchelonCVRPSols) {
			if (it.getSatelliteNode() == dropFromSat) {
				newTwoEchelonSol.insertSecondEchelonSolution(updatedDropFromCVRPSol);
				satToDemandMap.insert(std::pair<int, int>(it.getSatelliteNode(), updatedDropFromCVRPSol.getTotalDemandSatisfied()));
			}
			else if (it.getSatelliteNode() == addToSat) {
				newTwoEchelonSol.insertSecondEchelonSolution(updatedAddToCVRPSol);
				satToDemandMap.insert(std::pair<int, int>(it.getSatelliteNode(), updatedAddToCVRPSol.getTotalDemandSatisfied()));
			}
			else {
				newTwoEchelonSol.insertSecondEchelonSolution(it);
				satToDemandMap.insert(std::pair<int, int>(it.getSatelliteNode(), it.getTotalDemandSatisfied()));
			}
		}
		newTwoEchelonSol.populateTwoEchelonSolution(intermediateSolution.getCustomerToDemandMap(), satToDemandMap, intermediateSolution.getRoadDistanceMatrix(), intermediateSolution.getAerialDistanceMatrix(), intermediateSolution.getCustomerNodes(), intermediateSolution.getSatelliteNodes(), intermediateSolution.getCustomersDedicatedToDepot());
		//check with terminating conditions and/or update local solution
		if (newTwoEchelonSol.getSolutionFitness() < intermediateSolution.getSolutionFitness()) {
			currentBestSolution = newTwoEchelonSol;
			intermediateSolution = newTwoEchelonSol;
			maxPathIter += 1;
			if (maxPathIter >= 2) {
				break;
			}
			//currentBestSolution.showTwoEchelonSolution();
			//break;
		}
	}
}

//prints current best solution
void Pathrelinking::showCurrentBestSolution() {
	currentBestSolution.showTwoEchelonSolution();
}



