#include "Grasp.h"

//implement grasp here

//default constructor
CustomerToSatelliteAssignmentProbabilities::CustomerToSatelliteAssignmentProbabilities() {
	customerID = 0;
	firstChoiceSatID = 0;
	highestProb = 0;
}

//constructor
CustomerToSatelliteAssignmentProbabilities::CustomerToSatelliteAssignmentProbabilities(int cusID, std::map<int, double> satIDToProbMap) {
	customerID = cusID;
	satIDToProbabilityMap = satIDToProbMap;
	double prob = 0;
	int sat = 0;
	for (auto & it : satIDToProbabilityMap) {
		if (it.second > prob) {
			prob = it.second;
			sat = it.first;
		}
	}
	firstChoiceSatID = sat;
	highestProb = prob;
}

//copy constructor
CustomerToSatelliteAssignmentProbabilities::CustomerToSatelliteAssignmentProbabilities(const CustomerToSatelliteAssignmentProbabilities& cusToSatAssgnProb) {
	customerID = cusToSatAssgnProb.customerID;
	satIDToProbabilityMap = cusToSatAssgnProb.satIDToProbabilityMap;
	firstChoiceSatID = cusToSatAssgnProb.firstChoiceSatID;
	highestProb = cusToSatAssgnProb.highestProb;
}

//returns cusID
int CustomerToSatelliteAssignmentProbabilities::getCustomerID() {
	return customerID;
}

//returns first choice sat id
int CustomerToSatelliteAssignmentProbabilities::getFirstChoiceSatID() {
	return firstChoiceSatID;
}

//returns prob measure
double CustomerToSatelliteAssignmentProbabilities::getHighestProb() {
	return highestProb;
}

//returns assignment probability distribution map for each customers
std::map<int, double> CustomerToSatelliteAssignmentProbabilities::getSatIDToProbabilityMap() {
	return satIDToProbabilityMap;
}

//updates probability
void CustomerToSatelliteAssignmentProbabilities::updateCustomerToSatelliteProbabilityDistribution(int satID, double probability) {
	satIDToProbabilityMap.erase(satID);
	satIDToProbabilityMap.insert(std::pair<int, double>(satID, probability));
	double prob = 0;
	int sat = 0;
	for (auto& it : satIDToProbabilityMap) {
		if (it.second > prob) {
			prob = it.second;
			sat = it.first;
		}
	}
	firstChoiceSatID = sat;
	highestProb = prob;
}

//show assignment probs
void CustomerToSatelliteAssignmentProbabilities::showCusToSatAssngProbs() {
	std::cout << "Customer ID : " << customerID << std::endl;
	std::cout << "First choice satID : " << firstChoiceSatID << std::endl;
	std::cout << "Highest probability : " << highestProb << std::endl;
	for (auto & it : satIDToProbabilityMap) {
		std::cout << "sat ID : " << it.first << ", probability : " << it.second << std::endl;
	}
}

//default constructor
AdaptiveCustomersAssignmentProbabilityDistribution::AdaptiveCustomersAssignmentProbabilityDistribution() {
	std::cout << "default constructor has been called!" << std::endl;
}

//constructor
AdaptiveCustomersAssignmentProbabilityDistribution::AdaptiveCustomersAssignmentProbabilityDistribution(std::map<int, CustomerToSatelliteAssignmentProbabilities> cusToSatAssignmentMap) {
	customerToSatAssignmentMap = cusToSatAssignmentMap;
}

//copy constructor
AdaptiveCustomersAssignmentProbabilityDistribution::AdaptiveCustomersAssignmentProbabilityDistribution(const AdaptiveCustomersAssignmentProbabilityDistribution& adaptiveCusToSatProbDist) {
	customerToSatAssignmentMap = adaptiveCusToSatProbDist.customerToSatAssignmentMap;
}

//retuns customer to assignment probability distribution map
std::map<int, CustomerToSatelliteAssignmentProbabilities> AdaptiveCustomersAssignmentProbabilityDistribution::getCustomerToSatAssignmentMap() {
	return customerToSatAssignmentMap;
}

//updates
void AdaptiveCustomersAssignmentProbabilityDistribution::updateAdaptiveCustomersAssignmentProbabilityDistribution(CustomerToSatelliteAssignmentProbabilities cusToSatAssgnProb) {
	int cusID = cusToSatAssgnProb.getCustomerID();
	customerToSatAssignmentMap.erase(cusID);
	customerToSatAssignmentMap.insert(std::pair<int, CustomerToSatelliteAssignmentProbabilities>(cusID, cusToSatAssgnProb));
}

//populates initial customers' assignment probability distribution
void AdaptiveCustomersAssignmentProbabilityDistribution::populateAdaptiveCusAssgnProbDist(ProblemParameters probParams) {
	std::vector<std::vector<double>> distanceMat = probParams.getDistanceMatrix();
	std::set<int> customers = probParams.getCustomerNodes();
	std::set<int> satellites = probParams.getSatelliteNodes();
	std::set<int> dedicatedCusForFirstEchelon = probParams.getCustomersMustServeByFirstEchelon();
	customerToSatAssignmentMap.clear();
	//populate the map
	for (auto it : customers) {
		double totalDistance = 0;
		std::map<int, double> satToProbMap;
		for (auto itt : satellites) {
			totalDistance += distanceMat[it][itt];
		}
		//check dedicated condition
		auto iter = dedicatedCusForFirstEchelon.find(it);
		if (iter == dedicatedCusForFirstEchelon.end()) {
			for (auto itt : satellites) {
				double prob = (1 - (distanceMat[it][itt] / totalDistance)) / (satellites.size());
				satToProbMap.insert(std::pair<int, double>(itt, prob));
			}
		}
		else {
			for (auto itt : satellites) {
				if (itt == 0) {//depot is 0 by default which could be changed later
					satToProbMap.insert(std::pair<int, double>(itt, 1.0));
				}
				else {
					satToProbMap.insert(std::pair<int, double>(itt, 0.0));
				}
			}
		}
		CustomerToSatelliteAssignmentProbabilities cusToSatProb(it, satToProbMap);
		customerToSatAssignmentMap.insert(std::pair<int, CustomerToSatelliteAssignmentProbabilities>(it, cusToSatProb));
	}
	//show the probs
	/*
	for (auto &it : customerToSatAssignmentMap) {
		std::cout << "the customer id : " << it.first << std::endl;
		for (auto & itr : it.second.getSatIDToProbabilityMap()) {
			std::cout << "sat id : " << itr.first << " assignment prob : " << itr.second << std::endl;
		}
	}
	*/
}

//comparator
bool compareHighestProbs::operator()(CustomerToSatelliteAssignmentProbabilities& a, CustomerToSatelliteAssignmentProbabilities& b) {
	return (a.getHighestProb() < b.getHighestProb());
}

//default grasp constructor
Grasp::Grasp() {
	std::cout << "Grasp default constructor has been called!" << std::endl;
}

//constructor
Grasp::Grasp(AdaptiveCustomersAssignmentProbabilityDistribution cusAssgnDistribution, ProblemParameters probPrm) {
	cusAssignmentDistribution = cusAssgnDistribution;
	probParams = probPrm;
	for (auto it: probParams.getSatelliteNodes()) {
		satelliteToDemandMap.insert(std::pair<int, int>(it, 0));
	}
}

//copy constructor
Grasp::Grasp(const Grasp& grasp) {
	cusAssignmentDistribution = grasp.cusAssignmentDistribution;
	probParams = grasp.probParams;
	graspSolution = grasp.graspSolution;
	satelliteToCustomersMap = grasp.satelliteToCustomersMap;
	satelliteToDemandMap = grasp.satelliteToDemandMap;
}

//returns customers to satellite assignment probability distribution
AdaptiveCustomersAssignmentProbabilityDistribution Grasp::getCusAssigmentDistribution() {
	return cusAssignmentDistribution;
}

//returns satellite to customers' set map
std::multimap<int, int> Grasp::getSatelliteToCustomersMap() {
	return satelliteToCustomersMap;
}

//defines customers cluster to each satellite
void Grasp::makeCustomerToSatelliteAssignment() {
	//populate the priority queue
	for (auto & it : cusAssignmentDistribution.getCustomerToSatAssignmentMap()) {
		cusAssignQueue.push(it.second);
	}
	//assign customers to satellites
	while (!cusAssignQueue.empty()) {
		CustomerToSatelliteAssignmentProbabilities abc = cusAssignQueue.top();
		cusAssignQueue.pop();
		int cusID = abc.getCustomerID();
		int satID = abc.getFirstChoiceSatID();
		if (satID == 0) {
			satelliteToCustomersMap.insert(std::make_pair(satID, cusID));
			int demand = satelliteToDemandMap[satID]+probParams.getCustomerToDemandMap()[cusID];
			satelliteToDemandMap.erase(satID);
			satelliteToDemandMap.insert(std::pair<int, int>(satID, demand));
		}
		else {
			int demand = satelliteToDemandMap[satID];
			if ((demand + probParams.getCustomerToDemandMap()[cusID]) <= (probParams.getMaxNumberOfVehicleInSecondEchelon()*probParams.getSecondEchelonVehicleCapacityLimit())) {
				demand += probParams.getCustomerToDemandMap()[cusID];
				satelliteToDemandMap.erase(satID);
				satelliteToDemandMap.insert(std::pair<int, int>(satID, demand));
				satelliteToCustomersMap.insert(std::make_pair(satID, cusID));
			}
			else {
				abc.updateCustomerToSatelliteProbabilityDistribution(satID, 0.0);//default 0.0
				cusAssignmentDistribution.updateAdaptiveCustomersAssignmentProbabilityDistribution(abc);
				cusAssignQueue.push(abc);
			}
		}
		std::cout << "\n" << std::endl;
	}
	//show th assignment per satellite
	/*
	std::cout << "show satellite to customer assignment" << std::endl;
	for (auto & it: satelliteToCustomersMap) {
		std::cout << "sat ID : " << it.first << ", customer id : " << it.second << std::endl;
	}
	*/
	int depotDemand = satelliteToDemandMap[0];
	for (auto it: probParams.getSatelliteNodes()) {
		if (it != 0) {
			depotDemand += satelliteToDemandMap[it];
		}
	}
	satelliteToDemandMap.erase(0);
	satelliteToDemandMap.insert(std::pair<int, int>(0, depotDemand));
	/*
	std::cout << "show satellite to demand map" << std::endl;
	for (auto& it : satelliteToDemandMap) {
		std::cout << "sat ID : " << it.first << " demand : " << it.second << std::endl;
	}
	*/
	for (auto it : probParams.getSatelliteNodes()) {
		if (it != 0) {
			satelliteToCustomersMap.insert(std::make_pair(0, it));
		}
	}
}

//generates cvrp solution
void Grasp::generateCVRPSolutions() {
	std::vector<std::vector<double>> distance = probParams.getDistanceMatrix();
	int capLimitFirst = probParams.getFirstEchelonVehicleCapacityLimit();
	int capLimitSecond = probParams.getSecondEchelonVehicleCapacityLimit();
	std::map<int, int> cusToDemand = probParams.getCustomerToDemandMap();
	//update cusToDemand map
	for (auto& it : satelliteToDemandMap) {
		if (it.first != 0) {
			cusToDemand.insert(it);
		}
	}
	//solve cvrp problem
	graspSolution.clearSecondEchelonSolutionList();
	for (auto iter : probParams.getSatelliteNodes()) {
		int satID = iter;
		std::vector<int> customerSet;
		for (auto it = satelliteToCustomersMap.lower_bound(iter); it != satelliteToCustomersMap.upper_bound(iter); ++it) {
			customerSet.push_back((*it).second);
		}
		int capLimit = 0;
		satID == 0 ? capLimit = capLimitFirst : capLimit = capLimitSecond;
		std::cout << "\nsatellite : " << satID << " capacity limit for a route : " << capLimit << std::endl;
		std::cout << "\nThe customers assigned to this satellite are : " << std::endl;
		for (auto it : customerSet) {
			std::cout << "customer : " << it << " demand : " << cusToDemand[it] << std::endl;
		}
		Geneticalgorithm ga(satID, capLimit, cusToDemand, distance, customerSet);
		ga.runGeneticAlgorithm();
		chrom = ga.getGASolution();
		std::set<int> customers;
		for (auto it : customerSet) {
			customers.insert(it);
		}
		cvrpSol = CVRPSolution(chrom, customers, cusToDemand);
		satID == 0 ? graspSolution.insertFirstEchelonSolution(cvrpSol) : graspSolution.insertSecondEchelonSolution(cvrpSol);
	}
	
}

//generates two echelon solution
void Grasp::generateTwoEchelonSolution() {
	graspSolution.populateTwoEchelonSolution(probParams.getCustomerToDemandMap(), satelliteToDemandMap, probParams.getDistanceMatrix(), probParams.getCustomerNodes(), probParams.getSatelliteNodes());
}


//runs grasp procedure
void Grasp::runGrasp() {
	makeCustomerToSatelliteAssignment();
	generateCVRPSolutions();
	generateTwoEchelonSolution();
	graspSolution.showTwoEchelonSolution();
}

//returns two echelon solution 
TwoEchelonSolution Grasp::getGraspSolution() {
	return graspSolution;
}


