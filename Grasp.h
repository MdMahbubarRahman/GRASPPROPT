#include <iostream>

#include "Initialsolution.h"

#ifndef GRASP_H
#define GRASP_H 

//class 
class CustomerToSatelliteAssignmentProbabilities {
private:
	int customerID;
	std::map<int, double> satIDToProbabilityMap;
public:
	CustomerToSatelliteAssignmentProbabilities();
	CustomerToSatelliteAssignmentProbabilities(int customerID, std::map<int, double> satIDToProbabilityMap);
	CustomerToSatelliteAssignmentProbabilities(const CustomerToSatelliteAssignmentProbabilities & cusToSatAssgnProb);
	int getCustomerID();
	std::map<int, double> getSatIDToProbabilityMap();
	void updateCustomerToSatelliteProbabilityDistribution(int satID, double probability);
};

//class
class AdaptiveCustomersAssignmentProbabilityDistribution {
private:
	std::map<int, CustomerToSatelliteAssignmentProbabilities> customerToSatAssignmentMap;
public:
	AdaptiveCustomersAssignmentProbabilityDistribution();
	AdaptiveCustomersAssignmentProbabilityDistribution(std::map<int, CustomerToSatelliteAssignmentProbabilities> customerToSatAssignmentMap);
	AdaptiveCustomersAssignmentProbabilityDistribution(const AdaptiveCustomersAssignmentProbabilityDistribution & adaptiveCusToSatProbDist);
	std::map<int, CustomerToSatelliteAssignmentProbabilities> getCustomerToSatAssignmentMap();
	void updateAdaptiveCustomersAssignmentProbabilityDistribution(CustomerToSatelliteAssignmentProbabilities cusToSatAssgnProb);
};

//class
class Grasp{
private:
	AdaptiveCustomersAssignmentProbabilityDistribution cusAssignmentDistribution;
	std::map<int, std::set<int>> satelliteToCustomersMap;
	TwoEchelonSolution graspSolution;
	ProblemParameters probParams;
public:
	Grasp();
	Grasp(const Grasp & grasp);
	Grasp(AdaptiveCustomersAssignmentProbabilityDistribution cusAssigmentDistribution, ProblemParameters probParams);
	AdaptiveCustomersAssignmentProbabilityDistribution getCusAssigmentDistribution();
	std::map<int, std::set<int>> getSatelliteToCustomersMap();
	TwoEchelonSolution getGraspSolution();
	void assignCustomersToSatellites();
	void runGrasp();
};

#endif //GRASP_H

