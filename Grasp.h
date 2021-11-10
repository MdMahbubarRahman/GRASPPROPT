#include <iostream>

#include "Initialsolution.h"

#ifndef GRASP_H
#define GRASP_H 

//class 
class CustomerToSatelliteAssignmentProbabilities {
private:
	int customerID;
	int firstChoiceSatID;
	double highestProb;
	std::map<int, double> satIDToProbabilityMap;
public:
	CustomerToSatelliteAssignmentProbabilities();
	CustomerToSatelliteAssignmentProbabilities(int customerID, std::map<int, double> satIDToProbabilityMap);
	CustomerToSatelliteAssignmentProbabilities(const CustomerToSatelliteAssignmentProbabilities & cusToSatAssgnProb);
	int getCustomerID();
	int getFirstChoiceSatID();
	double getHighestProb();
	std::map<int, double> getSatIDToProbabilityMap();
	void updateCustomerToSatelliteProbabilityDistribution(int satID, double probability);
	void showCusToSatAssngProbs();
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
	void populateAdaptiveCusAssgnProbDist(ProblemParameters probParams);
};
 
//class 
class compareHighestProbs {
public:
	bool operator()(CustomerToSatelliteAssignmentProbabilities & a, CustomerToSatelliteAssignmentProbabilities & b);
};

//class
class Grasp{
private:
	Chromosome chrom;
	CVRPSolution cvrpSol;
	ProblemParameters probParams;
	TwoEchelonSolution graspSolution;
	std::map<int, int> satelliteToDemandMap;
	std::multimap<int, int> satelliteToCustomersMap;
	std::priority_queue<CustomerToSatelliteAssignmentProbabilities, std::vector<CustomerToSatelliteAssignmentProbabilities>, compareHighestProbs> cusAssignQueue;
	AdaptiveCustomersAssignmentProbabilityDistribution cusAssignmentDistribution;	
public:
	Grasp();
	Grasp(const Grasp & grasp);
	Grasp(AdaptiveCustomersAssignmentProbabilityDistribution cusAssigmentDistribution, ProblemParameters probParams);
	AdaptiveCustomersAssignmentProbabilityDistribution getCusAssigmentDistribution();
	std::multimap<int, int> getSatelliteToCustomersMap();
	TwoEchelonSolution getGraspSolution();
	void makeCustomerToSatelliteAssignment();
	void generateCVRPSolutions();
	void generateTwoEchelonSolution();
	void runGrasp();
};

#endif //GRASP_H



