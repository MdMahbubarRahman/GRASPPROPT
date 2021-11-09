#include "Grasp.h"

//implement grasp here

//default constructor
CustomerToSatelliteAssignmentProbabilities::CustomerToSatelliteAssignmentProbabilities() {
	customerID = 0;
}

//constructor
CustomerToSatelliteAssignmentProbabilities::CustomerToSatelliteAssignmentProbabilities(int cusID, std::map<int, double> satIDToProbMap) {
	customerID = cusID;
	satIDToProbabilityMap = satIDToProbMap;
}

//copy constructor
CustomerToSatelliteAssignmentProbabilities::CustomerToSatelliteAssignmentProbabilities(const CustomerToSatelliteAssignmentProbabilities& cusToSatAssgnProb) {
	customerID = cusToSatAssgnProb.customerID;
	satIDToProbabilityMap = cusToSatAssgnProb.satIDToProbabilityMap;
}

//returns cusID
int CustomerToSatelliteAssignmentProbabilities::getCustomerID() {
	return customerID;
}

//returns assignment probability distribution map for each customers
std::map<int, double> CustomerToSatelliteAssignmentProbabilities::getSatIDToProbabilityMap() {
	return satIDToProbabilityMap;
}

//updates probability
void CustomerToSatelliteAssignmentProbabilities::updateCustomerToSatelliteProbabilityDistribution(int satID, double probability) {
	satIDToProbabilityMap.erase(satID);
	satIDToProbabilityMap.insert(std::pair<int, double>(satID, probability));
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

//default grasp constructor
Grasp::Grasp() {
	std::cout << "Grasp default constructor has been called!" << std::endl;
}

//constructor
Grasp::Grasp(AdaptiveCustomersAssignmentProbabilityDistribution cusAssgnDistribution, ProblemParameters probPrm) {
	cusAssignmentDistribution = cusAssgnDistribution;
	probParams = probPrm;
}

//copy constructor
Grasp::Grasp(const Grasp& grasp) {
	cusAssignmentDistribution = grasp.cusAssignmentDistribution;
	probParams = grasp.probParams;
	graspSolution = grasp.graspSolution;
	satelliteToCustomersMap = grasp.satelliteToCustomersMap;
}

//returns customers to satellite assignment probability distribution
AdaptiveCustomersAssignmentProbabilityDistribution Grasp::getCusAssigmentDistribution() {
	return cusAssignmentDistribution;
}

//returns satellite to customers' set map
std::map<int, std::set<int>> Grasp::getSatelliteToCustomersMap() {
	return satelliteToCustomersMap;
}

//defines customers cluster to each satellite
void Grasp::assignCustomersToSatellites() {
	//needs implementation
}

//runs grasp procedure
void Grasp::runGrasp() {
	//needs implementation
}

//returns two echelon solution 
TwoEchelonSolution Grasp::getGraspSolution() {
	return graspSolution;
}

