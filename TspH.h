#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <Queue>
#include <map>

#ifndef TSPH_H
#define TSPH_H

class TspH{
private:
	int depotNode;
	double reducedCost;
	double baseValue;
	double initSolCost;
	double finalSolCost;
	bool tspOptimal;
	int reducedCostMatrixSize;
	struct IntersectionPoint {
		int row;
		int column;
		IntersectionPoint(int a, int b) {
			row = a;
			column = b;
		}
	};
	std::map<int, int> boxPoints;
	bool assignmentOptimal;
	std::list<std::vector<int>> listOfRoutes;
	std::vector<int> initialTsp;//tsp presentation node/node/.../node/depotnode
	std::vector<int> incumbentTsp;
	std::vector<int> optimalTsp;
	std::vector<std::vector<double>> costMatrix;
	std::vector<std::vector<double>> reducedCostMatrix;
public:
	TspH();
	TspH(std::vector<int> initialTsp, std::vector<std::vector<double>> costMatrix);
	void solveAssignmentProblem();
	void performRowColumnReduction();
	void coverMinimumValueAssignmentsAndCheckOptimality(bool &assignmentOptimal, std::set<int> &coveredRows, std::set<int> &coveredColumns, std::map<int, int> &boxPoints, std::list<IntersectionPoint> &twiceCoveredPoints);
	void performRowColumnReductionForUncoveredCells(bool& assignmentOptimal, std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints, std::list<IntersectionPoint>& twiceCoveredPoints);
	void generateUpperBoundTsp(std::list<std::vector<int>> listOfRoutes);
	void showAssignmentOutcomes();
	void populateListOfRoutes(std::map<int, int> boxPoints);
	void checkListOfRoutes(std::map<int, int> boxPoints, std::list<std::vector<int>> &routeList);
	void twoOptAlgorithm();
	void routeConstructionHeuristic();
};



#endif // !TSPH_H


