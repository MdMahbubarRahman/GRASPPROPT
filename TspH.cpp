#include "TspH.h"

//default constructor
TspH::TspH() {
	depotNode = INFINITY;
	reducedCost = 0.0;
	baseValue = 0.0;
	initSolCost = 0.0;
	finalSolCost = 0.0;
	tspOptimal = false;
	assignmentOptimal = false;
	reducedCostMatrixSize = 0;
}

//constructor
TspH::TspH(std::vector<int> initTsp, std::vector<std::vector<double>> cMatrix) {
	depotNode = 3;
	reducedCost = 0.0;
	baseValue = 0.0;
	initSolCost = 0.0;
	finalSolCost = 0.0;
	tspOptimal = false;
	assignmentOptimal = false;
	initialTsp = initTsp;
	costMatrix = cMatrix;
	//populate reduced const matrix
	reducedCostMatrixSize = initialTsp.size();
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		std::vector<double> costVec;
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			double val = 0;
			i == j? val = INFINITY : val = costMatrix[initialTsp.at(i)][initialTsp.at(j)];
			costVec.push_back(val);
		}
		reducedCostMatrix.push_back(costVec);
	}
	//initial tsp cost
	for (int i = 0; i < reducedCostMatrixSize-1; i++) {
		initSolCost += costMatrix[initialTsp.at(i)][initialTsp.at(i+1)];
	}
	initSolCost += costMatrix[initialTsp.at(reducedCostMatrixSize - 1)][initialTsp.at(0)];
}


//performs row reduction
void TspH::performRowColumnReduction() {
	std::cout << "\nShow reduced cost matrix ." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
	//perform row reduction
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		double minVal = INFINITY;
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			if (minVal > reducedCostMatrix[i][j]) {
				minVal = reducedCostMatrix[i][j];
			}
		}
		//subtract min val from each the corresponding row
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			reducedCostMatrix[i][j] = reducedCostMatrix[i][j]-minVal;
		}
		reducedCost += minVal;
	}
	std::cout << "\nShow reduced cost matrix after row reduction." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}


	//perform column reduction
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		double minVal = INFINITY;
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			if (minVal > reducedCostMatrix[j][i]) {
				minVal = reducedCostMatrix[j][i];
			}
		}
		//subtract min val from each the corresponding column
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			reducedCostMatrix[j][i] = reducedCostMatrix[j][i] - minVal;
		}
		reducedCost += minVal;
	}
	reducedCost == initSolCost ? tspOptimal = true : tspOptimal = false;
	std::cout << "\nShow reduced cost matrix after row column reduction." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
}


//cover zero values and check for optimal assignment problem solution
void TspH::coverMinimumValueAssignmentsAndCheckOptimality(bool &assignmentOptimal, std::set<int> &coveredRows, std::set<int> &coveredColumns, std::map<int, int> &boxPoints, std::list<IntersectionPoint> &twiceCoveredPoints) {
	//row scaning
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		int numOfZeros = 0;
		int colPos = 0;
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			auto iter = coveredColumns.find(j);
			if (iter == coveredColumns.end() && reducedCostMatrix[i][j] == 0) {
				numOfZeros += 1;
				colPos = j;
			}
		}
		if (numOfZeros == 1) {
			coveredColumns.insert(colPos);
			boxPoints.insert(std::pair<int, int>(i, colPos));
		}
	}
	//column scaning
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredColumns.find(i);
		if (it == coveredColumns.end()) {
			int numOfZeros = 0;
			int rowPos = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto iter = coveredRows.find(j);
				if (iter == coveredRows.end() && reducedCostMatrix[j][i] == 0) {
					numOfZeros += 1;
					rowPos = j;
				}
			}
			if (numOfZeros == 1) {
				coveredRows.insert(rowPos);
				boxPoints.insert(std::pair<int, int>(rowPos, i));
			}
		}
	}
	//check if there is any uncovered zero
	std::map<int, int> numZeroInRow;
	std::map<int, int> numZeroInCol;
	std::cout << "\nUpdate number of uncovered zeros in each row." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto itr = coveredRows.find(i);
		if (itr == coveredRows.end()) {
			int numZero = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto it = coveredColumns.find(j);
				if (it == coveredColumns.end() && reducedCostMatrix[i][j] == 0){
					numZero += 1;
				}
			}
			numZeroInRow.insert(std::pair<int, int>(i, numZero));
		}
	}
	std::cout << "\nUpdate number of uncovered zeros in each column." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto itr = coveredColumns.find(i);
		if (itr == coveredColumns.end()) {
			int numZero = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto it = coveredRows.find(j);
				if (it == coveredRows.end() && reducedCostMatrix[j][i] == 0) {
					numZero += 1;
				}
			}
			numZeroInCol.insert(std::pair<int, int>(i, numZero));
		}
	}
	//Number of uncovered zeros
	int numOfZeros = 0;
	for (auto &it: numZeroInCol) {
		numOfZeros += it.second;
	}
	if (numOfZeros == 0) {
		std::cout << "\nAll zeros are covered." << std::endl;
	}
	else {
		//cover uncovered zeros using diagonal procedure
		std::set<int> newlyCoveredCols;
		while (numOfZeros != 0) {
			std::set<int> newlyCoveredRows;
			int rowID = INFINITY;
			int init = 0;
			for (auto& it : numZeroInCol) {//accending through the column numbers
				auto col = newlyCoveredCols.find(it.first);
				if (col == newlyCoveredCols.end()) {
					if (init == 0) {
						if (it.second > 0) {
							bool flag = false;
							for (auto& itt : numZeroInRow) {
								if (flag == false) {
									if (reducedCostMatrix[itt.first][it.first] == 0) {
										coveredRows.insert(itt.first);
										newlyCoveredRows.insert(itt.first);
										newlyCoveredCols.insert(it.first);
										boxPoints.insert(std::pair<int, int>(itt.first, it.first));
										flag = true;
									}
								}
								else {
									rowID = itt.first;
									break;
								}
							}
						}
					}
					else {
						if (rowID != INFINITY) {
							if (reducedCostMatrix[rowID][it.first] == 0) {
								coveredRows.insert(rowID);
								newlyCoveredRows.insert(rowID);
								newlyCoveredCols.insert(it.first);
								boxPoints.insert(std::pair<int, int>(rowID, it.first));
								//find the next row id
								bool flg = false;
								for (auto& iter : numZeroInRow) {
									if (flg == false) {
										if (iter.first == rowID) {
											flg = true;
										}
									}
									else {
										rowID = iter.first;
										break;
									}
								}
							}
						}
					}
					init += 1;
				}
			}
			//delete newly covered rows
			for (auto it: newlyCoveredRows) {
				numZeroInRow.erase(it);
			}
			//again calculate number of uncovered zeros
			numOfZeros = 0;
			for (auto& it : numZeroInRow) {
				numOfZeros += it.second;
			}
		}
	}
	std::cout << "\nShow the intersection points" << std::endl;
	for (auto i : coveredRows) {
		for (auto j : coveredColumns) {
			IntersectionPoint point = IntersectionPoint(i,j);
			std::cout << "Row : " << i << " Column : " << j << std::endl;
			twiceCoveredPoints.push_back(point);
		}
	}
	//check assignment optimality
	std::cout << "\nNumber of rows covered : " << coveredRows.size() << ", Number of columns covered : " << coveredColumns.size() << std::endl;
	if ((coveredColumns.size()+coveredRows.size()) == reducedCostMatrixSize) {
		assignmentOptimal = true;
		std::cout << "\nAssignment optimal : " << assignmentOptimal << std::endl;
	}
}

//performs row column reduction for uncovered nonzero cells
void TspH::performRowColumnReductionForUncoveredCells(bool& assignmentOptimal, std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints, std::list<IntersectionPoint>& twiceCoveredPoints) {
	//perform row reduction
	std::set<int> uncoveredRows, uncoveredColumns;
	//update uncovered rows
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredRows.find(i);
		if (it == coveredRows.end()) {
			uncoveredRows.insert(i);
		}
	}
	//update uncovered columns
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredColumns.find(i);
		if (it == coveredColumns.end()) {
			uncoveredColumns.insert(i);
		}
	}
	//find minimum value from the uncovered cells
	double minVal = INFINITY;
	for (auto it: uncoveredRows) {
		for (auto itt: uncoveredColumns) {
			if (minVal > reducedCostMatrix[it][itt]) {
				minVal = reducedCostMatrix[it][itt];
			}
		}
	}
	//subtract min value from all uncovered cells
	for (auto it : uncoveredRows) {
		for (auto itt : uncoveredColumns) {
			reducedCostMatrix[it][itt] = reducedCostMatrix[it][itt] - minVal;
		}
		reducedCost += minVal;
	}
	//add min value to all intersection points
	for (auto &point : twiceCoveredPoints) {
		reducedCostMatrix[point.row][point.column] = reducedCostMatrix[point.row][point.column] + minVal;
	}
	std::cout << "\nShow reduced cost matrix after updating uncovered cells." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
}

//generates routes from assignment solution
void TspH::populateListOfRoutes(std::map<int, int> boxPoints) {
	std::cout << "\nPrint the contents of the newAssignment" << std::endl;
	for (auto& it : boxPoints) {
		std::wcout << it.first << " " << it.second << std::endl;
	}
	//start route from depot
	std::vector<int> route;
	int sourceNode = 0;
	int destinationNode = 0;
	int routeStartNode = 0;
	std::map<int, int> newAssignments = boxPoints;
	int i = 0;
	while (!newAssignments.empty()) {
		if (i == 0) {
			auto it = newAssignments.begin();
			routeStartNode = (*it).first;
			destinationNode = (*it).second;
			route.push_back(destinationNode);
			newAssignments.erase(routeStartNode);
			i += 1;
		}
		else {
			sourceNode = destinationNode;
			destinationNode = newAssignments[sourceNode];
			route.push_back(destinationNode);
			newAssignments.erase(sourceNode);
			if (destinationNode == routeStartNode) {
				listOfRoutes.push_back(route);
				route.clear();
				i = 0;
			}
			else {
				i += 1;
			}
		}
	}
}

//finds the number of routes generated from assignment solution
void TspH::checkListOfRoutes(std::map<int, int> boxPoints, std::list<std::vector<int>> &routeList) {
	//start route from depot
	routeList.clear();
	std::vector<int> route;
	int sourceNode = 0;
	int destinationNode = 0;
	int routeStartNode = 0;
	std::map<int, int> newAssignments = boxPoints;
	int i = 0;
	std::cout << "\nStart routes construction." << std::endl;
	while (!newAssignments.empty()) {
		if (i == 0) {
			auto it = newAssignments.begin();
			routeStartNode = (*it).first;
			destinationNode = (*it).second;
			route.push_back(destinationNode);
			newAssignments.erase(routeStartNode);
			i += 1;
		}
		else {
			sourceNode = destinationNode;
			destinationNode = newAssignments[sourceNode];
			route.push_back(destinationNode);
			newAssignments.erase(sourceNode);
			if (destinationNode == routeStartNode) {
				routeList.push_back(route);
				route.clear();
				i = 0;
			}
			else {
				i += 1;
			}
		}
	}
	std::cout << "\nEnd routes construction." << std::endl;
}

//generates a tsp solution which is a valid upper bound
void TspH::generateUpperBoundTsp(std::list<std::vector<int>> listOfRoutes) {
	//find the values in the reduced cost matrix and sort them in assending order in a list
	int iteration = 0;
	std::set<double> values;
	std::list<std::vector<int>> routeList = listOfRoutes;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			if (reducedCostMatrix[i][j] != INFINITY) {
				values.insert(reducedCostMatrix[i][j]);
			}
		}
	}
	auto it = values.begin();
	values.erase(*it);
	while (routeList.size() != 1) {
		auto val = values.begin();
		for (int i = 0; i < reducedCostMatrixSize; i++) {
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				if (reducedCostMatrix[i][j] == (*val)) {
					reducedCostMatrix[i][j] = 0.0;
				}
			}
		}
		values.erase(*val);
		//get assignment solution
		std::set<int> coveredRows;
		std::set<int> coveredColumns;
		std::list<IntersectionPoint> twiceCoveredPoints;
		boxPoints.clear();
		coverMinimumValueAssignmentsAndCheckOptimality(assignmentOptimal, coveredRows, coveredColumns, boxPoints, twiceCoveredPoints);
		checkListOfRoutes(boxPoints, routeList);
		iteration += 1;
	}
	//print the upper bound tsp route
	std::cout << "\nNumber of upperbound iterations : " << iteration << std::endl;
	std::cout << "\nTSP upperbound sol: " << std::endl;
	for (auto &it: routeList) {
		for (auto itt: it) {
			std::cout << itt << " ";
		}
		std::cout << ";" << std::endl;
	}
}

//solves relaxed tsp problem which is an assignment problem and generates lower bound of TSP
void TspH::solveAssignmentProblem() {
	std::set<int> coveredRows;
	std::set<int> coveredColumns;
	std::list<IntersectionPoint> twiceCoveredPoints;
	performRowColumnReduction();
	if (tspOptimal == true) {
		std::cout << "\nThe initial tsp solution is optimal! No need to look for optimal tsp solution." << std::endl;
	}
	else {
		while (assignmentOptimal != true) {
			coverMinimumValueAssignmentsAndCheckOptimality(assignmentOptimal, coveredRows, coveredColumns, boxPoints, twiceCoveredPoints);
			if (assignmentOptimal != true) {
				performRowColumnReductionForUncoveredCells(assignmentOptimal, coveredRows, coveredColumns, boxPoints, twiceCoveredPoints);
				coveredColumns.clear();
				coveredRows.clear();
				boxPoints.clear();
				twiceCoveredPoints.clear();
			}
			else {
				std::cout << "\nOptimum assignment is found!" << std::endl;
				populateListOfRoutes(boxPoints);
				if (listOfRoutes.size()>1) {
					generateUpperBoundTsp(listOfRoutes);
				}
			}
		}
	}
}

//two-opt heuristic algorithm for tsp route
void TspH::twoOptAlgorithm() {
	std::vector<int> initialTspSol;//representation is important==>>depoNode/node/node.../node/depoNode
	std::vector<int> finalTspSol;
	std::vector<std::vector<double>> costMatrix;
	double costA = 0;
	double costB = 0;
	double costC = 0;
	double costD = 0;
	double costDiff = 0;
	int c1 = 0;
	int c2 = 0;
	int c3 = 0;
	int c4 = 0;
	if (initialTspSol.size() >= 5) {
		for (int i = 0; i < initialTspSol.size()-1; i++) {
			for (int k = 0; k < i-1; k++) {
				costA = costMatrix[initialTspSol.at(i)][initialTspSol.at(i+1)];
				costB = costMatrix[initialTspSol.at(k)][initialTspSol.at(k+1)];
				costC = costMatrix[initialTspSol.at(i)][initialTspSol.at(k)];
				costD = costMatrix[initialTspSol.at(i+1)][initialTspSol.at(k + 1)];
				if (costDiff < (costA + costB - costC - costD)) {
					costDiff = (costA + costB - costC - costD);
					c1 = i;
					c2 = i+1;
					c3 = k;
					c4 = k+1;
				}
			}
			if (i == 0) {
				for (int j = i + 1; j < initialTspSol.size() - 2; j++) {
					costA = costMatrix[initialTspSol.at(i)][initialTspSol.at(i + 1)];
					costB = costMatrix[initialTspSol.at(j)][initialTspSol.at(j + 1)];
					costC = costMatrix[initialTspSol.at(i)][initialTspSol.at(j)];
					costD = costMatrix[initialTspSol.at(i + 1)][initialTspSol.at(j + 1)];
					if (costDiff < (costA + costB - costC - costD)) {
						costDiff = (costA + costB - costC - costD);
						c1 = i;
						c2 = i + 1;
						c3 = j;
						c4 = j + 1;
					}
				}			
			}
			else {
				for (int j = i + 1; j < initialTspSol.size() - 1; j++) {
					costA = costMatrix[initialTspSol.at(i)][initialTspSol.at(i + 1)];
					costB = costMatrix[initialTspSol.at(j)][initialTspSol.at(j + 1)];
					costC = costMatrix[initialTspSol.at(i)][initialTspSol.at(j)];
					costD = costMatrix[initialTspSol.at(i + 1)][initialTspSol.at(j + 1)];
					if (costDiff < (costA + costB - costC - costD)) {
						costDiff = (costA + costB - costC - costD);
						c1 = i;
						c2 = i + 1;
						c3 = j;
						c4 = j + 1;
					}
				}
			}
		}
		std::cout << "\nThe improved route is the following." << std::endl;
		int val1 = initialTspSol.at(c2);
		int val2 = initialTspSol.at(c3);
		initialTspSol[c2] = val2;
		initialTspSol[c3] = val1;
		std::cout << "\nThe improved route is :" << std::endl;
		for (auto it: initialTsp) {
			std::cout << it << " ";
		}
		std::cout << ";" << std::endl;
	}
	else {
		std::cout << "\nNumber of elements in this route is less than or equal to 3. And it is an optimal tsp route." << std::endl;
	}
}

//heuristic method to create tsp route from assinment solution
void TspH::routeConstructionHeuristic() {
	std::vector<int> initialTspSol;
	std::vector<int> finalTspSol;
	std::vector<std::vector<double>> costMatrix;
	if (initialTspSol.size() >= 4) {
		for (int i = 0; i < initialTspSol.size(); i++) {
			for (int j = 0; j < initialTspSol.size(); j++) {

			}
		}
	}
	else {

	}
	//needs work
}

//shows results
void TspH::showAssignmentOutcomes() {
	std::cout << "\nShow the reduced cost matrix : " << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
	std::cout << "\nShow assignments : " << std::endl;
	for (auto& it : boxPoints) {
		std::cout << it.first << " --> " << it.second << std::endl;
	}
	std::cout << "\nShow the routes : " << std::endl;
	for (auto& it : listOfRoutes) {
		for (auto itt : it) {
			std::cout << itt << " ";
		}
		std::cout << ";" << std::endl;
	}
}


