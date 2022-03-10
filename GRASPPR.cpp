#include "GRASPPR.h"

//default constructor is called
GRASPPR::GRASPPR() {
	std::cout << "Default GRASPPR constructor is called!\n" << std::endl;
}

//copy constructor is called
GRASPPR::GRASPPR(ProblemParameters pbParam, TwoEchelonSolution initSol, TwoEchelonSolution incumbentSol, TwoEchelonSolution finalSol) {
	probParam = pbParam;
	initialTwoEchelonSol = initSol;
	currentBestTwoEchelonSol = incumbentSol;
	finalTwoEchelonSol = finalSol;
}

//start with problem parameters
GRASPPR::GRASPPR(ProblemParameters pbParam) {
	probParam = pbParam;
}

//run grasp with path relinking algorithm
void GRASPPR::runGraspPr() {
	auto start = high_resolution_clock::now();
	using std::chrono::duration;
	//generate initial solution
	std::cout << "\nInitial solution process has been started" << std::endl;
	Initialsolution initSol(probParam);
	initSol.runInitialSolution();
	initialTwoEchelonSol = initSol.getTwoEchelonSolution();
	currentBestTwoEchelonSol = initialTwoEchelonSol;
	//completion of initial solution generation
	//initiate adaptive customer to satellite assignment probabililty distribution class
	AdaptiveCustomersAssignmentProbabilityDistribution abc = AdaptiveCustomersAssignmentProbabilityDistribution();
	abc.populateAdaptiveCusAssgnProbDist(probParam);
	//completion of adaptive customer to satellite assignment
	//start the multistart procedure for GRASP with Path Relinking
	int maxIter = 10;
	int iterNumber = 0;
	int localSrchFlag = 0;
	TwoEchelonSolution intermediateSolution;
	std::list<TwoEchelonSolution> graspSols;
	while (iterNumber < maxIter){
		iterNumber += 1;
		std::cout << "Get GRASP solution" << std::endl;
		Grasp grasp = Grasp(abc, probParam);
		grasp.runGrasp();
		intermediateSolution = grasp.getGraspSolution();
		//std::cout << "\nThe first GRASP solution is : \n" << std::endl;
		//intermediateSolution.showTwoEchelonSolution();
		std::cout << "Run FEASIBILITY search on the grasp solution" << std::endl;
		Feasibilitysearch febSrch(probParam, intermediateSolution);
		febSrch.runFeasibilitySearch();
		if (febSrch.getFeasibilityStatus() == 1) {
			intermediateSolution = febSrch.getFeasibilitySearchSolution();
			if (intermediateSolution.getSolutionFitness() <= currentBestTwoEchelonSol.getSolutionFitness()+5000) {  //NOW strictly looks for better solution
				std::cout << "Perform LOCAL search on the current solution" << std::endl;
				Localsearch localSrch(intermediateSolution, probParam);
				localSrch.createCustomersPriorityQueue();
				localSrch.runLocalSearch();
				if (localSrch.isCurrentSolutionImproved() == 1) {
					intermediateSolution = localSrch.getCurrentSolution();
					localSrchFlag += 1;
					std::cout << "Run Path Relinking with the current solution and currentBestSolution" << std::endl;
					Pathrelinking path(intermediateSolution, currentBestTwoEchelonSol);
					path.runPathRelinking();
					intermediateSolution = path.getCurrentBestSolution();
				}
			}
			if (intermediateSolution.getSolutionFitness() < currentBestTwoEchelonSol.getSolutionFitness()) {
				currentBestTwoEchelonSol = intermediateSolution;
			}
		}
		//store current grasp solution
		graspSols.push_back(intermediateSolution);
	}
	//finalTwoEchelonSol = currentBestTwoEchelonSol;
	auto stop = high_resolution_clock::now();
	duration<double, std::milli> ms_double = stop - start;
	//showTabuSolution();
	double secDuration = double(ms_double.count()) / 1000;
	std::cout<<"GRASPPR Duration: " << secDuration << " seconds" << std::endl;
	std::cout << "\nThe initial solution for the Two Echelon VRP is : " << std::endl;
	initialTwoEchelonSol.showTwoEchelonSolution();
	std::cout << "\nShow grasp solutions" << std::endl;
	int i = 0;
	for (auto &it: graspSols) {
		std::cout << "\nGRASP solution number : " << i << std::endl;
		it.showTwoEchelonSolution();
		i++;
	}
	std::cout << "\nNumber of local search is : " << localSrchFlag << std::endl;
	std::cout << "\nThe best solution is : " << std::endl;
	currentBestTwoEchelonSol.showTwoEchelonSolution();
}

//return the solution of the Algorithm
TwoEchelonSolution GRASPPR::getGraspPRSolution() {
	return finalTwoEchelonSol;
}

//return initial two echelon solution of the algorithm
TwoEchelonSolution GRASPPR::getInitialSolution() {
	return initialTwoEchelonSol;
}

