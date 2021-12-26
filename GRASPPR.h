#include <iostream>
#include <iterator>
#include <cmath>
#include <vector>
#include <list>
#include <map>

#include "Tabusearch.h"
#include "Geneticalgorithm.h"
#include "Initialsolution.h"
#include "Feasibilitysearch.h"
#include "Localsearch.h"
#include "Pathrelinking.h"
#include "Grasp.h"

#ifndef GRASPPR_H
#define GRASPPR_H 

class GRASPPR{
private:
	ProblemParameters probParam;
	TwoEchelonSolution initialTwoEchelonSol;
	TwoEchelonSolution currentBestTwoEchelonSol;
	TwoEchelonSolution finalTwoEchelonSol;
public:
	GRASPPR();
	GRASPPR(ProblemParameters probParam, TwoEchelonSolution initialSol, TwoEchelonSolution incumbentSol, TwoEchelonSolution finalSol);
	GRASPPR(ProblemParameters probParam);
	TwoEchelonSolution getInitialSolution();
	TwoEchelonSolution getGraspPRSolution();
	void runGraspPr();
};



#endif //GRASPPR_H
