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

const int POPULATION_SIZE = 50;
const int NUMBER_OF_NODES = 15;
const int DEPOT_NODE = 7;
const int CAPACITY_LIMIT = 10;
const int KCHAIN_LENGTH = 5;
const int SWAPCHAIN_LENGTH = 5;


int main()
{
    //random int generator
    std::random_device rd;
    std::mt19937 gen(rd());

    //populate demand vector
    std::map<int, int> demand;
    //demand.insert(std::pair<int, int>(0, 0));
    demand.insert(std::pair<int, int>(1, 2));
    demand.insert(std::pair<int, int>(2, 3));
    //demand.insert(std::pair<int, int>(3, 1));
    demand.insert(std::pair<int, int>(4, 2));
    demand.insert(std::pair<int, int>(5, 3));
    demand.insert(std::pair<int, int>(6, 2));
    demand.insert(std::pair<int, int>(7, 2));
    //demand.insert(std::pair<int, int>(8, 4));
    demand.insert(std::pair<int, int>(9, 2));
    demand.insert(std::pair<int, int>(10, 3));
    demand.insert(std::pair<int, int>(11, 2));
    demand.insert(std::pair<int, int>(12, 3));
    //demand.insert(std::pair<int, int>(13, 2));
    demand.insert(std::pair<int, int>(14, 3));
    demand.insert(std::pair<int, int>(15, 2));
   
    //0,3,8,13
    //populate distance matrix
    std::vector<std::vector<double>> distance;
    for (int i = 0; i <= NUMBER_OF_NODES; ++i) {
        std::vector<double> dist;
        for (int j = 0; j <= NUMBER_OF_NODES; ++j) {
            if (i == j) {
                dist.push_back(0);
            }
            else {
                std::uniform_int_distribution<> distr(1, 9);
                int node = distr(gen);
                dist.push_back(node);
            }
        }
        for (auto it: dist) {
            std::cout << it << " ";
        }
        std::cout << ";" << std::endl;
        distance.push_back(dist);
        dist.clear();
    }
    //customers
    std::set<int> cusCluster;
    cusCluster.insert(1);
    cusCluster.insert(2);
    cusCluster.insert(4);
    cusCluster.insert(5);
    cusCluster.insert(6);
    cusCluster.insert(7);
    cusCluster.insert(9);
    cusCluster.insert(10);
    cusCluster.insert(11);
    cusCluster.insert(12);
    cusCluster.insert(14);
    cusCluster.insert(15);
    //satellites
    std::set<int> satCluster;
    satCluster.insert(0);
    satCluster.insert(3);
    satCluster.insert(8);
    satCluster.insert(13);
    
    std::set<int> customersMustServeByFirstEchelon;


    ProblemParameters probParam(demand, distance, cusCluster, satCluster, customersMustServeByFirstEchelon, 12, 5, 5, 3);
    Initialsolution initSol(probParam);
    initSol.runInitialSolution();
    TwoEchelonSolution initsl = initSol.getTwoEchelonSolution();
    Localsearch localSrch(initSol.getTwoEchelonSolution(), initSol.getTwoEchelonSolution());
    std::cout << "the priority queue is being created" << std::endl;
    localSrch.createCustomersPriorityQueue();
    std::cout << "the local search is being done" << std::endl;
    localSrch.runLocalSearch();
    std::cout << "\nShow the improved solution" << std::endl;
    localSrch.showCurrentSolution();
    //std::cout << "\nShow the current solution" << std::endl;
    //localSrch.showBestSolution();

    Pathrelinking path(localSrch.getCurrentSolution(), localSrch.getCurrentSolution());
    path.runPathRelinking();
   
   
    return 0;
}


