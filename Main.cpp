#include <iostream>

#include "Tabusearch.h"
#include "Geneticalgorithm.h"
#include "Initialsolution.h"
#include "Feasibilitysearch.h"
#include "Localsearch.h"
#include "Pathrelinking.h"

const int NumberOfNodes = 100;
const int DepotNode = 10;
const int CapacityLimit = 20;
const int KChain = 5;
const int SwapChain = 5;


int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(1, 9);
    int node = distr(gen);

    std::list<int> container1;
    for (int i = 0; i < NumberOfNodes; ++i) {
        if (i% DepotNode ==0) {
            container1.push_back(DepotNode);
        }
        else {
            container1.push_back(i);
        }
    }
    container1.push_back(DepotNode);

    std::vector<int> demand;
    std::vector<std::vector<double>> distance;

    //populate demand vector
    for (int i = 0; i < NumberOfNodes; ++i) {
        demand.push_back(2);
    }
    //populate distance matrix
    for (int i = 0; i < NumberOfNodes; ++i) {
        std::vector<double> dist;
        for (int j = 0; j < NumberOfNodes; ++j) {
            if (i == 0 && j == 0) {
                dist.push_back(0);
            }
            else {
                int node = distr(gen);
                dist.push_back(node);
               // std::cout << i << " " << j << " " << node<< std::endl;
            }
        }
        distance.push_back(dist);
        dist.clear();
    }

    double cost = 0;
    int preVal = 0;
    int val = 0;
    int postVal = 0;
    for (auto &it : container1) {
        //std::cout << "The node is : " << it << std::endl;
        if (it != 4) {
            val = it;
        }
        else {
            val = 0;
        }
        cost = cost + distance[preVal][val];
        //std::cout << preVal << " " << val << " " << distance[preVal][val] << " " << std::endl;
        preVal = val;
    }
    cost = cost + distance[preVal][postVal];
    std::cout << "The total cost is " << cost << std::endl;

    
    FeasibleSolution sol(container1, cost, DepotNode, 1, 2);
    sol.showSolution();
    
    Tabusearch tabu(sol, demand, distance, KChain, SwapChain, CapacityLimit);
    tabu.runTabuSearch();

        
    return 0;
}
