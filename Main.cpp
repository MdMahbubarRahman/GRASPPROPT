#include <iostream>
#include <iterator>
#include <cmath>

#include "Tabusearch.h"
#include "Geneticalgorithm.h"
#include "Initialsolution.h"
#include "Feasibilitysearch.h"
#include "Localsearch.h"
#include "Pathrelinking.h"

const int POPULATION_SIZE = 50;
const int NUMBER_OF_NODES = 20;
const int DEPOT_NODE = 20;
const int CAPACITY_LIMIT = 10;
const int KCHAIN_LENGTH = 5;
const int SWAPCHAIN_LENGTH = 5;


int main()
{
    //random int generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    //populate demand vector
    std::vector<int> demand;
    for (int i = 0; i < NUMBER_OF_NODES; ++i) {
        if (i == 0) {
            demand.push_back(0);
        }
        else {
            demand.push_back(2);
        }
    }

    //populate distance matrix
    std::vector<std::vector<double>> distance;
    for (int i = 0; i < NUMBER_OF_NODES; ++i) {
        std::vector<double> dist;
        for (int j = 0; j < NUMBER_OF_NODES; ++j) {
            if (i == 0 && j == 0) {
                dist.push_back(0);
            }
            else {
                std::uniform_int_distribution<> distr(1, 9);
                int node = distr(gen);
                dist.push_back(node);
            }
        }
        distance.push_back(dist);
        dist.clear();
    }
    
         

    Geneticalgorithm ga(POPULATION_SIZE, NUMBER_OF_NODES, DEPOT_NODE, CAPACITY_LIMIT, KCHAIN_LENGTH, SWAPCHAIN_LENGTH, demand, distance);
    ga.runGeneticAlgorithm();
      

    return 0;
}
