#include <iostream>
#include <iterator>
#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <fstream>
#include <string>
#include <sstream>

#include "Tabusearch.h"
#include "Geneticalgorithm.h"
#include "Initialsolution.h"
#include "Feasibilitysearch.h"
#include "Localsearch.h"
#include "Pathrelinking.h"
#include "Grasp.h"
#include "GRASPPR.h"
#include "TspH.h"

//const int POPULATION_SIZE = 50;
const int NUMBER_OF_NODES = 60;//cus+sats
const int DEPOT_NODE = 0;
//const int CAPACITY_LIMIT = 10;
const int KCHAIN_LENGTH = 5;
const int SWAPCHAIN_LENGTH = 5;
const int TRUCK_CAP_LIMIT = 20;
const int UAV_CAP_LIMIT = 3;
const int MAX_NO_TRUCK = 5;
const int MAX_NO_UAV = 3;

int main(){
    
    std::fstream myFile;
    //open raod distance data file
    myFile.open("DistanceMatrix.txt", std::ios::in);//read
    std::vector<std::vector<double>> roadDistance;
    if (myFile.is_open()) {
        std::string line, word;
        std::istringstream iss;
        int rowNum = 0;
        while (std::getline(myFile, line)) {
            if (rowNum > 0) {
                std::vector<double> dist;
                iss.clear();
                iss.str(line);
                int colNum = 0;
                while (iss.good()) {
                    iss >> word;
                    //&& colNum <= NUMBER_OF_NODES
                    if (colNum > 0) {
                        double abc = std::stod(word);//abc == road distance
                        //std::cout << abc << std::endl;
                        dist.push_back(abc);
                    }
                    colNum += 1;
                }
                roadDistance.push_back(dist);
            }
            rowNum += 1;
        }
        myFile.close();
    }
    
    //std::cout << "Show road distance matrix" << std::endl;
    //std::cout << "size of the roadDistance matrix : " << roadDistance.size() << std::endl;
    //for (int i = 0; i < roadDistance.size(); i++) {
    //    for (int j = 0; j < roadDistance.size(); j++) {
    //        std::cout << roadDistance[i][j] << " ";
    //    }
    //    std::cout << ";" << std::endl;
    //}
    

    //open aerial distance data file
    myFile.open("GCDistanceMatrix.txt", std::ios::in);//read
    std::vector<std::vector<double>> aerialDistance;
    if (myFile.is_open()) {
        std::string line, word;
        std::istringstream iss;
        int rowNum = 0;
        while (std::getline(myFile, line)) {
            if (rowNum > 0) {
                std::vector<double> dist;
                iss.clear();
                iss.str(line);
                int colNum = 0;
                while (iss.good()) {
                    iss >> word;
                    if (colNum > 0) {
                        double abc = std::stod(word);//abc == road distance
                        //std::cout << abc << std::endl;
                        dist.push_back(abc);
                    }
                    colNum += 1;
                }
                aerialDistance.push_back(dist);
            }
            rowNum += 1;
        }
        myFile.close();
    }
    //
    //std::cout << "Show aerial distance matrix" << std::endl;
    //std::cout << "size of the aerialDistance matrix : " << aerialDistance.size() << std::endl;
    //for (int i = 0; i < aerialDistance.size(); i++) {
    //    for (int j = 0; j < aerialDistance.size(); j++) {
    //        std::cout << aerialDistance[i][j] << " ";
    //    }
    //    std::cout << ";" << std::endl;
    //}
    

    std::cout << " Populate customer cluster" << std::endl;
    std::set<int> cusCluster;
    for (int i = 11; i <= NUMBER_OF_NODES; i++) {
        cusCluster.insert(i);
    }
    std::cout << "Populate depot/satellite cluster" << std::endl;
    std::set<int> satCluster;
    for (int i = 0; i <= 10; i++) {
        satCluster.insert(i);
    }
    std::cout << "Populate demand for the customers" << std::endl;
    std::map<int, int> demand;
    for (int i = 11; i <= NUMBER_OF_NODES; i++) {
        demand.insert(std::pair<int, int>(i,1));//unit demand
    }

    std::cout << " Populate customers must serve by first echelon vehicle" << std::endl;
    std::set<int> customersMustServeByFirstEchelon;
    customersMustServeByFirstEchelon.insert(15);
    customersMustServeByFirstEchelon.insert(20);

    //std::cout << "\nRoad distance 60 to 4 is : " << roadDistance[4][60] << std::endl;
    
    ProblemParameters probParam(demand, aerialDistance, cusCluster, satCluster, customersMustServeByFirstEchelon, TRUCK_CAP_LIMIT, UAV_CAP_LIMIT, MAX_NO_TRUCK,  MAX_NO_UAV);
    std::cout << "\nRun the main GRASP with Path Relinking algorithm\n" << std::endl;
    GRASPPR Alg(probParam);
    Alg.runGraspPr();
    std::cout << "\nThe GRASP with Path Relinking solution : " << std::endl;
    Alg.getGraspPRSolution().showTwoEchelonSolution();
    


    /*
    std::vector<std::vector<double>> costMatrix;
    std::vector<double> vec{10,20,30,12,15};
    std::vector<double> vec1{5,21,10,10,25};
    std::vector<double> vec2{20,10,40,10,5};
    std::vector<double> vec3{5,11,20,8,5};
    std::vector<double> vec4{15,15,25,18,10};
    costMatrix.push_back(vec);
    costMatrix.push_back(vec1);
    costMatrix.push_back(vec2);
    costMatrix.push_back(vec3);
    costMatrix.push_back(vec4);
    */
    /*
    for (int i = 0; i < 5; i++) {
        std::vector<double> vec;
        for (int j = 0; j < 5; j++) {
            int val = 0;
            if (i == j) {
                val = 100000;
            }
            else {
                val = (i + j) * 4;
            }
            vec.push_back(val);
        }
        costMatrix.push_back(vec);
        vec.clear();
    }
    */
    /*
    std::cout << "\nShow the cost matrix." << std::endl;
    for (int i = 0; i < costMatrix.size(); i++) {
        for (int j = 0; j < costMatrix.size(); j++) {
            std::cout << costMatrix[i][j] << " ";
        }
        std::cout << ";" << std::endl;
    }
    
    std::vector<int> tsp;
    tsp.push_back(1);
    tsp.push_back(2);
    tsp.push_back(4);
    tsp.push_back(3);
    tsp.push_back(0);
    
    
    TspH assignment = TspH(tsp, costMatrix);
    assignment.solveAssignmentProblem();
    assignment.showAssignmentOutcomes();
    */

    return 0;
}



