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
    
    ProblemParameters probParam(demand, roadDistance, aerialDistance, cusCluster, satCluster, customersMustServeByFirstEchelon, TRUCK_CAP_LIMIT, UAV_CAP_LIMIT, MAX_NO_TRUCK,  MAX_NO_UAV);
    //Initialsolution init(probParam);
    //init.runInitialSolution();
    
    std::cout << "\nRun the main GRASP with Path Relinking algorithm\n" << std::endl;
    GRASPPR Alg(probParam);
    Alg.runGraspPr();
    /*
    std::cout << "\nThe GRASP with Path Relinking solution : " << std::endl;
    Alg.getGraspPRSolution().showTwoEchelonSolution();
    */
    return 0;
}



