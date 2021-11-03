#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <iterator>
#include <cmath>
#include <set>

#include "Tabusearch.h"

#ifndef GENETICALGOIRTHM_H
#define GENETICALGOIRTHM_H

/*
TODO:Based on "EVE-OPT: a hybrid algorithm for the capacitated vehicle routing problem"
paper implement following functions and algorithm:
	1. Elitism function
	2. Replication function
	3. Chain Mutation algorithm
	4. Partial Crossover algorihtm
	5. Generate initial population based on "Scheduling of vehicles from a central depot
	to a number of delivery points by G Clarke, JW Wright " paper.
	6. Genetic algorithm incorporating Tabu search algorithm

More knowledge to develop above functions and algorithms can be found from following papers:
	1. [BOOK] An introduction to genetic algorithms by M Mitchell.
	2. A simple and effective evolutionary algorithm for the vehicle routing problem by C Prins 
	3. Active guided evolution strategies for large-scale vehicle routing problems with time windows by D Mester, O Bräysy 
*/

/*
The initial populatioon is made up of an individual obtained from the Clarke and Wright heuristic (1964) and a number of individuals with 
chromosomes randomly generated. During the whole evolutionary process the population size remains constant.
*/

//chromosome
class Chromosome {
private:
	std::vector<int> chromosome;
	double fitness;
	int sepInt;
	bool isFeasible;
	int maxRouteCapacity;
public:
	Chromosome();
	Chromosome(const Chromosome & chrom);
	Chromosome(std::vector<int> chromosome,	double fitness,	int sepInt,	bool isFeasible, int maxRouteCapacity);
	std::vector<int> getChromosomeRepresentation();
	double getFitness();
	int getSize();
	int getSepInt();
	bool getFeasibilityStatus();
	int getMaxRouteCapacity();
	void showChromosome();
	void updateToFeasibleChromosome(std::map<int, int> demand, std::vector<std::vector<double>> distance);
};


//cross over operator
class CrossOver {
private:
	Chromosome parent1;
	Chromosome parent2;
	Chromosome offspring;
public:
	CrossOver();
	CrossOver(const CrossOver & crossr);
	CrossOver(Chromosome parent1, Chromosome parent2);
	void showParent1();
	void showParent2();
	void showOffspring();
	Chromosome getParent1();
	Chromosome getParent2();
	Chromosome getOffspring();
	void performPertiallyMappedCrossover(std::map<int, int> demand, std::vector<std::vector<double>> distance);
};


//mutation operator
class Mutation {
	Chromosome originalChrom;
	Chromosome mutatedChrom;
public:
	Mutation();
	Mutation(const Mutation & mutatn);
	Mutation(Chromosome originalChrom);
	Chromosome getOriginalChromosome();
	Chromosome getMutatedOffspring();
	void showOriginalChromosome();
	void showMutatedChromosome();
	void performChainMutation(std::map<int, int> demand, std::vector<std::vector<double>> distance);
}; 


//comparator for reproduction probability
class Comparator {
public:
	bool operator()(Chromosome &a, Chromosome &b);
};

//population of chromosomes
class Population {
	int diversitySize;
	int populationSize;
	int numberOfNodes;
	int depotNode;
	int capacityLimit;
	int kChainLength;
	int sWapLength;
	std::map<int, int> demand;
	std::vector<std::vector<double>> distance;
	//std::priority_queue<Chromosome, std::vector<Chromosome>, Comparator> population;
	std::list<Chromosome> population;
	std::vector<int> customerCluster;
	Chromosome offspring;
	Chromosome populationBest;
	Chromosome crossOverChild;
	Chromosome mutationChild;
public:
	Population();
	Population(const Population & poplatn);
	Population(int populationSize,	int numberOfNodes,	int depotNode,	int capacityLimit,	int kChainLength,	int sWapLength, std::map<int, int> demand, std::vector<std::vector<double>> distance, std::vector<int> customerCluster);
	//std::priority_queue<Chromosome, std::vector<Chromosome>, Comparator> getPopulation();
	std::list<Chromosome> getPopulation();
	Chromosome getOffspring();
	Chromosome getCrossOverChild();
	Chromosome getMutationChild();
	Chromosome getBestChromosome();
	Chromosome generateRandomChromosome();
	void populatePopulation();
	void showPopulation();
	void updatePopulationWithOffspring();
	void manageClones();
	void performCrossOver();
	void performMutation();
	void measureDiversitySize();
	int getDiversitySize();
	FeasibleSolution getFeasibleSolutionFromChromosome(Chromosome chrom);
	Chromosome getChromosomeFromFeasibleSolution(FeasibleSolution febSol);
};

//tournament

//genetic algorithm
class Geneticalgorithm{
private:
	int maxIterations;
	int populationSize;
	int numberOfNodes;
	int depotNode;
	int capacityLimit;
	int kChainLength;
	int sWapLength;
	std::vector<Chromosome> generationBestSolutions;
	std::vector<Chromosome> generationalOffsprings;//need change later
	std::map<int, int> demand;
	std::vector<std::vector<double>> distance;
	std::vector<int> customerCluster;
	Population ppl;
	Chromosome initialSolution;
	Chromosome incumbentSolution;
	Chromosome bestSolution;
public:
	Geneticalgorithm();
	Geneticalgorithm(const Geneticalgorithm & ga);	
	Geneticalgorithm(int populationSize, int numberOfNodes, int depotNode, int capacityLimit, int kChainLength, int sWapLength, std::map<int, int> demand, std::vector<std::vector<double>> distance, std::vector<int> customerCluster);
	Geneticalgorithm(int depotNode, int capacityLimit, std::map<int, int> demand, std::vector<std::vector<double>> distance, std::vector<int> customerCluster);
	void populateInitialGeneration();
	void showCurrentGenerationBestChromosome();
	std::vector<Chromosome> getGenerationBestSolutions();
	std::vector<Chromosome> getGenerationalOffsprings();
	Chromosome getGASolution();
	void runGeneticAlgorithm();
	void showGASolution();
	void showCurrentGeneration();
};



#endif // !GENETICALGOIRTHM_H



