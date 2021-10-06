#include <iostream>
#include <vector>
#include <list>

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
//chromosomes are fixed size
class Chromosome {
private:
	std::list<int> chromosome;
	double fitness;
	int separatorInt;
	bool isFeasible;
	double reProductionProbability;
	double reInsertionProbabililty;
	double mutationProbability;
	int size;//size of the chromosome 
public:
	Chromosome();//default constructor
	Chromosome(const Chromosome & chrom);//copy constructor
	std::list<int> getChromosome();
	double getFitness();
	void showFitness();
	int getSize();
	void updateSize(int size);
	void updateToFeasibleChromosome();
	int getSeparatorInt();
	void showSeparatorInt();
	bool getFeasibilityStatus();
	void checkFeasiblity();
	double getReproductionProbability();
	double getReinsertionProbability();
	double getMutationProbability();
	void updateMutationProbability();
	void updateReproductionAndReinsertionProbabilities();
};



//class for elitism operator
class Elitism {
	Chromosome elit;
public:
	Elitism();//default constructor
	Elitism(const Elitism & elit);//copy constructor
};

//class for replication operator
class Replication {
	Chromosome replica;
public:
	Replication();//default constructor
	Replication(const Replication & replica);//copy constructor
};

//class for cross over operator
class CrossOver {
private:
	Chromosome parent1;
	Chromosome parent2;
	Chromosome offspring;
public:
	CrossOver();//default constructor
	CrossOver(const CrossOver & crossr);//copy constructor
	void performPertiallyMappedCrossover();
	void showParent1();
	void showParent2();
	void showOffspring();
	Chromosome getParent1();
	Chromosome getParent2();
	Chromosome getOffspring();
};

//class for mutation operator
class Mutation {
	Chromosome originalChrom;
	Chromosome mutatedChrom;
public:
	Mutation();//default constructor
	Mutation(const Mutation & mutatn);//copy constructor
	Chromosome getOriginalChromosome();
	Chromosome getMutatedChromosome();
	void showOriginalChromosome();
	void showMutatedChromosome();
	void performMutation();
};

//class for a of population
class Population {
	std::list<std::list<int>> population;
public:
	Population();//default constructor
	Population(const Population & poplatn);//copy constructor
	std::list<std::list<int>> getPopulation();
	void showPopulation();
};

//class for tournament
class Tournament {
private:
	Population population;
public:
	Tournament();//default constructor
	Tournament(const Tournament & tourmnt);//copy constructor
	void arrangeTournament();
};

//class for genetic algorithm
class Geneticalgorithm{
private:
	Population currentPopulation;
	Population futurePopulation;
	Elitism elitism;
	Replication replication;
	Mutation mutation;
	CrossOver crossover;
	Tournament tournament;
public:
	Geneticalgorithm();//default constructor
	Geneticalgorithm(const Geneticalgorithm & ga);//copy constructor	
	void selectElitChromosome();
	void arrangeTournament();
	void performPartiallyMappedCrossover();
	void performChainMutation();
	void showCurrentPopulation();
	void showFuturePopulation();
	Population getCurrentPopulation();
	Population getFuturePopulation();
	void runGeneticAlgorithm();
};


#endif // !GENETICALGOIRTHM_H



