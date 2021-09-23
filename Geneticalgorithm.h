#include <iostream>

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

class Geneticalgorithm{
	public:
		Geneticalgorithm();
		Geneticalgorithm(int, std::string);

		std::string show_name();
		int show_number();

	private:
		int number;
		std::string name;
};


#endif // !GENETICALGOIRTHM_H



