#include "Geneticalgorithm.h"



//default constructor
Chromosome::Chromosome() {
	fitness    = 0.0;
	separatorInt  = 0;
	isFeasible    = false;
	reProductionProbability   = 0.0;
	reInsertionProbabililty   = 0.0;
	mutationProbability       = 0.0;
	size = 0;
}

//copy constructor
Chromosome::Chromosome(const Chromosome& chrom) {
	chromosome = chrom.chromosome;
	fitness = chrom.fitness;
	separatorInt = chrom.separatorInt;
	isFeasible = chrom.isFeasible;
	reProductionProbability = chrom.reProductionProbability;
	reInsertionProbabililty = chrom.reInsertionProbabililty;
	mutationProbability = chrom.mutationProbability;
	size = chrom.size;
}

//returns chromosome
std::list<int> Chromosome::getChromosome() {

}

//returns fitness value
double Chromosome::getFitness() {

}

//prints fitness value
void Chromosome::showFitness() {

}

//updates to feasible chromosome
void Chromosome::updateToFeasibleChromosome() {

}

//returns seperator integer value
int Chromosome::getSeparatorInt() {

}

//prints separator integer
void Chromosome::showSeparatorInt() {

}

//returns feasibility status
bool Chromosome::getFeasibilityStatus() {

}

//checks whether the chromosome is feasible
void Chromosome::checkFeasiblity() {

}

//returns reproduction probability
double Chromosome::getReproductionProbability() {

}

//returns reinsertion probability
double Chromosome::getReinsertionProbability() {

}

//returns mutation probability
double Chromosome::getMutationProbability() {

}

//updates reproduction and reinsertion probabilities
void Chromosome::updateReproductionAndReinsertionProbabilities() {

}

//updates mutation probability
void Chromosome::updateMutationProbability() {

}

//returns size of the chromosome
int Chromosome::getSize() {

}

//updates size of the chromosome
void Chromosome::updateSize(int size) {

}

//default constructor
Elitism::Elitism() {

}
	
//copy constructor
Elitism::Elitism(const Elitism & elit) {

}

//default constructor
Replication::Replication() {

}

//copy constructor
Replication::Replication(const Replication & replica) {

}

//default constructor
CrossOver::CrossOver() {

}

//copy constructor
CrossOver::CrossOver(const CrossOver & crossr) {

}

//performs partially mapped crossover operation
void CrossOver::performPertiallyMappedCrossover() {

}

//prints parent 1
void CrossOver::showParent1() {

}

//prints parent 2
void CrossOver::showParent2() {

}

//prints offspring
void CrossOver::showOffspring() {

}

//returns parent 1
Chromosome CrossOver::getParent1() {

}

//returns parent 2
Chromosome CrossOver::getParent2() {

}

//returns offspring
Chromosome CrossOver::getOffspring() {

}

//default constructor
Mutation::Mutation() {

}
	
//copy constructor
Mutation::Mutation(const Mutation & mutatn){

}

//returns premutation chromosome
Chromosome Mutation::getOriginalChromosome() {

}

//returns mutated chromosome
Chromosome Mutation::getMutatedChromosome() {

}

//prints premutation chromosome
void Mutation::showOriginalChromosome() {

}

//prints mutated chromosome
void Mutation::showMutatedChromosome() {

}

//performs chain mutation
void Mutation::performMutation() {

}


//default constructor
Population::Population() {

}

//copy constructor
Population::Population(const Population & gentn) {

}

//default constructor
Tournament::Tournament() {

}
	
//copy constructor
Tournament::Tournament(const Tournament & tourmnt){

}

//arrange tournament
void Tournament::arrangeTournament() {

}

//default constructor
Geneticalgorithm::Geneticalgorithm() {

}

//copy constructor	
Geneticalgorithm::Geneticalgorithm(const Geneticalgorithm & ga) {

}

//selects elit chromosome from current generation
void Geneticalgorithm::selectElitChromosome() {

}

//arranges tournament with current generation population
void Geneticalgorithm::arrangeTournament() {

}

//performs partially mapped crossover
void Geneticalgorithm::performPartiallyMappedCrossover() {

}

//performs mutation operation
void Geneticalgorithm::performChainMutation() {

}

//prints current generation population
void Geneticalgorithm::showCurrentPopulation() {

}

//prints future generation population
void Geneticalgorithm::showFuturePopulation() {

}

//returns current generation population
Population Geneticalgorithm::getCurrentPopulation() {
	return currentPopulation;
}

//returns future generation population
Population Geneticalgorithm::getFuturePopulation() {
	return futurePopulation;
}

//runs genetic algorithm 
void Geneticalgorithm::runGeneticAlgorithm() {

}







